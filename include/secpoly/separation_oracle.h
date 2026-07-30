// Separation from optimization, by away-step Frank-Wolfe. See DESIGN.md 4.
//
// The reduction: the point of Sigma(A) nearest to y gives a separating direction
// gamma = y - x*, with margin exactly delta^2. Frank-Wolfe computes that
// projection using ONLY the linear optimization oracle -- its linear subproblem
// IS the oracle -- and the separation certificate appears long before the
// projection converges.

#ifndef SECPOLY_SEPARATION_ORACLE_H
#define SECPOLY_SEPARATION_ORACLE_H

#include <secpoly/affine_hull.h>
#include <secpoly/bootstrap.h>
#include <secpoly/opt_oracle.h>
#include <secpoly/rationalize.h>
#include <secpoly/types.h>

#include <deque>
#include <memory>

namespace secpoly {

class SeparationOracle {
public:
    SeparationOracle(const SecondaryOptOracle* oracle,
                     AffineHull hull,
                     BootstrapResult boot,
                     Tolerances tol = {})
        : oracle_(oracle), hull_(std::move(hull)), boot_(std::move(boot)),
          tol_(tol) {
        hull_.set_origin(boot_.x0);
        v_norm_ = oracle_->normalized_volume();

        // Seed the atom cache with the bootstrap vertices -- they are already
        // paid-for oracle answers.
        for (const auto& g : boot_.vertices) add_atom(g);
    }

    int n() const { return hull_.n(); }
    const AffineHull& hull() const { return hull_; }
    NT rho_simp() const { return boot_.rho_simp; }

    /// Decide whether the intrinsic point y lies in Q = {y : x0 + Uy in Sigma(A)}.
    SeparationResult separate(const VecX& y_intrinsic) const {
        SeparationResult res;
        const VecX x = hull_.to_ambient(y_intrinsic);

        // --- layer 0: ambient box, free ---------------------------------------
        // Every GKZ coordinate lies in [0, V_norm] (draft Lemma 6), so a point
        // outside that box is outside Sigma(A) with no oracle call at all.
        for (Eigen::Index i = 0; i < x.size(); ++i) {
            if (x(i) < -kBoxTol || x(i) > static_cast<NT>(v_norm_) + kBoxTol) {
                res.verdict = Verdict::Outside;
                res.cut_c.assign(static_cast<std::size_t>(x.size()), 0);
                if (x(i) < 0) {
                    res.cut_c[static_cast<std::size_t>(i)] = -1;
                    res.cut_h = 0;                 // max of -x_i over Sigma(A) is 0
                } else {
                    res.cut_c[static_cast<std::size_t>(i)] = 1;
                    res.cut_h = v_norm_;           // max of x_i is V_norm
                }
                return res;
            }
        }

        // --- layer 1: cached halfspaces, free ---------------------------------
        // Every cached (c, h) satisfies c^T z <= h for ALL z in Sigma(A), exactly.
        // Uses a POSITIVE tolerance so a near-boundary point is never declared
        // outside by cache alone; it falls through to Frank-Wolfe instead.
        for (const auto& hs : cuts_) {
            if (mixed_dot(hs.c, x) > static_cast<NT>(hs.h) + kCutTol) {
                res.verdict = Verdict::Outside;
                res.cut_c = hs.c;
                res.cut_h = hs.h;
                return res;
            }
        }

        // --- layer 2: Frank-Wolfe ---------------------------------------------
        return frank_wolfe(y_intrinsic, x, res);
    }

    std::size_t num_atoms() const { return atoms_.size(); }
    std::size_t num_cuts() const { return cuts_.size(); }
    std::size_t oracle_calls() const { return oracle_->num_calls(); }
    std::size_t lazy_hits() const { return lazy_hits_; }

private:
    struct Halfspace { IVec c; long long h; };

    static constexpr NT kBoxTol = 1e-9;
    static constexpr NT kCutTol = 1e-9;

    void add_atom(const IVec& g) const {
        for (const auto& a : atoms_exact_)
            if (a == g) return;
        atoms_exact_.push_back(g);
        VecX x(static_cast<Eigen::Index>(g.size()));
        for (std::size_t i = 0; i < g.size(); ++i)
            x(static_cast<Eigen::Index>(i)) = static_cast<NT>(g[i]);
        atoms_.push_back(hull_.to_intrinsic(x));
    }

    /// The loop of DESIGN.md 4.7.
    SeparationResult frank_wolfe(const VecX& y, const VecX& /*x_ambient*/,
                                 SeparationResult res) const {
        // Deterministic start: the bootstrap relative-interior point, i.e. the
        // chart origin. Different starts give different -- equally valid --
        // hyperplanes and call counts, so pinning it makes traces reproducible.
        VecX cur = VecX::Zero(hull_.n());

        // Active set as (atom index, weight); cur == sum w_i * atoms_[i].
        std::vector<std::pair<int, NT>> active;
        {
            // start at the barycentre of the bootstrap simplex == chart origin
            const NT w = NT(1) / static_cast<NT>(atoms_.size());
            for (std::size_t i = 0; i < atoms_.size(); ++i)
                active.emplace_back(static_cast<int>(i), w);
            cur = VecX::Zero(hull_.n());
            for (const auto& [idx, wt] : active) cur += wt * atoms_[idx];
        }

        for (int iter = 0; iter < tol_.max_fw_iters; ++iter) {
            res.fw_iterations = iter + 1;

            const VecX gamma = y - cur;                       // steepest direction
            const NT dist = gamma.norm();

            if (dist <= tol_.membership_tol) {
                res.verdict = Verdict::Inside;
                res.witness = cur;
                res.distance = dist;
                return res;
            }

            const VecX gamma_ambient = hull_.U() * gamma;
            const IVec c = rationalize(gamma_ambient, hull_);
            if (is_zero(c)) {
                res.verdict = Verdict::Ambiguous;
                res.distance = dist;
                return res;
            }

            // --- lazy step: try the cached atoms before paying for respol -----
            // Cached atoms are genuine vertices of Sigma(A), so the max of gamma
            // over them is a LOWER bound on the true max. If that lower bound
            // already reaches gamma^T y, no hyperplane in this direction can
            // separate, and we may take the Frank-Wolfe step toward the cached
            // atom without an oracle call. This is the lazified conditional
            // gradient pattern (Braun-Pokutta-Zink); the caching is where the
            // speedup comes from.
            //
            // It never affects correctness: an OUTSIDE verdict still requires a
            // real oracle answer, since only that establishes the TRUE maximum.
            {
                const NT gy = gamma.dot(y);
                int best = -1;
                NT best_val = -std::numeric_limits<NT>::infinity();
                for (std::size_t i = 0; i < atoms_.size(); ++i) {
                    const NT v = gamma.dot(atoms_[i]);
                    if (v > best_val) { best_val = v; best = static_cast<int>(i); }
                }
                if (best >= 0 && best_val >= gy) {
                    ++lazy_hits_;
                    if (!fw_step(best, gamma, cur, active)) {
                        res.verdict = (dist <= tol_.membership_tol * 100)
                                          ? Verdict::Inside : Verdict::Ambiguous;
                        res.witness = cur;
                        res.distance = dist;
                        return res;
                    }
                    continue;   // no oracle call was needed
                }
            }

            IVec gkz, c_used;
            if (!oracle_->maximize_robust(c, gkz, c_used)) {
                res.verdict = Verdict::Ambiguous;   // degeneracy retries exhausted
                res.distance = dist;
                return res;
            }
            ++res.oracle_calls;

            // --- the separation test ------------------------------------------
            // Test with c_used -- the vector the oracle ACTUALLY maximized -- not
            // with gamma. Using gamma would silently reintroduce a rounding-sized
            // tolerance into what is an unconditional proof. Both sides here are
            // exact integer arithmetic.
            const long long cs = exact_dot(c_used, gkz);
            const IVec y_amb_int = ambient_int(y);
            const bool exact_y = !y_amb_int.empty();

            bool separated;
            if (exact_y) {
                separated = cs < exact_dot(c_used, y_amb_int);
            } else {
                separated = static_cast<NT>(cs) <
                            mixed_dot(c_used, hull_.to_ambient(y)) - kCutTol;
            }

            // Every oracle answer yields an exactly valid halfspace, whatever
            // c_used turned out to be. Cache it unconditionally.
            remember_cut(c_used, cs);
            add_atom(gkz);

            if (separated) {
                res.verdict = Verdict::Outside;
                res.cut_c = c_used;
                res.cut_h = cs;
                return res;
            }

            // --- Frank-Wolfe step ---------------------------------------------
            if (!fw_step(atom_index(gkz), gamma, cur, active)) {
                // Stall: the rounded direction disagrees with the true steepest
                // direction by more than the remaining gap, which can only happen
                // at numerical convergence.
                res.verdict = (dist <= tol_.membership_tol * 100)
                                  ? Verdict::Inside : Verdict::Ambiguous;
                res.witness = cur;
                res.distance = dist;
                return res;
            }
        }

        res.verdict = Verdict::Ambiguous;
        res.distance = (y - cur).norm();
        res.witness = cur;
        return res;
    }

    /// One AWAY-STEP Frank-Wolfe iteration. Returns false on stall.
    ///
    /// Plain forward-only Frank-Wolfe has an O(1/k) rate, so separating a point at
    /// distance delta from the body needs O(D^2/delta^2) iterations. That is
    /// hopeless near the boundary -- at D ~ 10, delta ~ 0.03 it is ~1e5 iterations
    /// -- and the practical consequence is that near-boundary EXTERIOR points come
    /// back Ambiguous and get accepted, putting samples measurably outside
    /// Sigma(A). Away steps restore linear convergence, O(log(D/delta)), which is
    /// what makes the budget sufficient (Lacoste-Julien & Jaggi 2015).
    ///
    /// The idea: as well as moving TOWARD the best atom, allow moving AWAY from
    /// the worst atom currently carrying weight. Forward-only FW zigzags when the
    /// optimum lies on a face; away steps let it back off an atom it
    /// over-committed to, straightening the path.
    bool fw_step(int s_idx, const VecX& gamma, VecX& cur,
                 std::vector<std::pair<int, NT>>& active) const {
        const VecX& s = atoms_[static_cast<std::size_t>(s_idx)];

        const VecX d_fw = s - cur;
        const NT gain_fw = gamma.dot(d_fw);          // == the Frank-Wolfe gap, >= 0

        // Worst active atom: the one least aligned with gamma.
        int away_pos = -1;
        NT worst = std::numeric_limits<NT>::infinity();
        for (std::size_t i = 0; i < active.size(); ++i) {
            const NT v = gamma.dot(atoms_[static_cast<std::size_t>(active[i].first)]);
            if (v < worst) { worst = v; away_pos = static_cast<int>(i); }
        }

        bool use_away = false;
        NT alpha_max = 1;
        VecX dir = d_fw;

        if (away_pos >= 0 && active.size() > 1) {
            const auto& av = active[static_cast<std::size_t>(away_pos)];
            const VecX d_away = cur - atoms_[static_cast<std::size_t>(av.first)];
            const NT gain_away = gamma.dot(d_away);
            if (gain_away > gain_fw) {
                use_away = true;
                dir = d_away;
                // Cannot push the away atom's weight below zero.
                alpha_max = (av.second < 1 - 1e-12)
                                ? av.second / (1 - av.second)
                                : NT(1e12);
            }
        }

        const NT dd = dir.squaredNorm();
        if (dd < 1e-300) return false;

        NT alpha = gamma.dot(dir) / dd;              // exact line search
        if (alpha < 0) alpha = 0;
        if (alpha > alpha_max) alpha = alpha_max;
        if (alpha <= 0) return false;

        if (!use_away) {
            for (auto& kv : active) kv.second *= (1 - alpha);
            bool found = false;
            for (auto& kv : active)
                if (kv.first == s_idx) { kv.second += alpha; found = true; break; }
            if (!found) active.emplace_back(s_idx, alpha);
        } else {
            for (auto& kv : active) kv.second *= (1 + alpha);
            active[static_cast<std::size_t>(away_pos)].second -= alpha;
        }

        // Drop atoms that have fallen to zero weight (a "drop step").
        active.erase(std::remove_if(active.begin(), active.end(),
                                    [](const std::pair<int, NT>& p) {
                                        return p.second < 1e-14;
                                    }),
                     active.end());

        cur += alpha * dir;
        return true;
    }

    /// Exact integer ambient image of an intrinsic point, when one exists.
    /// Returns empty when the point is not integral, which is the normal case for
    /// a sampled point -- the caller then falls back to the tolerance-guarded
    /// floating comparison.
    IVec ambient_int(const VecX& y) const {
        const VecX x = hull_.to_ambient(y);
        IVec out(static_cast<std::size_t>(x.size()));
        for (Eigen::Index i = 0; i < x.size(); ++i) {
            const NT r = std::round(x(i));
            if (std::fabs(x(i) - r) > 1e-7) return {};
            out[static_cast<std::size_t>(i)] = static_cast<long long>(r);
        }
        return out;
    }

    int atom_index(const IVec& g) const {
        for (std::size_t i = 0; i < atoms_exact_.size(); ++i)
            if (atoms_exact_[i] == g) return static_cast<int>(i);
        add_atom(g);
        return static_cast<int>(atoms_exact_.size() - 1);
    }

    void remember_cut(const IVec& c, long long h) const {
        // Purging is provably SAFE: dropping a cut only loses a fast path, since
        // cuts are used exclusively to declare OUTSIDE and every one of them is
        // exactly valid. It can cost oracle calls, never correctness.
        if (cuts_.size() >= kMaxCuts) cuts_.pop_front();
        cuts_.push_back({c, h});
    }

    static constexpr std::size_t kMaxCuts = 4096;

    const SecondaryOptOracle* oracle_;
    AffineHull hull_;
    BootstrapResult boot_;
    Tolerances tol_;
    long long v_norm_ = 0;

    // mutable: is_in() must be const for volesti's BallWalk (apply takes
    // `BallPolytope const&`), but every query legitimately grows these caches.
    mutable std::vector<IVec> atoms_exact_;
    mutable std::vector<VecX> atoms_;
    mutable std::deque<Halfspace> cuts_;
    mutable std::size_t lazy_hits_ = 0;
};

}  // namespace secpoly

#endif  // SECPOLY_SEPARATION_ORACLE_H
