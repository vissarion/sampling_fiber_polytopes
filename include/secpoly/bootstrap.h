// Draft Lemma 4: a CERTIFIED relative-interior point, obtained with oracle calls
// only and no genericity assumption.
//
// Yields four things the sampler needs and cannot otherwise get:
//   * n+1 affinely independent GKZ vertices,
//   * x0 = their barycentre, provably in relint Sigma(A),
//   * the initial atom cache (free -- they are oracle answers already),
//   * rho_simp, a certified inner-ball radius (draft Remark 2).

#ifndef SECPOLY_BOOTSTRAP_H
#define SECPOLY_BOOTSTRAP_H

#include <secpoly/affine_hull.h>
#include <secpoly/opt_oracle.h>
#include <secpoly/rationalize.h>
#include <secpoly/types.h>

#include <random>
#include <stdexcept>

namespace secpoly {

struct BootstrapResult {
    std::vector<IVec> vertices;   ///< n+1 affinely independent GKZ vectors (ambient)
    VecX x0;                      ///< barycentre, in relint Sigma(A)
    NT   rho_simp = 0;            ///< certified inner radius about x0, in the chart
    int  oracle_calls = 0;
};

/// Grow a subspace E from {0} by probing directions orthogonal to it.
///
/// Lemma 4's argument: for u in L0 \ {0}, u^T x is non-constant on Sigma(A), so its
/// max and min differ; since every point of v0 + E has the same u-value, at least
/// one of argmax(+u), argmax(-u) lies off v0 + E and raises dim E by one. That
/// gives <= 1 + 2n calls deterministically.
///
/// Perturbation (needed when a lifting is degenerate) breaks the exact
/// orthogonality that argument relies on, so we fall back to redrawing u at random
/// from L0 ∩ E^perp -- which works with probability 1. That makes this Las Vegas
/// rather than deterministic, which is the honest cost of degeneracy-robustness.
inline BootstrapResult bootstrap(const SecondaryOptOracle& oracle,
                                 const AffineHull& hull,
                                 std::uint64_t seed = 20260728u) {
    const int n = hull.n();
    const int m = hull.m();
    BootstrapResult out;

    std::mt19937_64 rng(seed);
    std::normal_distribution<NT> gauss(0.0, 1.0);

    auto call = [&](const VecX& dir_intrinsic, IVec& gkz) -> bool {
        const VecX c_ambient = hull.U() * dir_intrinsic;
        const IVec c = rationalize(c_ambient, hull);
        if (is_zero(c)) return false;
        IVec used;
        const bool ok = oracle.maximize_robust(c, gkz, used);
        ++out.oracle_calls;
        return ok;
    };

    // Intrinsic images of the vertices found so far, relative to vertices[0].
    // Their span is E; we stop when rank(E) == n.
    // NOTE the explicit `-> VecX`. Without it `auto` deduces an Eigen EXPRESSION
    // TEMPLATE (Eigen::Product<...>) that stores references to its operands --
    // including the local `x` -- and every later use reads freed stack memory.
    // The symptom is a bad_alloc at an unrelated site, which is a miserable bug
    // to chase. Always name the return type when returning an Eigen expression.
    auto intrinsic_of = [&](const IVec& g) -> VecX {
        VecX x(m);
        for (int i = 0; i < m; ++i) x(i) = static_cast<NT>(g[i]);
        return hull.U().transpose() * x;   // origin-independent: differences only
    };

    // --- first vertex ---------------------------------------------------------
    {
        VecX u(n);
        for (int i = 0; i < n; ++i) u(i) = gauss(rng);
        IVec g;
        if (!call(u, g)) throw std::runtime_error("bootstrap: first oracle call failed");
        out.vertices.push_back(g);
    }
    const VecX base = intrinsic_of(out.vertices[0]);

    MatX diffs(n, 0);   // columns: intrinsic(v_i) - intrinsic(v_0)
    auto current_rank = [&]() -> int {
        if (diffs.cols() == 0) return 0;
        Eigen::JacobiSVD<MatX> svd(diffs);
        const auto& sv = svd.singularValues();
        if (sv.size() == 0 || sv(0) <= 0) return 0;
        int r = 0;
        for (int i = 0; i < sv.size(); ++i)
            if (sv(i) > sv(0) * 1e-10) ++r;   // RELATIVE threshold, never absolute
        return r;
    };

    // A direction in E^perp: the last right-singular vectors of diffs^T, or a
    // random direction when E is trivial.
    auto orthogonal_direction = [&]() -> VecX {
        VecX u(n);
        for (int i = 0; i < n; ++i) u(i) = gauss(rng);
        const int r = current_rank();
        if (diffs.cols() == 0 || r == 0) return u.normalized();
        // remove the component lying inside E
        Eigen::HouseholderQR<MatX> qr(diffs);
        const MatX Q = qr.householderQ();          // n x n
        const MatX Ein = Q.leftCols(r);            // n x r, orthonormal basis of E
        u -= Ein * (Ein.transpose() * u);
        const NT nrm = u.norm();
        if (nrm < 1e-12) return VecX::Zero(n);
        return u / nrm;
    };

    int guard = 0;
    while (current_rank() < n) {
        if (++guard > 50 * (n + 1))
            throw std::runtime_error(
                "bootstrap: could not reach n affinely independent vertices");

        const VecX u = orthogonal_direction();
        if (u.norm() < 1e-12) continue;

        const int rank_before = current_rank();
        bool progressed = false;

        for (int sign = 0; sign < 2 && !progressed; ++sign) {
            IVec g;
            if (!call(sign == 0 ? VecX(u) : VecX(-u), g)) continue;

            // Does this vertex enlarge span{v_i - v_0}? Try it, keep it only if so.
            //
            // Built by explicit block assignment, NOT Eigen's comma initializer:
            // `trial << saved, vec` misbehaves when `saved` has zero columns,
            // which is exactly the state on the first iteration.
            const MatX saved = diffs;
            MatX trial(n, saved.cols() + 1);
            if (saved.cols() > 0) trial.leftCols(saved.cols()) = saved;
            trial.col(saved.cols()) = intrinsic_of(g) - base;
            diffs = trial;

            if (current_rank() > rank_before) {
                out.vertices.push_back(g);
                progressed = true;
            } else {
                diffs = saved;
            }
        }
        // If neither sign helped, the loop simply redraws a fresh random u.
    }

    if (static_cast<int>(out.vertices.size()) != n + 1)
        throw std::runtime_error("bootstrap: wrong number of vertices");

    // --- x0 = barycentre, certified in relint Sigma(A) ------------------------
    out.x0 = VecX::Zero(m);
    for (const auto& g : out.vertices)
        for (int i = 0; i < m; ++i) out.x0(i) += static_cast<NT>(g[i]);
    out.x0 /= static_cast<NT>(n + 1);

    // --- rho_simp (draft Remark 2) -------------------------------------------
    // For the simplex S = conv(v_0..v_n) with barycentre x0,
    //     B(x0, (1/(n+1)) * min_i h_i) is contained in S,
    // where h_i is the altitude from v_i to the opposite facet, measured inside
    // the affine hull. Unconditional, though it can be very small.
    {
        std::vector<VecX> V;
        V.reserve(n + 1);
        for (const auto& g : out.vertices) {
            VecX x(m);
            for (int i = 0; i < m; ++i) x(i) = static_cast<NT>(g[i]);
            V.push_back(hull.U().transpose() * x);   // intrinsic, R^n
        }

        NT min_h = std::numeric_limits<NT>::infinity();
        for (int i = 0; i <= n; ++i) {
            // Opposite facet spanned by {V[j] - V[ref] : j != i, j != ref}.
            const int ref = (i == 0) ? 1 : 0;
            MatX W(n, n - 1);
            int col = 0;
            for (int j = 0; j <= n; ++j) {
                if (j == i || j == ref) continue;
                W.col(col++) = V[j] - V[ref];
            }
            // Normal to that facet = the direction orthogonal to all of W.
            VecX normal;
            if (n == 1) {
                normal = VecX::Ones(1);
            } else {
                Eigen::JacobiSVD<MatX> svd(W.transpose(),
                                           Eigen::ComputeFullV);
                normal = svd.matrixV().col(n - 1);   // smallest singular direction
            }
            const NT nn = normal.norm();
            if (nn < 1e-14) { min_h = 0; break; }
            const NT h = std::fabs(normal.dot(V[i] - V[ref])) / nn;
            min_h = std::min(min_h, h);
        }
        out.rho_simp = (min_h > 0 && std::isfinite(min_h))
                           ? min_h / static_cast<NT>(n + 1)
                           : NT(0);
    }

    return out;
}

}  // namespace secpoly

#endif  // SECPOLY_BOOTSTRAP_H
