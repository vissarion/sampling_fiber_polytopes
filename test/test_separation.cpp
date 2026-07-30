// Validation of the Frank-Wolfe separation oracle against exact ground truth.
//
// For small instances respol lets us enumerate ALL vertices of Sigma(A), so we
// can decide membership independently (by exact convex-hull membership over the
// complete vertex list) and check every verdict.

#include <secpoly/bootstrap.h>
#include <secpoly/opt_oracle.h>
#include <secpoly/separation_oracle.h>

#include <cstdio>
#include <random>
#include <set>
#include <vector>

static int failures = 0;

#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, (msg));    \
            ++failures;                                                    \
        }                                                                  \
    } while (0)

using namespace secpoly;

static std::set<IVec> enumerate_vertices(const SecondaryOptOracle& o, int trials) {
    std::set<IVec> found;
    std::mt19937_64 rng(4242);
    std::uniform_int_distribution<long long> dist(-2000, 2000);
    IVec c(o.ambient_dim()), gkz, used;
    for (int t = 0; t < trials; ++t) {
        for (auto& x : c) x = dist(rng);
        if (o.maximize_robust(c, gkz, used)) found.insert(gkz);
    }
    return found;
}

int main() {
    // Example 1: A = {1,2,4,6}, n = 2. All four vertices are known exactly.
    const std::vector<std::vector<long long>> A{{1}, {2}, {4}, {6}};
    SecondaryOptOracle oracle(1, A);
    AffineHull hull(1, A);

    const BootstrapResult boot = bootstrap(oracle, hull);
    CHECK(static_cast<int>(boot.vertices.size()) == hull.n() + 1,
          "bootstrap must return n+1 affinely independent vertices");
    CHECK(boot.rho_simp > 0, "rho_simp must be strictly positive");
    std::printf("bootstrap: %zu vertices, rho_simp=%.6f, %d oracle calls\n",
                boot.vertices.size(), boot.rho_simp, boot.oracle_calls);

    auto sep = std::make_shared<SeparationOracle>(&oracle, hull, boot);

    // --- ground truth: the complete vertex set --------------------------------
    const std::set<IVec> verts = enumerate_vertices(oracle, 300);
    CHECK(verts.size() == 4, "Example 1 has exactly 4 vertices");

    AffineHull h2(1, A);
    h2.set_origin(boot.x0);
    std::vector<VecX> V;
    for (const auto& g : verts) {
        VecX x(4);
        for (int i = 0; i < 4; ++i) x(i) = static_cast<NT>(g[i]);
        V.push_back(h2.to_intrinsic(x));
    }

    // Exact membership in conv(V) for n=2, by checking the point against every
    // edge of the (convex) polygon -- independent of the FW machinery.
    auto in_hull = [&](const VecX& y) {
        // order the 4 planar vertices by angle about their centroid
        VecX ctr = VecX::Zero(2);
        for (const auto& v : V) ctr += v;
        ctr /= static_cast<NT>(V.size());
        std::vector<VecX> P = V;
        std::sort(P.begin(), P.end(), [&](const VecX& a, const VecX& b) {
            return std::atan2(a(1) - ctr(1), a(0) - ctr(0)) <
                   std::atan2(b(1) - ctr(1), b(0) - ctr(0));
        });
        for (std::size_t i = 0; i < P.size(); ++i) {
            const VecX& p = P[i];
            const VecX& q = P[(i + 1) % P.size()];
            const NT cross = (q(0) - p(0)) * (y(1) - p(1)) -
                             (q(1) - p(1)) * (y(0) - p(0));
            if (cross < -1e-9) return false;   // right of an edge => outside
        }
        return true;
    };

    // --- soundness: OUTSIDE must never fire for an interior point -------------
    std::mt19937_64 rng(2026);
    std::uniform_real_distribution<NT> u(0, 1);
    int inside_tested = 0, sound = 0;
    for (int t = 0; t < 400; ++t) {
        // random convex combination => provably inside
        NT w[4];
        NT s = 0;
        for (int i = 0; i < 4; ++i) { w[i] = u(rng); s += w[i]; }
        VecX y = VecX::Zero(2);
        for (int i = 0; i < 4; ++i) y += (w[i] / s) * V[static_cast<std::size_t>(i)];

        const SeparationResult r = sep->separate(y);
        ++inside_tested;
        if (r.verdict != Verdict::Outside) ++sound;
    }
    CHECK(sound == inside_tested,
          "SOUNDNESS: the OUTSIDE branch must be unreachable for interior points");
    std::printf("soundness: %d/%d interior points not declared outside\n",
                sound, inside_tested);

    // --- every OUTSIDE certificate must be a genuine proof --------------------
    int outside_tested = 0, valid_cuts = 0;
    std::uniform_real_distribution<NT> wide(-40, 40);
    for (int t = 0; t < 400; ++t) {
        VecX y(2);
        y << wide(rng), wide(rng);
        if (in_hull(y)) continue;

        const SeparationResult r = sep->separate(y);
        if (r.verdict != Verdict::Outside) continue;
        ++outside_tested;

        // The cut must dominate EVERY vertex of Sigma(A) and be violated by y.
        bool ok = true;
        for (const auto& g : verts)
            if (exact_dot(r.cut_c, g) > r.cut_h) ok = false;
        const VecX x = h2.to_ambient(y);
        if (mixed_dot(r.cut_c, x) <= static_cast<NT>(r.cut_h)) ok = false;
        if (ok) ++valid_cuts;
    }
    CHECK(outside_tested > 0, "should have tested some exterior points");
    CHECK(valid_cuts == outside_tested,
          "every OUTSIDE cut must dominate all vertices AND be violated by y");
    std::printf("certificates: %d/%d cuts are genuine proofs\n",
                valid_cuts, outside_tested);

    // --- agreement with ground truth on a mixed sample ------------------------
    int agree = 0, decided = 0;
    for (int t = 0; t < 400; ++t) {
        VecX y(2);
        y << wide(rng) * 0.35, wide(rng) * 0.35;   // straddles the boundary
        const bool truth = in_hull(y);
        const SeparationResult r = sep->separate(y);
        if (r.verdict == Verdict::Ambiguous) continue;
        ++decided;
        if ((r.verdict == Verdict::Inside) == truth) ++agree;
    }
    CHECK(decided > 0, "should have decided some queries");
    CHECK(agree == decided, "verdicts must agree with exact hull membership");
    std::printf("agreement: %d/%d verdicts match exact convex-hull membership\n",
                agree, decided);

    std::printf("atoms=%zu cuts=%zu total_oracle_calls=%zu\n",
                sep->num_atoms(), sep->num_cuts(), sep->oracle_calls());

    if (failures == 0) std::printf("test_separation: PASS\n");
    return failures == 0 ? 0 : 1;
}
