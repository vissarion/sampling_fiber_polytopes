// End-to-end: volesti's BallWalk sampling Sigma(A) through the oracle.
//
// Ground truth is available because respol can enumerate ALL vertices for small
// instances. We build a completely independent uniform sampler by rejection
// inside the bounding box, deciding membership against the complete vertex list,
// and compare distributions. That reference shares no code with the ball walk or
// with Frank-Wolfe.

#include <secpoly/bootstrap.h>
#include <secpoly/opt_oracle.h>
#include <secpoly/secondary_body.h>
#include <secpoly/separation_oracle.h>

#include "cartesian_geom/cartesian_kernel.h"
#include "generators/boost_random_number_generator.hpp"
#include "sampling/sphere.hpp"                  // GetPointInDsphere, used by BallWalk
#include "random_walks/uniform_ball_walk.hpp"
#include "volume/sampling_policies.hpp"         // PushBackWalkPolicy
#include "sampling/random_point_generators.hpp" // RandomPointGenerator
#include "sampling/sampling.hpp"                // uniform_sampling

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <list>
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
using Kernel = Cartesian<double>;
using Point = typename Kernel::Point;
using RNGType = BoostRandomNumberGenerator<boost::mt19937, double, 42>;

int main() {
    // Example 1: A = {1,2,4,6}. Sigma(A) is a quadrilateral -- 2-dimensional,
    // living in R^4.
    const std::vector<std::vector<long long>> A{{1}, {2}, {4}, {6}};
    SecondaryOptOracle oracle(1, A);
    AffineHull hull(1, A);
    const BootstrapResult boot = bootstrap(oracle, hull);

    auto sep = std::make_shared<SeparationOracle>(&oracle, hull, boot);
    SecondaryBody<Point> body(sep, boot.rho_simp);

    // THE critical assertion: volesti must be told n, never m. Sigma(A) has
    // measure zero in R^4, so a ball walk there would accept nothing.
    CHECK(body.dimension() == 2u,
          "body.dimension() must be n = m-d-1 = 2, NOT the ambient m = 4");

    // --- ground truth: all vertices, then exact planar hull membership --------
    std::set<IVec> verts;
    {
        std::mt19937_64 rng(31337);
        std::uniform_int_distribution<long long> di(-2000, 2000);
        IVec c(oracle.ambient_dim()), g, u;
        for (int t = 0; t < 400; ++t) {
            for (auto& x : c) x = di(rng);
            if (oracle.maximize_robust(c, g, u)) verts.insert(g);
        }
    }
    CHECK(verts.size() == 4, "Example 1 has 4 vertices");

    AffineHull h2(1, A);
    h2.set_origin(boot.x0);
    std::vector<VecX> V;
    for (const auto& g : verts) {
        VecX x(4);
        for (int i = 0; i < 4; ++i) x(i) = static_cast<NT>(g[i]);
        V.push_back(h2.to_intrinsic(x));
    }
    VecX ctr = VecX::Zero(2);
    for (const auto& v : V) ctr += v;
    ctr /= static_cast<NT>(V.size());
    std::sort(V.begin(), V.end(), [&](const VecX& a, const VecX& b) {
        return std::atan2(a(1) - ctr(1), a(0) - ctr(0)) <
               std::atan2(b(1) - ctr(1), b(0) - ctr(0));
    });
    auto in_hull = [&](const VecX& y) {
        for (std::size_t i = 0; i < V.size(); ++i) {
            const VecX& p = V[i];
            const VecX& q = V[(i + 1) % V.size()];
            if ((q(0) - p(0)) * (y(1) - p(1)) - (q(1) - p(1)) * (y(0) - p(0)) < 0)
                return false;
        }
        return true;
    };

    // --- reference sampler: rejection, independent of everything above --------
    VecX lo = V[0], hi = V[0];
    for (const auto& v : V) { lo = lo.cwiseMin(v); hi = hi.cwiseMax(v); }
    std::vector<VecX> ref;
    {
        std::mt19937_64 rng(2718);
        std::uniform_real_distribution<NT> ux(lo(0), hi(0)), uy(lo(1), hi(1));
        while (ref.size() < 20000) {
            VecX p(2);
            p << ux(rng), uy(rng);
            if (in_hull(p)) ref.push_back(p);
        }
    }

    // --- the ball walk --------------------------------------------------------
    const unsigned int walk_len = 5, num = 20000, burns = 500;
    RNGType rng(body.dimension());
    Point start(body.dimension());          // zero == x0, certified in relint
    std::list<Point> pts;

    // Explicit step size: rho_simp is a simplex inradius and badly
    // under-estimates the body's, so the derived 4*rho/sqrt(n) would be far too
    // small. Tuned below by acceptance rate.
    const double L = 2.5;
    BallWalk walk(L);
    body.reset_counters();
    uniform_sampling(pts, body, rng, walk, walk_len, num, start, burns);

    std::printf("ball walk: %zu samples, acceptance=%.3f, oracle_calls=%zu, "
                "lazy_hits=%zu, atoms=%zu\n",
                pts.size(), body.acceptance_rate(), sep->oracle_calls(),
                sep->lazy_hits(), sep->num_atoms());

    CHECK(pts.size() == num, "should have produced the requested sample count");
    CHECK(body.acceptance_rate() > 0.02,
          "acceptance rate must not collapse -- a near-zero rate is the "
          "signature of sampling in the wrong dimension");

    const NT span = (hi - lo).norm();   // the body's own length scale

    // --- containment: measure DEPTH of violation, not a raw count -------------
    // A sample sitting 1e-12 outside is not a defect, it is the GLS weak
    // membership slack the whole method operates under. What matters is that no
    // sample is MEANINGFULLY outside, so we measure how far outside the worst one
    // is, relative to the body's own scale.
    auto depth_outside = [&](const VecX& y) {
        NT worst = 0;   // >0 means outside, in units of length
        for (std::size_t i = 0; i < V.size(); ++i) {
            const VecX& p = V[i];
            const VecX& q = V[(i + 1) % V.size()];
            const VecX e = q - p;
            const NT len = e.norm();
            if (len < 1e-300) continue;
            const NT signed_dist =
                ((q(0) - p(0)) * (y(1) - p(1)) - (q(1) - p(1)) * (y(0) - p(0))) / len;
            if (-signed_dist > worst) worst = -signed_dist;
        }
        return worst;
    };

    int inside = 0;
    NT max_depth = 0;
    std::vector<VecX> got;
    for (const auto& p : pts) {
        VecX y(2);
        y << p[0], p[1];
        got.push_back(y);
        if (in_hull(y)) ++inside;
        max_depth = std::max(max_depth, depth_outside(y));
    }
    std::printf("containment: %d/%zu strictly inside; worst excursion %.3e "
                "(%.2e of span)\n",
                inside, got.size(), max_depth, max_depth / span);
    CHECK(max_depth < 1e-4 * span,
          "no sample may lie MEANINGFULLY outside Sigma(A)");

    // --- distribution comparison ---------------------------------------------
    auto mean_of = [](const std::vector<VecX>& S) {
        VecX m = VecX::Zero(2);
        for (const auto& s : S) m += s;
        return VecX(m / static_cast<NT>(S.size()));
    };
    const VecX m_ref = mean_of(ref), m_got = mean_of(got);
    std::printf("mean:  reference (%.4f,%.4f)  sampled (%.4f,%.4f)  |diff|=%.4f "
                "(%.2f%% of span)\n",
                m_ref(0), m_ref(1), m_got(0), m_got(1), (m_ref - m_got).norm(),
                100.0 * (m_ref - m_got).norm() / span);
    CHECK((m_ref - m_got).norm() < 0.05 * span,
          "sampled mean must match the reference mean");

    // chi-square on a coarse grid over the bounding box
    const int G = 6;
    std::vector<int> cr(G * G, 0), cg(G * G, 0);
    auto bin = [&](const VecX& p) {
        int i = static_cast<int>((p(0) - lo(0)) / (hi(0) - lo(0)) * G);
        int j = static_cast<int>((p(1) - lo(1)) / (hi(1) - lo(1)) * G);
        i = std::min(std::max(i, 0), G - 1);
        j = std::min(std::max(j, 0), G - 1);
        return i * G + j;
    };
    for (const auto& p : ref) ++cr[static_cast<std::size_t>(bin(p))];
    for (const auto& p : got) ++cg[static_cast<std::size_t>(bin(p))];

    NT chi2 = 0;
    int cells = 0;
    for (int k = 0; k < G * G; ++k) {
        const NT e = static_cast<NT>(cr[static_cast<std::size_t>(k)]);
        const NT o = static_cast<NT>(cg[static_cast<std::size_t>(k)]);
        if (e < 20) continue;                       // pool away tiny cells
        chi2 += (o - e) * (o - e) / e;
        ++cells;
    }
    std::printf("chi-square: %.1f over %d cells (~%d dof)\n", chi2, cells,
                cells - 1);
    // Generous threshold: MCMC samples are autocorrelated, so the nominal
    // chi-square distribution understates the variance. This catches gross
    // distributional error, not fine-grained mismatch.
    CHECK(chi2 < 6.0 * static_cast<NT>(cells),
          "sampled histogram must not differ grossly from the reference");

    if (failures == 0) std::printf("test_sampler_e2e: PASS\n");
    return failures == 0 ? 0 : 1;
}
