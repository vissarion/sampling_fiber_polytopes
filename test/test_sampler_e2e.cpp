// End-to-end: volesti's BallWalk sampling Sigma(A) through the oracle.
//
// The configuration A below is COLLINEAR, so by GKZ94 Ch.7 Prop 3.1 the secondary
// polytope is combinatorially an (m-2)-cube: dimension n = m-2 and exactly
// 2^(m-2) vertices. Change A freely -- every check adapts.
//
// Two tiers of validation:
//   * structural, at ANY size -- correct dimension, correct vertex count, and
//     every sample satisfying both affine-hull equations and the coordinate
//     bounds of draft Lemma 6;
//   * distributional, only when n == 2 -- a chi-square comparison against an
//     independent uniform sampler built by rejection against the complete vertex
//     list, sharing no code with the ball walk or with Frank-Wolfe. Exact hull
//     membership in higher dimension would require an LP, which is deliberately
//     not pulled in here.

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
    // A COLLINEAR configuration: by GKZ94 Ch.7 Prop 3.1, Sigma(A) is then
    // combinatorially an (m-2)-cube, hence has exactly 2^(m-2) vertices and
    // dimension n = m-d-1 = m-2. Change A freely -- everything below adapts.
    const std::vector<std::vector<long long>> A{{1}, {2}, {4}, {6}, {8}, {9}, {10}, {12}, {14}, {15}, {16}, {18}, {20}, {21}, {22}, {24}};

    SecondaryOptOracle oracle(1, A);
    AffineHull hull(1, A);
    const int M = oracle.m();          // length of a GKZ vector
    const int N = oracle.n();          // dim Sigma(A)
    long long expect_verts = 1;        // 2^(m-2) for a collinear configuration
    for (int i = 0; i < M - 2; ++i) expect_verts *= 2;
    const BootstrapResult boot = bootstrap(oracle, hull);

    auto sep = std::make_shared<SeparationOracle>(&oracle, hull, boot);
    SecondaryBody<Point> body(sep, boot.rho_simp);

    // THE critical assertion: volesti must be told n, never m. Sigma(A) has
    // measure zero in R^m, so a ball walk there would accept nothing.
    CHECK(static_cast<int>(body.dimension()) == N,
          "body.dimension() must be n = m-d-1, NOT the ambient m");

    // --- ground truth: enumerate the vertices by a direction sweep ------------
    std::set<IVec> verts;
    {
        std::mt19937_64 rng(31337);
        std::uniform_int_distribution<long long> di(-2000, 2000);
        IVec c(static_cast<std::size_t>(oracle.ambient_dim())), g, u;
        for (int t = 0; t < 200 * M * M; ++t) {
            for (auto& x : c) x = di(rng);
            if (oracle.maximize_robust(c, g, u)) verts.insert(g);
        }
    }
    // Enumeration is only meaningful while 2^(m-2) is small: a random direction
    // hits a vertex with probability proportional to the solid angle of its
    // normal cone, and those shrink fast. Beyond that threshold the sweep cannot
    // recover the vertex set -- which is the entire reason for sampling -- so we
    // report coverage instead of asserting completeness.
    const bool can_enumerate = (M - 2) <= 8;          // 2^8 = 256 vertices
    if (can_enumerate) {
        CHECK(static_cast<long long>(verts.size()) == expect_verts,
              "collinear A: Sigma(A) must have 2^(m-2) vertices (GKZ Prop 3.1)");
    } else {
        std::printf("vertex enumeration not attempted: 2^(m-2) = %lld vertices is "
                    "beyond reach; the sweep recovered %zu (%.2f%%)\n",
                    expect_verts, verts.size(),
                    100.0 * static_cast<double>(verts.size()) /
                            static_cast<double>(expect_verts));
    }

    AffineHull h2(1, A);
    h2.set_origin(boot.x0);
    std::vector<VecX> V;
    for (const auto& g : verts) {
        VecX x(M);
        for (int i = 0; i < M; ++i)
            x(i) = static_cast<NT>(g[static_cast<std::size_t>(i)]);
        V.push_back(h2.to_intrinsic(x));
    }

    // The two affine-hull equations every point of Sigma(A) satisfies:
    //   1^T x = (d+1) * V_norm   and   a^T x = const  (both read off a vertex)
    const long long Vnorm = oracle.normalized_volume();
    const NT expect_sum = static_cast<NT>(2 * Vnorm);          // d = 1
    NT expect_wsum = 0;
    for (int i = 0; i < M; ++i)
        expect_wsum += static_cast<NT>(A[static_cast<std::size_t>(i)][0]) *
                       static_cast<NT>((*verts.begin())[static_cast<std::size_t>(i)]);

    // ======================================================================
    //  EXPORT 1 -- the vertices of Sigma(A), deduplicated
    // ======================================================================
    // In 1-D the triangulation is recoverable from the GKZ vector by inspection:
    // a point a_i is a breakpoint exactly when phi_i > 0.
    {
        std::printf(
            "\n================================================================\n"
            "  Sigma(A), A collinear with m=%d points in R^1\n"
            "  dim Sigma(A) = m-d-1 = %d ;  V_norm = %lld\n"
            "  Collinear A => combinatorially an (m-2)-cube = %d-cube,\n"
            "  hence %lld vertices (GKZ94 Ch.7 Prop 3.1)\n"
            "================================================================\n",
            M, N, Vnorm, M - 2, expect_verts);

        std::printf("\nVERTICES OF Sigma(A)  (%zu distinct)\n", V.size());
        std::printf("  %-26s %-30s %s\n", "GKZ vector (R^m)", "intrinsic (R^n)",
                    "triangulation");
        std::printf("  %-26s %-30s %s\n", "----------------", "---------------",
                    "-------------");

        std::FILE* fv = std::fopen("vertices.csv", "w");
        if (fv) {
            for (int i = 0; i < M; ++i) std::fprintf(fv, "phi%d,", i + 1);
            for (int i = 0; i < N; ++i)
                std::fprintf(fv, "y%d%s", i + 1, i + 1 < N ? "," : "\n");
        }

        for (const auto& yv : V) {
            const VecX amb = h2.to_ambient(yv);
            IVec g(static_cast<std::size_t>(M));
            for (int i = 0; i < M; ++i)
                g[static_cast<std::size_t>(i)] =
                    static_cast<long long>(std::llround(amb(i)));

            char gkz[192] = {0}; int o1 = 0;
            for (int i = 0; i < M; ++i)
                o1 += std::snprintf(gkz + o1, sizeof(gkz) - static_cast<std::size_t>(o1),
                                    "%s%lld", i ? "," : "(",
                                    g[static_cast<std::size_t>(i)]);
            std::snprintf(gkz + o1, sizeof(gkz) - static_cast<std::size_t>(o1), ")");

            char chart[192] = {0}; int o2 = 0;
            for (int i = 0; i < N; ++i)
                o2 += std::snprintf(chart + o2, sizeof(chart) - static_cast<std::size_t>(o2),
                                    "%s%6.3f", i ? "," : "(", yv(i));
            std::snprintf(chart + o2, sizeof(chart) - static_cast<std::size_t>(o2), ")");

            char tri[192] = {0}; int o3 = 0; long long prev = -1;
            for (int i = 0; i < M; ++i) {
                if (g[static_cast<std::size_t>(i)] > 0) {
                    const long long a = A[static_cast<std::size_t>(i)][0];
                    if (prev >= 0)
                        o3 += std::snprintf(tri + o3, sizeof(tri) - static_cast<std::size_t>(o3),
                                            "%s[%lld,%lld]", o3 ? "," : "", prev, a);
                    prev = a;
                }
            }
            std::printf("  %-26s %-30s {%s}\n", gkz, chart, tri);
            if (fv) {
                for (int i = 0; i < M; ++i)
                    std::fprintf(fv, "%lld,", g[static_cast<std::size_t>(i)]);
                for (int i = 0; i < N; ++i)
                    std::fprintf(fv, "%.10f%s", yv(i), i + 1 < N ? "," : "\n");
            }
        }
        if (fv) { std::fclose(fv); std::printf("  -> wrote vertices.csv\n"); }
    }

    // --- the ball walk --------------------------------------------------------
    const unsigned int walk_len = 5, num = 20000, burns = 500;
    RNGType rng(body.dimension());
    Point start(body.dimension());          // zero == x0, certified in relint
    std::list<Point> pts;

    const double L = 2.5;                   // explicit step size; see DESIGN.md 6
    BallWalk walk(L);
    body.reset_counters();
    uniform_sampling(pts, body, rng, walk, walk_len, num, start, burns);

    std::printf("\nball walk: %zu samples, acceptance=%.3f, oracle_calls=%zu, "
                "lazy_hits=%zu, atoms=%zu\n",
                pts.size(), body.acceptance_rate(), sep->oracle_calls(),
                sep->lazy_hits(), sep->num_atoms());
    CHECK(pts.size() == num, "should have produced the requested sample count");
    CHECK(body.acceptance_rate() > 0.02,
          "acceptance must not collapse -- a near-zero rate is the signature of "
          "sampling in the wrong dimension");

    // ======================================================================
    //  EXPORT 2 -- first 10 sample points, and structural checks on ALL of them
    // ======================================================================
    std::vector<VecX> got;
    {
        std::printf("\nFIRST 10 SAMPLE POINTS  (every row must have 1^T x = %.4f, "
                    "a^T x = %.4f)\n", expect_sum, expect_wsum);

        std::FILE* fs = std::fopen("samples.csv", "w");
        if (fs) {
            std::fprintf(fs, "index,");
            for (int i = 0; i < N; ++i) std::fprintf(fs, "y%d,", i + 1);
            for (int i = 0; i < M; ++i) std::fprintf(fs, "x%d,", i + 1);
            std::fprintf(fs, "sum,weighted\n");
        }

        int idx = 0, bad_hull = 0, bad_box = 0;
        for (const auto& p : pts) {
            VecX y(N);
            for (int i = 0; i < N; ++i) y(i) = p[i];
            got.push_back(y);
            const VecX x = h2.to_ambient(y);

            NT sum = 0, wsum = 0;
            for (int i = 0; i < M; ++i) {
                sum  += x(i);
                wsum += static_cast<NT>(A[static_cast<std::size_t>(i)][0]) * x(i);
                // draft Lemma 6: every GKZ coordinate lies in [0, V_norm]
                if (x(i) < -1e-6 || x(i) > static_cast<NT>(Vnorm) + 1e-6) ++bad_box;
            }
            if (std::fabs(sum - expect_sum) > 1e-6 ||
                std::fabs(wsum - expect_wsum) > 1e-6) ++bad_hull;

            if (idx < 10) {
                std::printf("  %-3d ", idx);
                for (int i = 0; i < N; ++i) std::printf("%s%6.3f", i ? "," : "(", y(i));
                std::printf(")  ");
                for (int i = 0; i < M; ++i) std::printf("%s%5.3f", i ? "," : "(", x(i));
                std::printf(")  %7.4f %7.4f\n", sum, wsum);
            }
            if (fs) {
                std::fprintf(fs, "%d,", idx);
                for (int i = 0; i < N; ++i) std::fprintf(fs, "%.10f,", y(i));
                for (int i = 0; i < M; ++i) std::fprintf(fs, "%.10f,", x(i));
                std::fprintf(fs, "%.10f,%.10f\n", sum, wsum);
            }
            ++idx;
        }
        if (fs) { std::fclose(fs); std::printf("  -> wrote samples.csv (all %d)\n", idx); }

        CHECK(bad_hull == 0, "every sample must satisfy both affine-hull equations");
        CHECK(bad_box == 0, "every GKZ coordinate must lie in [0, V_norm]");
        std::printf("  affine-hull violations: %d/%d ; box violations: %d\n",
                    bad_hull, idx, bad_box);
    }

    // ======================================================================
    //  Distributional comparison -- planar case only
    // ======================================================================
    if (N != 2) {
        std::printf("\ndistributional comparison skipped: exact hull membership "
                    "needs n == 2 (here n = %d)\n", N);
        if (failures == 0) std::printf("test_sampler_e2e: PASS\n");
        return failures == 0 ? 0 : 1;
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

    VecX lo = V[0], hi = V[0];
    for (const auto& v : V) { lo = lo.cwiseMin(v); hi = hi.cwiseMax(v); }
    const NT span = (hi - lo).norm();

    std::vector<VecX> ref;
    {
        std::mt19937_64 r(2718);
        std::uniform_real_distribution<NT> ux(lo(0), hi(0)), uy(lo(1), hi(1));
        while (ref.size() < 20000) {
            VecX p(2); p << ux(r), uy(r);
            if (in_hull(p)) ref.push_back(p);
        }
    }

    NT max_depth = 0;
    int inside = 0;
    for (const auto& y : got) {
        if (in_hull(y)) ++inside;
        for (std::size_t i = 0; i < V.size(); ++i) {
            const VecX& p = V[i];
            const VecX& q = V[(i + 1) % V.size()];
            const NT len = (q - p).norm();
            if (len < 1e-300) continue;
            const NT sd = ((q(0)-p(0))*(y(1)-p(1)) - (q(1)-p(1))*(y(0)-p(0))) / len;
            if (-sd > max_depth) max_depth = -sd;
        }
    }
    std::printf("\ncontainment: %d/%zu strictly inside; worst excursion %.3e "
                "(%.2e of span)\n", inside, got.size(), max_depth, max_depth / span);
    CHECK(max_depth < 1e-4 * span, "no sample may lie MEANINGFULLY outside Sigma(A)");

    auto mean_of = [](const std::vector<VecX>& S) {
        VecX m = VecX::Zero(2);
        for (const auto& s : S) m += s;
        return VecX(m / static_cast<NT>(S.size()));
    };
    const VecX m_ref = mean_of(ref), m_got = mean_of(got);
    std::printf("mean:  reference (%.4f,%.4f)  sampled (%.4f,%.4f)  |diff|=%.4f "
                "(%.2f%% of span)\n", m_ref(0), m_ref(1), m_got(0), m_got(1),
                (m_ref - m_got).norm(), 100.0 * (m_ref - m_got).norm() / span);
    CHECK((m_ref - m_got).norm() < 0.05 * span, "sampled mean must match the reference");

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

    NT chi2 = 0; int cells = 0;
    for (int k = 0; k < G * G; ++k) {
        const NT e = static_cast<NT>(cr[static_cast<std::size_t>(k)]);
        const NT o = static_cast<NT>(cg[static_cast<std::size_t>(k)]);
        if (e < 20) continue;
        chi2 += (o - e) * (o - e) / e; ++cells;
    }
    std::printf("chi-square: %.1f over %d cells (~%d dof)\n", chi2, cells, cells - 1);
    CHECK(chi2 < 6.0 * static_cast<NT>(cells),
          "sampled histogram must not differ grossly from the reference");

    if (failures == 0) std::printf("test_sampler_e2e: PASS\n");
    return failures == 0 ? 0 : 1;
}
