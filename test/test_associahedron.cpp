// The associahedron: a whole FAMILY of instances with closed-form ground truth.
//
// GKZ94 Ch. 7 section 3B: when A is the vertex set of a convex n-gon, Sigma(A) is
// the associahedron (Stasheff polytope). Its vertices biject with the
// triangulations of the n-gon, equivalently with the parenthesizations of a
// product x_1...x_{n-1}, so their number is the Catalan number C(n-2) -- a fact
// established long before secondary polytopes existed (Euler counted these).
//
//     dim Sigma(A) = m - d - 1 = n - 3
//     #vertices    = Catalan(n-2) = 1, 2, 5, 14, 42, 132, ...
//
// This is the most demanding validation in the suite, for three reasons:
//   * it is EXACT and independent -- Catalan numbers owe nothing to respol,
//     to our chart, or to Frank-Wolfe;
//   * it is a family, so an off-by-one or a boundary case shows up as a pattern
//     rather than as a single lucky pass;
//   * it exercises d = 2, whereas Example 1 is d = 1 and cube3 is d = 3.
//
// GKZ94 works the small cases explicitly: n=4 gives an interval, n=5 a plane
// pentagon, n=6 the 3-dimensional polytope with 14 vertices.

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

/// Catalan(k) = (1/(k+1)) * binom(2k, k), built iteratively to stay integral.
static long long catalan(int k) {
    long long c = 1;
    for (int i = 0; i < k; ++i) c = c * 2 * (2 * i + 1) / (i + 2);
    return c;
}

/// Vertices of a convex n-gon with INTEGER coordinates.
///
/// Points on the parabola y = x^2 are the natural choice: the parabola is
/// strictly convex, so all n points are genuine vertices of their hull (no three
/// collinear, nothing in the interior), the coordinates stay integral as respol
/// requires, and index order i = 0..n-1 is already the order around the
/// perimeter, which is what GKZ's associahedron statement assumes.
static std::vector<std::vector<long long>> ngon(int n) {
    std::vector<std::vector<long long>> A;
    A.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        A.push_back({i, static_cast<long long>(i) * i});
    return A;
}

static std::set<IVec> sweep(const SecondaryOptOracle& o, int trials,
                            std::uint64_t seed) {
    std::set<IVec> found;
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<long long> dist(-5000, 5000);
    IVec c(static_cast<std::size_t>(o.ambient_dim())), gkz, used;
    for (int t = 0; t < trials; ++t) {
        for (auto& x : c) x = dist(rng);
        if (o.maximize_robust(c, gkz, used)) found.insert(gkz);
    }
    return found;
}

static void run(int n, int trials) {
    const auto A = ngon(n);
    SecondaryOptOracle oracle(2, A);

    const int expect_dim = n - 3;
    const long long expect_verts = catalan(n - 2);

    CHECK(oracle.m() == n, "m must equal the number of polygon vertices");
    CHECK(oracle.n() == expect_dim, "dim Sigma(A) must be n-3");

    const std::set<IVec> V = sweep(oracle, trials, 1000u + static_cast<unsigned>(n));

    CHECK(static_cast<long long>(V.size()) == expect_verts,
          "vertex count must equal Catalan(n-2)");

    // Every returned point must be a legitimate GKZ vector: correct length, in
    // the affine hull, and coordinatewise within [0, V_norm] (draft Lemma 6).
    const long long Vn = oracle.normalized_volume();
    for (const auto& g : V) {
        CHECK(static_cast<int>(g.size()) == n, "GKZ vectors live in R^m");
        long long s = 0;
        for (long long x : g) {
            s += x;
            CHECK(x >= 0 && x <= Vn, "GKZ coordinates lie in [0, V_norm]");
        }
        CHECK(s == 3 * Vn, "1^T Phi == (d+1) * V_norm, with d = 2");
    }

    std::printf("%2d-gon: dim=%d (expect %d)  vertices=%3zu (expect Catalan(%d)=%lld)  "
                "V_norm=%lld  calls=%zu  degenerate=%zu  %s\n",
                n, oracle.n(), expect_dim, V.size(), n - 2, expect_verts, Vn,
                oracle.num_calls(), oracle.num_degenerate(),
                (static_cast<long long>(V.size()) == expect_verts &&
                 oracle.n() == expect_dim)
                    ? "OK"
                    : "MISMATCH");
}

/// The pentagon case is 2-dimensional, so we can additionally drive the full
/// separation stack over it and confirm the geometry is consistent.
static void pentagon_separation() {
    const auto A = ngon(5);
    SecondaryOptOracle oracle(2, A);
    AffineHull hull(2, A);
    const BootstrapResult boot = bootstrap(oracle, hull);

    CHECK(static_cast<int>(boot.vertices.size()) == hull.n() + 1,
          "bootstrap returns n+1 vertices on the pentagon");
    CHECK(boot.rho_simp > 0, "rho_simp > 0 on the pentagon");

    auto sep = std::make_shared<SeparationOracle>(&oracle, hull, boot);

    // Interior points: the OUTSIDE branch must be structurally unreachable.
    std::mt19937_64 rng(99);
    std::uniform_real_distribution<NT> u(0, 1);
    std::vector<VecX> atoms;
    for (const auto& g : boot.vertices) {
        VecX x(hull.m());
        for (int i = 0; i < hull.m(); ++i) x(i) = static_cast<NT>(g[static_cast<std::size_t>(i)]);
        AffineHull h2(2, A);
        h2.set_origin(boot.x0);
        atoms.push_back(h2.to_intrinsic(x));
    }

    int sound = 0, tested = 0;
    for (int t = 0; t < 200; ++t) {
        NT w[8], s = 0;
        for (std::size_t i = 0; i < atoms.size(); ++i) { w[i] = u(rng); s += w[i]; }
        VecX y = VecX::Zero(hull.n());
        for (std::size_t i = 0; i < atoms.size(); ++i) y += (w[i] / s) * atoms[i];
        ++tested;
        if (sep->separate(y).verdict != Verdict::Outside) ++sound;
    }
    CHECK(sound == tested,
          "SOUNDNESS on the pentagon: interior points are never declared outside");

    // Far exterior points must be separated, and each cut must be a real proof.
    const std::set<IVec> allv = sweep(oracle, 400, 5u);
    int cuts = 0, valid = 0;
    std::uniform_real_distribution<NT> far(30, 90);
    for (int t = 0; t < 200; ++t) {
        VecX y(hull.n());
        for (int i = 0; i < hull.n(); ++i) y(i) = far(rng) * (u(rng) < 0.5 ? -1 : 1);
        const SeparationResult r = sep->separate(y);
        if (r.verdict != Verdict::Outside) continue;
        ++cuts;
        bool ok = true;
        for (const auto& g : allv)
            if (exact_dot(r.cut_c, g) > r.cut_h) ok = false;
        if (ok) ++valid;
    }
    CHECK(cuts > 0, "far exterior points should be separated");
    CHECK(valid == cuts, "every pentagon cut must dominate all vertices");
    std::printf("pentagon: soundness %d/%d, certificates %d/%d, atoms=%zu\n",
                sound, tested, valid, cuts, sep->num_atoms());
}

int main() {
    // Trials scale with the vertex count: a random direction hits a vertex with
    // probability proportional to the solid angle of its normal cone, and those
    // cones get small as the family grows.
    run(4, 300);     // interval,          Catalan(2) = 2
    run(5, 800);     // plane pentagon,    Catalan(3) = 5
    run(6, 4000);    // 3-dimensional,     Catalan(4) = 14
    run(7, 14000);   // 4-dimensional,     Catalan(5) = 42

    pentagon_separation();

    if (failures == 0) std::printf("test_associahedron: PASS\n");
    return failures == 0 ? 0 : 1;
}
