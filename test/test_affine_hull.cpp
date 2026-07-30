// The intrinsic chart and the exact/inexact boundary.
//
// Includes a NEGATIVE test pinning the failure mode the chart exists to prevent:
// a body that reports dimension() = m instead of n accepts nothing, because
// Sigma(A) has measure zero in R^m.

#include <secpoly/affine_hull.h>
#include <secpoly/opt_oracle.h>
#include <secpoly/rationalize.h>

#include <cstdio>
#include <random>
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

static const std::vector<std::vector<long long>> kExample1{{1}, {2}, {4}, {6}};

// The draft's four GKZ vectors for A = {1,2,4,6}.
static const std::vector<std::vector<long long>> kGkz{
    {5, 0, 0, 5}, {1, 5, 0, 4}, {3, 0, 5, 2}, {1, 3, 4, 2}};

static void test_chart_structure() {
    AffineHull h(1, kExample1);
    CHECK(h.m() == 4, "m == 4");
    CHECK(h.n() == 2, "n == m-d-1 == 2");
    CHECK(h.U().rows() == 4 && h.U().cols() == 2, "U is m x n");

    // U^T U == I_n : the chart is an isometry, which is what preserves the body's
    // shape and lets the Euclidean radius bounds transfer.
    const MatX gram = h.U().transpose() * h.U();
    CHECK((gram - MatX::Identity(2, 2)).cwiseAbs().maxCoeff() < 1e-12,
          "U^T U must equal I_n");

    // M U == 0 : U really spans ker(M), not just some 2-dimensional subspace.
    const MatX MU = h.M() * h.U();
    CHECK(MU.cwiseAbs().maxCoeff() < 1e-10, "M U must be 0 (U spans ker M)");
}

static void test_roundtrip_on_real_gkz_vectors() {
    AffineHull h(1, kExample1);

    // x0 = barycentre of the four true GKZ vectors: genuinely in relint Sigma(A).
    VecX x0 = VecX::Zero(4);
    for (const auto& g : kGkz)
        for (int i = 0; i < 4; ++i) x0(i) += static_cast<NT>(g[i]) / 4.0;
    h.set_origin(x0);

    // b = M x0. Every GKZ vector satisfies M x = b, with b = (10, 35) here.
    const VecX b = h.M() * x0;
    CHECK(std::fabs(b(0) - 10.0) < 1e-9, "1^T Phi should be 10");
    CHECK(std::fabs(b(1) - 35.0) < 1e-9, "A^T Phi should be 35");

    for (const auto& g : kGkz) {
        VecX x(4);
        for (int i = 0; i < 4; ++i) x(i) = static_cast<NT>(g[i]);

        CHECK(h.hull_residual(x, b) < 1e-12,
              "every GKZ vector must lie in the affine hull");

        const VecX y = h.to_intrinsic(x);
        CHECK(y.size() == 2, "intrinsic coordinates live in R^n = R^2");
        const VecX back = h.to_ambient(y);
        CHECK((back - x).cwiseAbs().maxCoeff() < 1e-9,
              "to_ambient(to_intrinsic(x)) must be the identity on the hull");
    }
}

static void test_cost_reduction_is_free() {
    // Reducing a cost modulo the affine functions must not change WHICH GKZ
    // vector maximizes it -- that is the claim licensing the reduction. Check it
    // against the live oracle rather than by re-deriving the algebra.
    AffineHull h(1, kExample1);
    SecondaryOptOracle o(1, kExample1);

    VecX x0 = VecX::Zero(4);
    for (const auto& g : kGkz)
        for (int i = 0; i < 4; ++i) x0(i) += static_cast<NT>(g[i]) / 4.0;
    h.set_origin(x0);

    std::mt19937_64 rng(999);
    std::uniform_real_distribution<NT> dist(-100.0, 100.0);

    int agree = 0, tested = 0;
    for (int t = 0; t < 40; ++t) {
        VecX c(4);
        for (int i = 0; i < 4; ++i) c(i) = dist(rng);

        const IVec raw = rationalize(c, h);       // reduce_cost is applied inside
        const VecX reduced = h.reduce_cost(c);
        // reduce is a projection, hence idempotent
        CHECK((h.reduce_cost(reduced) - reduced).cwiseAbs().maxCoeff() < 1e-9,
              "reduce_cost must be idempotent (it is a projection)");
        // and its image lies in ker(M)
        CHECK((h.M() * reduced).cwiseAbs().maxCoeff() < 1e-8,
              "reduce_cost must land in ker(M) = L0");

        // Unreduced vs reduced cost must select the same vertex.
        IVec unreduced(4);
        {
            const NT s = c.cwiseAbs().maxCoeff();
            for (int i = 0; i < 4; ++i)
                unreduced[i] = std::llround(c(i) / s * 1e9);
        }
        IVec g1, g2, u1, u2;
        const bool ok1 = o.maximize_robust(raw, g1, u1);
        const bool ok2 = o.maximize_robust(unreduced, g2, u2);
        if (ok1 && ok2) {
            ++tested;
            if (g1 == g2) ++agree;
        }
    }
    CHECK(tested > 0, "cost-reduction cross-check should have run");
    CHECK(agree == tested,
          "reducing a cost mod the affine functions must not change the argmax");
    std::printf("cost reduction: %d/%d agree with the unreduced cost\n",
                agree, tested);
}

static void test_exact_dot() {
    const IVec a{1, -2, 3, -4}, b{10, 20, 30, 40};
    CHECK(exact_dot(a, b) == 10 - 40 + 90 - 160, "exact_dot arithmetic");

    bool threw = false;
    try {
        const IVec big{1LL << 62, 1LL << 62}, two{4, 4};
        (void)exact_dot(big, two);
    } catch (const std::overflow_error&) {
        threw = true;
    }
    CHECK(threw, "exact_dot must detect int64 overflow rather than wrap silently");
}

static void test_measure_zero_failure_mode() {
    // The reason the chart exists. Sigma(A) is a 2-dimensional quadrilateral
    // sitting in R^4, so a ball of ANY radius drawn in R^4 hits it with
    // probability 0. Simulate the mistake of sampling in the ambient dimension:
    // membership must essentially never fire.
    AffineHull h(1, kExample1);
    VecX x0 = VecX::Zero(4);
    for (const auto& g : kGkz)
        for (int i = 0; i < 4; ++i) x0(i) += static_cast<NT>(g[i]) / 4.0;
    h.set_origin(x0);
    const VecX b = h.M() * x0;

    std::mt19937_64 rng(7);
    std::normal_distribution<NT> gauss(0.0, 1.0);

    int on_hull = 0;
    for (int t = 0; t < 20000; ++t) {
        VecX x = x0;
        for (int i = 0; i < 4; ++i) x(i) += gauss(rng);     // ambient perturbation
        if (h.hull_residual(x, b) < 1e-9) ++on_hull;
    }
    CHECK(on_hull == 0,
          "ambient sampling must NEVER land on the affine hull -- this is why "
          "dimension() must return n, not m");

    // ... whereas the intrinsic chart lands on it every single time.
    int intrinsic_on_hull = 0;
    for (int t = 0; t < 2000; ++t) {
        VecX y(2);
        y << gauss(rng), gauss(rng);
        if (h.hull_residual(h.to_ambient(y), b) < 1e-9) ++intrinsic_on_hull;
    }
    CHECK(intrinsic_on_hull == 2000,
          "intrinsic points must land in the affine hull by construction");
    std::printf("measure-zero check: ambient 0/20000 on hull, "
                "intrinsic %d/2000 on hull\n", intrinsic_on_hull);
}

static void test_rank_deficient_rejected() {
    // A that fails to affinely span R^2: three collinear points in the plane.
    // Then dim Sigma(A) != m-d-1 and the encoding would be wrong, so construction
    // must fail loudly instead of silently sampling the wrong-dimensional body.
    const std::vector<std::vector<long long>> flat{{0, 0}, {1, 1}, {2, 2}, {3, 3}};
    bool threw = false;
    try {
        AffineHull h(2, flat);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    CHECK(threw, "a configuration not affinely spanning R^d must be rejected");
}

int main() {
    test_chart_structure();
    test_roundtrip_on_real_gkz_vectors();
    test_cost_reduction_is_free();
    test_exact_dot();
    test_measure_zero_failure_mode();
    test_rank_deficient_rejected();

    if (failures == 0) std::printf("test_affine_hull: PASS\n");
    return failures == 0 ? 0 : 1;
}
