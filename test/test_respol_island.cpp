// Step-1 smoke test: proves the respol island builds, links and runs.
//
// This exists because the plan's biggest unknown was whether respol's patched
// CGAL Kernel_d + LEDA could be driven from a foreign CMake project at all. It
// was first validated in a scratch directory; this file makes that check
// reproducible in-tree, so it is a test rather than an assertion in a document.
//
// Companion check (see test/CMakeLists.txt): `nm` must report zero exported
// D/CD/PD symbols from libsecpoly_optoracle.a.

#include <secpoly/opt_oracle.h>

#include <cstdio>
#include <numeric>
#include <vector>

static int failures = 0;

#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, (msg));    \
            ++failures;                                                    \
        }                                                                  \
    } while (0)

int main() {
    // The draft's Example 1: A = {1,2,4,6} in R^1. m=4, d=1, n = m-d-1 = 2.
    const std::vector<std::vector<long long>> A{{1}, {2}, {4}, {6}};
    secpoly::SecondaryOptOracle oracle(/*d=*/1, A);

    CHECK(oracle.m() == 4, "m should be 4");
    CHECK(oracle.d() == 1, "d should be 1");
    CHECK(oracle.ambient_dim() == 4, "GKZ vectors live in R^4");
    CHECK(oracle.n() == 2, "dim Sigma(A) should be m-d-1 = 2");

    // V_norm = d! * vol(conv A) = 1 * (6-1) = 5.
    CHECK(oracle.normalized_volume() == 5, "V_norm should be 5");

    // One oracle call must return a GKZ vector satisfying the affine-hull
    // constraints of Sigma(A): 1^T x = (d+1) * V_norm = 10, and
    // (1,2,4,6)^T x = (d+1) * V_norm * barycentre-weighted sum = 35.
    std::vector<long long> c{0, 1, 0, 0}, gkz, used;
    const bool ok = oracle.maximize_robust(c, gkz, used);
    CHECK(ok, "maximize_robust should succeed after perturbation retries");

    if (ok) {
        CHECK(gkz.size() == 4, "GKZ vector should have 4 entries");
        const long long s = std::accumulate(gkz.begin(), gkz.end(), 0LL);
        CHECK(s == 10, "1^T Phi should equal (d+1)*V_norm = 10");

        long long w = 0;
        const long long a[4] = {1, 2, 4, 6};
        for (int i = 0; i < 4; ++i) w += a[i] * gkz[i];
        CHECK(w == 35, "A^T Phi should equal 35");

        for (long long v : gkz)
            CHECK(v >= 0 && v <= oracle.normalized_volume(),
                  "every GKZ coordinate must lie in [0, V_norm] (draft Lemma 6)");

        std::printf("GKZ = (%lld,%lld,%lld,%lld)  calls=%zu degenerate=%zu\n",
                    gkz[0], gkz[1], gkz[2], gkz[3],
                    oracle.num_calls(), oracle.num_degenerate());
    }

    if (failures == 0) std::printf("test_respol_island: PASS\n");
    return failures == 0 ? 0 : 1;
}
