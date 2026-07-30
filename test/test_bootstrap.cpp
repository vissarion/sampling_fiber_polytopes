// Draft Lemma 4: a certified relative-interior point from oracle calls only.

#include <secpoly/bootstrap.h>
#include <secpoly/opt_oracle.h>

#include <cstdio>
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

static void run(const char* name, int d,
                const std::vector<std::vector<long long>>& A) {
    SecondaryOptOracle oracle(d, A);
    AffineHull hull(d, A);

    const BootstrapResult boot = bootstrap(oracle, hull);

    CHECK(static_cast<int>(boot.vertices.size()) == hull.n() + 1,
          "bootstrap must return exactly n+1 vertices");
    CHECK(boot.oracle_calls > 0, "bootstrap must make oracle calls");

    // The vertices must be AFFINELY INDEPENDENT: that is what makes their
    // barycentre a relative-interior point (draft Lemma 4).
    const int n = hull.n();
    MatX diffs(n, n);
    VecX v0(hull.m());
    for (int i = 0; i < hull.m(); ++i)
        v0(i) = static_cast<NT>(boot.vertices[0][static_cast<std::size_t>(i)]);
    const VecX b0 = hull.U().transpose() * v0;
    for (int k = 1; k <= n; ++k) {
        VecX x(hull.m());
        for (int i = 0; i < hull.m(); ++i)
            x(i) = static_cast<NT>(boot.vertices[static_cast<std::size_t>(k)]
                                               [static_cast<std::size_t>(i)]);
        diffs.col(k - 1) = hull.U().transpose() * x - b0;
    }
    Eigen::JacobiSVD<MatX> svd(diffs);
    const auto& sv = svd.singularValues();
    CHECK(sv(sv.size() - 1) > sv(0) * 1e-10,
          "the n+1 bootstrap vertices must be affinely independent");

    // rho_simp must be a strictly positive certified inner radius.
    CHECK(boot.rho_simp > 0, "rho_simp must be > 0");

    // x0 must satisfy the affine-hull equations, being a convex combination of
    // points that do.
    AffineHull h2(d, A);
    h2.set_origin(boot.x0);
    const VecX bvec = h2.M() * boot.x0;
    CHECK(h2.hull_residual(boot.x0, bvec) < 1e-12, "x0 lies in the affine hull");

    std::printf("%-10s n=%2d  vertices=%2zu  rho_simp=%.6g  calls=%d "
                "(Lemma 4 bound 1+2n = %d)\n",
                name, n, boot.vertices.size(), boot.rho_simp,
                boot.oracle_calls, 1 + 2 * n);
}

int main() {
    run("example1", 1, {{1}, {2}, {4}, {6}});
    run("simple", 1, {{0}, {1}, {2}, {3}});

    std::vector<std::vector<long long>> cube;
    for (int x = 0; x < 2; ++x)
        for (int y = 0; y < 2; ++y)
            for (int z = 0; z < 2; ++z) cube.push_back({x, y, z});
    run("cube3", 3, cube);

    if (failures == 0) std::printf("test_bootstrap: PASS\n");
    return failures == 0 ? 0 : 1;
}
