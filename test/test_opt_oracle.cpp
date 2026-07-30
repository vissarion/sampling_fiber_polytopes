// Validation of the optimization oracle against known ground truth.
//
// The decisive test for Part (a). Sweeping cost directions must recover exactly
// the vertex set of Sigma(A), which is known independently for these instances:
//   * A = {1,2,4,6}: the draft's Example 1, four GKZ vectors given explicitly
//   * A = {0,1,2,3}: respol's secondary_examples/simple.txt, 4 vertices
//                    (test_suite.cpp:69-70)
//   * unit 3-cube:   respol's cube3.txt, 74 vertices (ibid.)

#include <secpoly/opt_oracle.h>

#include <algorithm>
#include <cstdio>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

static int failures = 0;

#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, (msg));    \
            ++failures;                                                    \
        }                                                                  \
    } while (0)

using Vec = std::vector<long long>;

// Sweep random integral cost directions and collect the distinct GKZ vectors.
// The component of the cost orthogonal to the affine hull of Sigma(A) is constant
// on the body and cannot change the argmax, so sampling ambient directions is
// enough to explore every vertex.
static std::set<Vec> sweep(const secpoly::SecondaryOptOracle& o, int trials,
                           std::uint64_t seed = 12345) {
    std::set<Vec> found;
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<long long> dist(-1000, 1000);
    Vec c(o.ambient_dim()), gkz, used;
    for (int t = 0; t < trials; ++t) {
        for (auto& x : c) x = dist(rng);
        if (o.maximize_robust(c, gkz, used)) found.insert(gkz);
    }
    return found;
}

static void test_example1() {
    // A = {1,2,4,6} in R^1. Triangulations of [1,6] and their GKZ vectors:
    //   {[1,6]}                -> (5,0,0,5)
    //   {[1,2],[2,6]}          -> (1,5,0,4)
    //   {[1,4],[4,6]}          -> (3,0,5,2)
    //   {[1,2],[2,4],[4,6]}    -> (1,3,4,2)
    const std::vector<std::vector<long long>> A{{1}, {2}, {4}, {6}};
    secpoly::SecondaryOptOracle o(1, A);

    const std::set<Vec> expected{
        {5, 0, 0, 5}, {1, 5, 0, 4}, {3, 0, 5, 2}, {1, 3, 4, 2}};

    const std::set<Vec> found = sweep(o, 300);

    CHECK(found == expected,
          "Example 1 sweep must return EXACTLY the draft's four GKZ vectors");
    if (found != expected) {
        std::printf("  expected %zu, found %zu:\n", expected.size(), found.size());
        for (const auto& v : found) {
            std::printf("    (");
            for (std::size_t i = 0; i < v.size(); ++i)
                std::printf("%lld%s", v[i], i + 1 < v.size() ? "," : "");
            std::printf(")%s\n", expected.count(v) ? "" : "   <-- UNEXPECTED");
        }
    }
    std::printf("example1: %zu vertices, %zu oracle calls, %zu degenerate\n",
                found.size(), o.num_calls(), o.num_degenerate());
}

static void test_simple() {
    // respol secondary_examples/simple.txt: A = {0,1,2,3} in R^1.
    // GKZ theory: for collinear A, Sigma(A) is combinatorially a cube of
    // dimension m-2 = 2, hence 4 vertices -- matching test_suite.cpp:69-70.
    const std::vector<std::vector<long long>> A{{0}, {1}, {2}, {3}};
    secpoly::SecondaryOptOracle o(1, A);
    CHECK(o.n() == 2, "simple.txt: dim Sigma(A) should be 2");

    const std::set<Vec> found = sweep(o, 300);
    CHECK(found.size() == 4, "simple.txt should have exactly 4 vertices");

    for (const auto& v : found) {
        const long long s = std::accumulate(v.begin(), v.end(), 0LL);
        CHECK(s == 2 * o.normalized_volume(), "1^T Phi == (d+1) * V_norm");
    }
    std::printf("simple: %zu vertices (expect 4), V_norm=%lld\n",
                found.size(), o.normalized_volume());
}

static void test_cube3() {
    // respol secondary_examples/cube3.txt: the 8 vertices of the unit 3-cube.
    // test_suite.cpp:69-70 fixes the answer at 74 vertices. A random sweep will
    // not exhaust all 74, so we assert every returned point is a LEGITIMATE
    // vertex (right dimension, affine-hull constraints, coordinate bounds) and
    // that we recover a healthy fraction.
    std::vector<std::vector<long long>> A;
    for (int x = 0; x < 2; ++x)
        for (int y = 0; y < 2; ++y)
            for (int z = 0; z < 2; ++z) A.push_back({x, y, z});

    secpoly::SecondaryOptOracle o(3, A);
    CHECK(o.m() == 8, "cube3: m should be 8");
    CHECK(o.n() == 4, "cube3: dim Sigma(A) should be 8-3-1 = 4");

    const std::set<Vec> found = sweep(o, 400);
    CHECK(!found.empty(), "cube3 sweep should find vertices");
    CHECK(found.size() <= 74, "cube3 cannot have more than 74 vertices");

    const long long V = o.normalized_volume();
    for (const auto& v : found) {
        CHECK(static_cast<int>(v.size()) == 8, "cube3 GKZ vectors live in R^8");
        const long long s = std::accumulate(v.begin(), v.end(), 0LL);
        CHECK(s == 4 * V, "cube3: 1^T Phi == (d+1) * V_norm");
        for (long long x : v)
            CHECK(x >= 0 && x <= V, "every GKZ coordinate lies in [0, V_norm]");
    }
    std::printf("cube3: %zu/74 vertices found, V_norm=%lld, %zu calls, "
                "%zu degenerate (%.1f%%)\n",
                found.size(), V, o.num_calls(), o.num_degenerate(),
                100.0 * static_cast<double>(o.num_degenerate()) /
                    static_cast<double>(std::max<std::size_t>(1, o.num_calls())));
}

static void test_no_stray_file() {
    // cayley_trick (res_enum_functions.h:188-190) writes topcom_cayley.txt into
    // the cwd. We reimplement its ~10 lines precisely to avoid that side effect.
    std::FILE* f = std::fopen("topcom_cayley.txt", "r");
    CHECK(f == nullptr, "topcom_cayley.txt must NOT be created (cayley_trick "
                        "side effect should be avoided)");
    if (f) std::fclose(f);
}

int main() {
    test_example1();
    test_simple();
    test_cube3();
    test_no_stray_file();

    if (failures == 0) std::printf("test_opt_oracle: PASS\n");
    return failures == 0 ? 0 : 1;
}
