// ============================================================================
// THE ONLY TRANSLATION UNIT IN THIS PROJECT THAT INCLUDES respol.
// ============================================================================
//
// res_enum_functions.h can be included exactly once program-wide: no include
// guards, non-inline globals at namespace scope (`double conv_time = 0;` :117),
// reads globals D/CD/PD that the includer must define, and #includes a .cpp
// (tropli/tropli_disc.cpp) at :27. A grep-based ctest enforces the "exactly one"
// part.
//
// D/CD/PD are declared `static` here so they get INTERNAL linkage and never
// escape into the object file's symbol table -- respol's own res_enum_d.cpp:28-30
// declares them non-static, which would export symbols literally named D, CD, PD.
// test_respol_island asserts `nm` reports none of them.

#include <secpoly/opt_oracle.h>

static int D, CD, PD;

#include <res_enum_functions.h>

#include <cassert>
#include <chrono>
#include <numeric>
#include <random>
#include <stdexcept>

namespace secpoly {

namespace {

/// respol reads D/CD/PD as globals, so every entry point must (re)install this
/// instance's values before touching respol. Cheap, and it makes multiple live
/// oracle instances safe.
struct GlobalsGuard {
    GlobalsGuard(int d, int cd, int pd) { D = d; CD = cd; PD = pd; }
};

}  // namespace

struct SecondaryOptOracle::Impl {
    int d_ = 0;
    int m_ = 0;
    int cd_ = 0;
    int pd_ = 0;
    int rd_ = 0;

    std::vector<std::vector<Field>> pointset;  // Cayley-augmented, m+d points
    std::vector<int> mi;                       // [m, 1, ..., 1], d+1 entries
    std::vector<int> proj;

    std::unique_ptr<HD> dets;
    std::unique_ptr<HD> pdets;                 // unused by compute_res_vertex2
    std::unique_ptr<Triangulation> res;        // unused by compute_res_vertex2
    std::unique_ptr<CTriangulation> tri;

    ResPol::config conf{};

    long long v_norm = 0;

    mutable std::size_t calls = 0;
    mutable std::size_t degenerate = 0;
    mutable double seconds = 0.0;
    mutable std::mt19937_64 rng{0xC0FFEEu};
};

SecondaryOptOracle::SecondaryOptOracle(int d,
                                       std::vector<std::vector<long long>> points,
                                       std::vector<int> proj)
    : p_(std::make_unique<Impl>()) {
    if (d < 1) throw std::invalid_argument("d must be >= 1");
    if (points.empty()) throw std::invalid_argument("empty point set");

    const int m = static_cast<int>(points.size());
    for (const auto& pt : points)
        if (static_cast<int>(pt.size()) != d)
            throw std::invalid_argument("every point must have exactly d coordinates");

    p_->d_ = d;
    p_->m_ = m;

    if (proj.empty()) {
        proj.resize(m);
        std::iota(proj.begin(), proj.end(), 0);
    }
    std::sort(proj.begin(), proj.end());
    p_->proj = proj;
    p_->pd_ = static_cast<int>(proj.size());
    p_->cd_ = 2 * d + 1;

    // ---- secondary encoding: A, then d copies of the origin -----------------
    // respol expresses Sigma(A) as a DEGENERATE RESULTANT. It offers no other
    // way: read_pointset hardcodes CD=2*D+1 (parse_functions.h:116) and requires
    // exactly d+1 supports (:124, exit(-1) at :130-135), and res_enum_d.cpp:130
    // applies the Cayley trick unconditionally. The shipped
    // secondary_examples/cube3.txt ("8 1 1 1" for d=3) is exactly this shape.
    //
    // With A_1 = A and A_2..A_{d+1} = {0}, the Minkowski sum is A itself, so
    // regular mixed subdivisions of the sum are regular subdivisions of A. Every
    // full-dimensional Cayley simplex must contain all d singleton pads (the only
    // points with nonzero Cayley coordinates), leaving exactly d+1 points of A --
    // a genuine cell of a triangulation of A. See DESIGN.md 3.
    p_->mi.assign(1, m);
    p_->mi.insert(p_->mi.end(), d, 1);

    p_->pointset.reserve(m + d);
    for (const auto& pt : points) {
        std::vector<Field> q;
        q.reserve(2 * d);
        for (long long v : pt) q.emplace_back(v);
        p_->pointset.push_back(std::move(q));
    }
    for (int j = 0; j < d; ++j)
        p_->pointset.emplace_back(std::vector<Field>(d, Field(0)));

    // ---- Cayley coordinates, inlined -----------------------------------------
    // Mirrors cayley_trick (res_enum_functions.h:159-185) but WITHOUT its side
    // effect: the original writes topcom_cayley.txt into the cwd at :188-190.
    // Support 0 (A's m points) gets d zeros; the j-th pad gets e_{j-1}.
    for (int i = 0; i < m; ++i)
        p_->pointset[i].insert(p_->pointset[i].end(), d, Field(0));
    for (int j = 0; j < d; ++j) {
        std::vector<Field> cayley(d, Field(0));
        cayley[j] = Field(1);
        auto& q = p_->pointset[m + j];
        q.insert(q.end(), cayley.begin(), cayley.end());
    }

    p_->rd_ = (m + d) - 2 * d - 1;  // == m - d - 1 == n
    if (p_->rd_ < 1)
        throw std::invalid_argument("m - d - 1 < 1: Sigma(A) is not full-dimensional");

    GlobalsGuard g(p_->d_, p_->cd_, p_->pd_);

    p_->conf.verbose = 0;
    p_->conf.read_from_file = false;
    p_->conf.output_f_vector = false;
    p_->conf.polytope_type = 1;  // 1 == SECONDARY (0 resultant, 2 discriminant)

    p_->dets = std::make_unique<HD>(p_->pointset.begin(), p_->pointset.end());
    p_->pdets = std::make_unique<HD>();
    p_->res = std::make_unique<Triangulation>(p_->pd_);
    p_->tri = std::make_unique<CTriangulation>(p_->cd_);

    StaticTriangulation(p_->pointset, p_->proj, *p_->tri, *p_->dets);

    // V_norm from any GKZ vector: sum of coordinates == (d+1) * V_norm.
    std::vector<long long> probe(p_->pd_, 0), gkz, used;
    probe[0] = 1;
    if (maximize_robust(probe, gkz, used)) {
        long long s = 0;
        for (long long v : gkz) s += v;
        p_->v_norm = s / (d + 1);
    }
}

SecondaryOptOracle::~SecondaryOptOracle() = default;
SecondaryOptOracle::SecondaryOptOracle(SecondaryOptOracle&&) noexcept = default;
SecondaryOptOracle& SecondaryOptOracle::operator=(SecondaryOptOracle&&) noexcept = default;

int SecondaryOptOracle::m() const { return p_->m_; }
int SecondaryOptOracle::d() const { return p_->d_; }
int SecondaryOptOracle::ambient_dim() const { return p_->pd_; }
int SecondaryOptOracle::n() const { return p_->rd_; }
long long SecondaryOptOracle::normalized_volume() const { return p_->v_norm; }
std::size_t SecondaryOptOracle::num_calls() const { return p_->calls; }
std::size_t SecondaryOptOracle::num_degenerate() const { return p_->degenerate; }
double SecondaryOptOracle::total_seconds() const { return p_->seconds; }

bool SecondaryOptOracle::maximize(const std::vector<long long>& c,
                                  std::vector<long long>& gkz_out) const {
    if (static_cast<int>(c.size()) != p_->pd_)
        throw std::invalid_argument("cost vector has wrong length");

    GlobalsGuard g(p_->d_, p_->cd_, p_->pd_);

    std::vector<Field> cf;
    cf.reserve(c.size());
    for (long long v : c) cf.emplace_back(v);
    PVector_d nli(p_->pd_, cf.begin(), cf.end());

    const auto t0 = std::chrono::steady_clock::now();
    // Res and Pdets are UNUSED in the body of compute_res_vertex2 (:865-905 --
    // the only occurrences of "Res" there are inside the debug string
    // "new Res vertex (up)"). Passing empty dummies is what keeps the
    // n-dimensional hull out of memory entirely.
    std::vector<Field> v = compute_res_vertex2<Field>(
        p_->pointset, p_->mi, p_->rd_, p_->proj,
        *p_->dets, *p_->pdets, *p_->res, *p_->tri, nli, p_->conf);
    p_->seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    ++p_->calls;

    // Empty return == degenerate lifting: the induced regular subdivision was not
    // a triangulation (project_upper_hull_r, res_enum_functions.h:518-527).
    if (v.empty()) {
        ++p_->degenerate;
        return false;
    }

    gkz_out.clear();
    gkz_out.reserve(v.size());
    for (const Field& f : v) {
        // GKZ entries are integral: the /d! normalization is commented out
        // (res_enum_functions.h:570-575), so these are un-normalized volumes.
        const CGAL::Gmpq q = f;
        if (q.denominator() != 1)
            throw std::runtime_error("GKZ coordinate is not integral");
        gkz_out.push_back(static_cast<long long>(
            CGAL::to_double(q.numerator())));
    }
    return true;
}

bool SecondaryOptOracle::maximize_robust(const std::vector<long long>& c,
                                         std::vector<long long>& gkz_out,
                                         std::vector<long long>& c_used,
                                         int max_attempts) const {
    c_used = c;
    if (maximize(c_used, gkz_out)) return true;

    // Perturb randomly, shrinking geometrically. The degenerate set is a finite
    // union of hyperplanes (the walls of the secondary fan), so a random
    // direction misses it with probability 1. Validity of the resulting halfspace
    // does NOT depend on the perturbation size -- for any c'' whatsoever,
    // (c'', c''^T v) is exactly valid. Only the FW convergence rate is affected;
    // see DESIGN.md 5.
    long long scale = 1;
    for (long long x : c) scale = std::max(scale, std::llabs(x));
    long long eps = std::max<long long>(1, scale >> 12);

    for (int attempt = 1; attempt < max_attempts; ++attempt) {
        std::uniform_int_distribution<long long> dist(-eps, eps);
        c_used = c;
        for (auto& x : c_used) x += dist(p_->rng);
        if (maximize(c_used, gkz_out)) return true;
        if (eps > 1) eps >>= 1;
    }
    return false;
}

}  // namespace secpoly
