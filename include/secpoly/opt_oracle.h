// Public facade for the respol-backed optimization oracle over Sigma(A).
//
// THIS HEADER MUST CONTAIN NO CGAL, NO LEDA, NO EIGEN, NO respol TYPES.
// It is the boundary that lets the rest of the project (and volesti) compile
// without ever seeing respol's patched CGAL Kernel_d. Enforced by
// test/test_no_cgal_leak.cpp and a grep-based ctest.

#ifndef SECPOLY_OPT_ORACLE_H
#define SECPOLY_OPT_ORACLE_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace secpoly {

/// Linear optimization oracle over the secondary polytope Sigma(A) of a point
/// configuration A subset Z^d.
///
/// Sigma(A) is never materialized. The only operation is: given an integral cost
/// vector c, return the GKZ vector of a triangulation maximizing c. That GKZ
/// vector is a vertex of Sigma(A) whenever c induces a genuine triangulation.
///
/// Coordinates are integral throughout: respol's `/d!` normalization is commented
/// out (res_enum_functions.h:570-575), so GKZ entries are un-normalized volumes.
class SecondaryOptOracle {
public:
    /// \param d      ambient dimension of the point configuration
    /// \param points the m points of A, each of length d. Must affinely span R^d.
    /// \param proj   which coordinates of the GKZ vector to keep. Pass an empty
    ///               vector for the full secondary polytope ([0..m-1]). A strict
    ///               subset yields a *projection* of Sigma(A) -- respol's
    ///               machinery is already a projection oracle, so this is kept
    ///               settable to leave the fiber-polytope case reachable.
    SecondaryOptOracle(int d,
                       std::vector<std::vector<long long>> points,
                       std::vector<int> proj = {});

    ~SecondaryOptOracle();
    SecondaryOptOracle(SecondaryOptOracle&&) noexcept;
    SecondaryOptOracle& operator=(SecondaryOptOracle&&) noexcept;
    SecondaryOptOracle(const SecondaryOptOracle&) = delete;
    SecondaryOptOracle& operator=(const SecondaryOptOracle&) = delete;

    int m() const;          ///< number of points in A
    int d() const;          ///< ambient dimension of A
    int ambient_dim() const;///< length of a GKZ vector == proj.size()
    int n() const;          ///< dim Sigma(A) == m - d - 1

    /// Maximize c over Sigma(A).
    ///
    /// \param c       cost vector, length ambient_dim()
    /// \param gkz_out on success, the maximizing GKZ vector (length ambient_dim())
    /// \return true on success. **false means the lifting c was degenerate** --
    ///         i.e. it induced a regular subdivision that is not a triangulation,
    ///         and respol returned the empty vector
    ///         (project_upper_hull_r, res_enum_functions.h:518-527). The caller
    ///         must perturb c and retry; see maximize_robust().
    ///
    /// When it returns true, gkz_out is the EXACT maximizer of c. Therefore
    /// (c, c^T gkz_out) is an exactly valid halfspace for Sigma(A) for ANY c
    /// whatsoever -- which is what makes the separation certificates unconditional
    /// regardless of how c was rounded or perturbed.
    bool maximize(const std::vector<long long>& c,
                  std::vector<long long>& gkz_out) const;

    /// maximize() with automatic perturb-and-retry on degeneracy.
    ///
    /// Perturbs by a random integral vector with geometrically shrinking
    /// magnitude. Randomized rather than a fixed direction because the degenerate
    /// set is a finite union of hyperplanes (the walls of the secondary fan),
    /// which a random direction misses with probability 1.
    ///
    /// \param c_used on success, the cost vector actually used -- possibly a
    ///        perturbation of c. Callers MUST run their separation test against
    ///        c_used, not against c, to keep the test exact.
    bool maximize_robust(const std::vector<long long>& c,
                         std::vector<long long>& gkz_out,
                         std::vector<long long>& c_used,
                         int max_attempts = 30) const;

    /// Normalized volume V_norm = d! * vol(conv A), computed exactly.
    /// Every GKZ coordinate lies in [0, V_norm] (draft Lemma 6), which gives the
    /// cheap ambient box test and the outer radius bound.
    long long normalized_volume() const;

    std::size_t num_calls() const;
    std::size_t num_degenerate() const;
    double total_seconds() const;

private:
    struct Impl;
    std::unique_ptr<Impl> p_;
};

}  // namespace secpoly

#endif  // SECPOLY_OPT_ORACLE_H
