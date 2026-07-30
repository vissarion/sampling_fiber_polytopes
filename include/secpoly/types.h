#ifndef SECPOLY_TYPES_H
#define SECPOLY_TYPES_H

#include <Eigen/Dense>
#include <cstdint>
#include <vector>

namespace secpoly {

using NT = double;
using VecX = Eigen::VectorXd;
using MatX = Eigen::MatrixXd;
using IVec = std::vector<long long>;   ///< exact ambient GKZ vector / cost vector

/// Result of a membership query.
///
/// OUTSIDE carries an unconditional proof: `cut_c` is the integral cost vector
/// actually sent to the oracle and `cut_h = max_{Sigma(A)} cut_c^T x`, both exact,
/// so `cut_c^T x <= cut_h` holds for every point of Sigma(A) while the query point
/// violates it.
///
/// INSIDE is weak in the GLS sense: the witness is an explicit convex combination
/// of GKZ vectors lying within `tol` of the query. That is the model the draft's
/// Theorem 8 already assumes.
enum class Verdict { Inside, Outside, Ambiguous };

struct SeparationResult {
    Verdict verdict = Verdict::Ambiguous;

    IVec   cut_c;      ///< OUTSIDE: separating direction (ambient, integral)
    long long cut_h = 0;  ///< OUTSIDE: max of cut_c over Sigma(A), exact

    VecX   witness;    ///< INSIDE: a point of Sigma(A) within tol of the query
    NT     distance = 0;  ///< INSIDE: ||witness - query||

    int    oracle_calls = 0;
    int    fw_iterations = 0;
};

struct Tolerances {
    /// INSIDE is declared once Frank-Wolfe drives ||y - x|| below this.
    ///
    /// Do NOT tighten this to something like 1e-9: Frank-Wolfe is a first-order
    /// method, so demanding 9 digits costs a great many iterations, and any point
    /// that exhausts the iteration budget comes back Ambiguous. Since interior
    /// points near the boundary are exactly the slow ones, that turns into a
    /// systematic REJECTION BIAS and the sampler stops being uniform -- observed
    /// as a chi-square of ~4000 against a reference sampler.
    ///
    /// The bias from a looser tolerance is only O(tol * surface/volume), i.e.
    /// utterly negligible next to MCMC error.
    NT membership_tol = 1e-7;
    int max_fw_iters  = 500;
};

}  // namespace secpoly

#endif  // SECPOLY_TYPES_H
