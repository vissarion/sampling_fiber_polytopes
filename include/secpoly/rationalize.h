// The double <-> exact-integer boundary.
//
// The sampler runs in double; respol computes in exact rationals. Cost vectors
// cross outward, GKZ vectors cross inward. Keeping the outward direction INTEGRAL
// is what makes every separation certificate exact: whatever integer vector c we
// hand respol, it returns the true maximizer of c, so (c, c^T v) is an exactly
// valid halfspace for Sigma(A) -- unconditionally, however c was rounded or
// perturbed. See DESIGN.md 5.

#ifndef SECPOLY_RATIONALIZE_H
#define SECPOLY_RATIONALIZE_H

#include <secpoly/affine_hull.h>
#include <secpoly/types.h>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace secpoly {

/// Bits of headroom for the rounded cost vector. 2^50 leaves room for exact
/// int64 dot products against GKZ coordinates bounded by V_norm.
inline constexpr int kDefaultCostBits = 50;

/// Convert a real cost direction into an integral one.
///
/// Integer lifts keep every determinant respol computes integral, which is by far
/// the fastest path under USE_HASHED_DETERMINANTS. Exact dyadic conversion would
/// be error-free but produces enormous mixed denominators for no benefit.
///
/// \param c    real cost vector (ambient, length m)
/// \param hull used to reduce c modulo the affine functions first -- exactly free,
///             and it shrinks the dynamic range before rounding
/// \param bits target magnitude, 2^bits
inline IVec rationalize(const VecX& c, const AffineHull& hull,
                        int bits = kDefaultCostBits) {
    const VecX r = hull.reduce_cost(c);
    const NT scale = r.cwiseAbs().maxCoeff();

    IVec out(static_cast<std::size_t>(r.size()), 0);
    if (!(scale > 0) || !std::isfinite(scale)) return out;  // caller checks all-zero

    // Target max|out| ~ 2^bits.
    const NT factor = std::ldexp(NT(1), bits) / scale;
    for (Eigen::Index i = 0; i < r.size(); ++i) {
        const NT v = r(i) * factor;
        if (!std::isfinite(v) || std::fabs(v) > 9.2e18)
            throw std::overflow_error("cost vector overflows int64 after scaling");
        out[static_cast<std::size_t>(i)] = std::llround(v);
    }
    return out;
}

inline bool is_zero(const IVec& c) {
    for (long long x : c)
        if (x != 0) return false;
    return true;
}

/// Exact integer inner product, with overflow detection.
///
/// Used for the separation test c^T s < c^T y and for the cached-halfspace check,
/// both of which must be EXACT for the OUTSIDE verdict to be a genuine proof
/// rather than a tolerance-dependent guess.
inline long long exact_dot(const IVec& a, const IVec& b) {
    if (a.size() != b.size())
        throw std::invalid_argument("exact_dot: size mismatch");
    long long acc = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        long long term;
        if (__builtin_mul_overflow(a[i], b[i], &term) ||
            __builtin_add_overflow(acc, term, &acc))
            throw std::overflow_error(
                "exact_dot overflowed int64 -- reduce the cost bit width or "
                "rescale the configuration (see DESIGN.md 5)");
    }
    return acc;
}

/// Inner product of an integral cost with a real ambient point.
/// Inexact by nature: the query point is a double. Used only where a margin is
/// checked, never as the basis of an exactness claim.
inline NT mixed_dot(const IVec& c, const VecX& x) {
    NT acc = 0;
    for (std::size_t i = 0; i < c.size(); ++i)
        acc += static_cast<NT>(c[i]) * x(static_cast<Eigen::Index>(i));
    return acc;
}

inline VecX to_vec(const IVec& v) {
    VecX out(static_cast<Eigen::Index>(v.size()));
    for (std::size_t i = 0; i < v.size(); ++i)
        out(static_cast<Eigen::Index>(i)) = static_cast<NT>(v[i]);
    return out;
}

}  // namespace secpoly

#endif  // SECPOLY_RATIONALIZE_H
