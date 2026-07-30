// The intrinsic chart. This is what stops volesti from ever seeing R^m.
//
// Sigma(A) has dimension n = m-d-1 but lives in R^m -- for A = {1,2,4,6} it is a
// quadrilateral in R^4, of Lebesgue measure zero. Handing a ball walk an
// m-dimensional body is not merely slow: it proposes points in an m-ball and
// accepts iff is_in, which happens with probability 0, so the chain never moves
// and silently returns copies of the start. See DESIGN.md 6.

#ifndef SECPOLY_AFFINE_HULL_H
#define SECPOLY_AFFINE_HULL_H

#include <secpoly/types.h>

#include <stdexcept>

namespace secpoly {

/// Orthonormal chart for the affine hull of Sigma(A).
///
/// Every GKZ vector satisfies  1^T x = (d+1) V  and  A_mat x = (d+1) V abar,
/// i.e. M x = b with M = [1^T ; A_mat] of shape (d+1) x m. So Sigma(A) lies in an
/// affine subspace with direction space L0 = ker(M), of dimension m-(d+1) = n.
///
/// U is an ORTHONORMAL basis of L0, which matters for two reasons -- and the
/// obvious one is not among them. Uniformity does NOT require orthonormality: any
/// affine bijection has constant Jacobian, so uniform pushes forward to uniform
/// under any linear chart. Orthonormality makes y -> x0 + U y an ISOMETRY onto the
/// affine hull, which (a) preserves the body's shape, so the ball walk is not
/// crawling across an artificially elongated body, and (b) lets the Euclidean
/// radius bounds rho_simp and V*sqrt(m) transfer unchanged.
class AffineHull {
public:
    /// \param d      ambient dimension of the point configuration
    /// \param points the m points of A (each of length d)
    /// \throws std::invalid_argument if A does not affinely span R^d, since then
    ///         dim Sigma(A) != m-d-1 and the whole encoding is wrong. Fail loudly
    ///         rather than sample a body of the wrong dimension.
    AffineHull(int d, const std::vector<std::vector<long long>>& points) {
        m_ = static_cast<int>(points.size());
        d_ = d;
        n_ = m_ - d_ - 1;
        if (n_ < 1)
            throw std::invalid_argument("m - d - 1 < 1: Sigma(A) is not full-dimensional");

        // M = [1^T ; A_mat],  (d+1) x m
        M_.resize(d_ + 1, m_);
        for (int j = 0; j < m_; ++j) {
            M_(0, j) = 1.0;
            for (int i = 0; i < d_; ++i)
                M_(i + 1, j) = static_cast<NT>(points[j][i]);
        }

        // Rank check with a RELATIVE threshold -- an absolute one is meaningless
        // once coordinates are large.
        Eigen::JacobiSVD<MatX> svd(M_);
        const auto& sv = svd.singularValues();
        const NT tol = sv(0) * 1e-10;
        int rank = 0;
        for (int i = 0; i < sv.size(); ++i)
            if (sv(i) > tol) ++rank;
        if (rank != d_ + 1)
            throw std::invalid_argument(
                "A does not affinely span R^d (rank[1;A] != d+1), so "
                "dim Sigma(A) != m-d-1");

        // Orthonormal basis of ker(M) = the LAST n columns of the full Q factor
        // of M^T. For M^T = QR with M^T of shape m x (d+1) and full column rank,
        // the first d+1 columns of Q span the row space of M, so the remaining
        // m-(d+1) = n columns span its kernel.
        Eigen::HouseholderQR<MatX> qr(M_.transpose());
        const MatX Q = qr.householderQ();
        U_ = Q.rightCols(n_);
    }

    int m() const { return m_; }
    int d() const { return d_; }
    int n() const { return n_; }

    const MatX& U() const { return U_; }
    const MatX& M() const { return M_; }

    void set_origin(const VecX& x0) { x0_ = x0; }
    const VecX& origin() const { return x0_; }

    /// ambient (R^m) -> intrinsic (R^n)
    VecX to_intrinsic(const VecX& x) const { return U_.transpose() * (x - x0_); }

    /// intrinsic (R^n) -> ambient (R^m). Lands in the affine hull BY
    /// CONSTRUCTION, so no drift can accumulate.
    VecX to_ambient(const VecX& y) const { return x0_ + U_ * y; }

    /// Project a cost vector onto L0, i.e. reduce it modulo the affine functions.
    ///
    /// Adding an affine function to a lifting provably does not change the induced
    /// regular subdivision, and the affine functions evaluated on A are exactly
    /// the row space of M. So this is EXACTLY free -- it cannot change which GKZ
    /// vector maximizes the cost -- and it shrinks the dynamic range before
    /// rounding, which is why rationalize() calls it first.
    VecX reduce_cost(const VecX& c) const { return U_ * (U_.transpose() * c); }

    /// Diagnostic: how far a supposedly-ambient point strays from the affine hull.
    /// Should sit at ~1e-15 for anything produced by to_ambient().
    NT hull_residual(const VecX& x, const VecX& b) const {
        return (M_ * x - b).cwiseAbs().maxCoeff() /
               std::max(NT(1), b.cwiseAbs().maxCoeff());
    }

private:
    int m_ = 0, d_ = 0, n_ = 0;
    MatX M_, U_;
    VecX x0_;
};

}  // namespace secpoly

#endif  // SECPOLY_AFFINE_HULL_H
