// The volesti-facing convex body: Sigma(A) seen through the intrinsic chart.
//
// This is the ONLY place the m <-> n boundary is crossed. volesti never sees
// R^m, because Sigma(A) has measure zero there and a ball walk in R^m would
// accept nothing (see DESIGN.md 6).
//
// Modelled on volesti/include/convex_bodies/ball.h, the minimal example of the
// body concept. BallWalk (uniform_ball_walk.hpp:17-94) touches exactly three
// methods and all must be const, because apply() takes `BallPolytope const& P`.

#ifndef SECPOLY_SECONDARY_BODY_H
#define SECPOLY_SECONDARY_BODY_H

#include <secpoly/separation_oracle.h>
#include <secpoly/types.h>

#include <memory>
#include <utility>

namespace secpoly {

/// \tparam Point a volesti point type (cartesian_geom/point.h)
template <typename Point>
class SecondaryBody {
public:
    using PointType = Point;
    using NTt = typename Point::FT;
    using VT = Eigen::Matrix<NTt, Eigen::Dynamic, 1>;
    using MT = Eigen::Matrix<NTt, Eigen::Dynamic, Eigen::Dynamic>;

    SecondaryBody() = default;

    /// Shared, not owned: is_in() is const while the oracle's caches mutate.
    SecondaryBody(std::shared_ptr<const SeparationOracle> sep, NTt inner_radius)
        : sep_(std::move(sep)), radius_(inner_radius) {}

    /// n = m - d - 1. NEVER m -- this is the whole point of the chart.
    /// Every vector volesti creates is sized by this: uniform_ball_walk.hpp:78
    /// calls GetPointInDsphere(P.dimension(), ...), which builds Point p(dim).
    unsigned int dimension() const {
        return static_cast<unsigned int>(sep_->n());
    }

    /// (centre, radius) in the chart. The chart's origin IS x0, the bootstrap
    /// relative-interior point, so the centre is the zero vector.
    /// Only consulted by BallWalk's default constructor to derive a step size;
    /// production runs pass an explicit tuned L instead.
    std::pair<Point, NTt> InnerBall() const {
        return {Point(dimension()), radius_};
    }

    /// volesti convention: -1 means INSIDE, 0 means outside.
    int is_in(Point const& p, NTt /*tol*/ = NTt(0)) const {
        VecX y(sep_->n());
        for (int i = 0; i < sep_->n(); ++i) y(i) = p[i];

        const SeparationResult r = sep_->separate(y);
        ++queries_;

        // Ambiguous counts as INSIDE, and the asymmetry is deliberate.
        //
        // OUTSIDE is the verdict that carries a proof: it is returned only when
        // the oracle exhibits a hyperplane dominating all of Sigma(A) that y
        // violates. Ambiguous means no such proof was found within the iteration
        // budget -- which, for a point the box test and the cut cache have already
        // passed, is evidence of being inside, not outside.
        //
        // Treating Ambiguous as outside instead makes the sampler NON-UNIFORM:
        // the slow-to-certify points are precisely those near the boundary, so
        // rejecting them carves a systematic bias out of the body. That is
        // measurable -- it showed up as chi-square ~4000 against an independent
        // reference sampler, versus ~30 here.
        //
        // The cost is that a point in the tol-thin shell just OUTSIDE the body may
        // be accepted. That is the GLS weak-membership slack the draft's Theorem 8
        // already assumes, and the end-to-end test asserts containment stays at
        // essentially 100%.
        if (r.verdict != Verdict::Outside) { ++accepted_; return -1; }
        return 0;
    }

    /// Map an intrinsic sample back to a GKZ vector in R^m.
    VecX to_gkz(Point const& p) const {
        VecX y(sep_->n());
        for (int i = 0; i < sep_->n(); ++i) y(i) = p[i];
        return sep_->hull().to_ambient(y);
    }

    const SeparationOracle& separation() const { return *sep_; }

    // Free instrumentation: is_in is called exactly once per ball-walk step, so
    // this is the acceptance rate.
    std::size_t queries() const { return queries_; }
    std::size_t accepted() const { return accepted_; }
    double acceptance_rate() const {
        return queries_ ? static_cast<double>(accepted_) /
                              static_cast<double>(queries_)
                        : 0.0;
    }
    void reset_counters() const { queries_ = accepted_ = 0; }

private:
    std::shared_ptr<const SeparationOracle> sep_;
    NTt radius_ = NTt(1);
    mutable std::size_t queries_ = 0;
    mutable std::size_t accepted_ = 0;
};

}  // namespace secpoly

#endif  // SECPOLY_SECONDARY_BODY_H
