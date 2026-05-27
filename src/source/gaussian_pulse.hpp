#pragma once
//
// GaussianPulse<D>: short broad-spectrum probe (overview §4).
//
//     S(t) = A * exp(-(t - t0)² / (2 σ_t²)) * sin(2π f t)
//
// The envelope makes this an excellent broadband impulse — the time-width
// σ_t inversely controls the frequency span (Δf ~ 1/(2π σ_t)). Good for
// excitation that probes multiple spatial scales in one shot.
//

#include "source.hpp"
#include "core/types.hpp"

#include <cmath>
#include <numbers>

namespace wavelab {

template <int D>
class GaussianPulse : public Source<D> {
public:
    GaussianPulse(IVec<D> loc,
                  Real    amplitude,
                  Real    frequency,
                  Real    t0,
                  Real    sigma_t,
                  bool    soft = false) noexcept
        : loc_(loc), amp_(amplitude), freq_(frequency),
          t0_(t0), sigma_(sigma_t), soft_(soft) {}

    void inject(Field<Real, D>& u_next, Real t, Real dt) const override {
        Real const omega = static_cast<Real>(2.0 * std::numbers::pi) * freq_;
        Real const dt_   = (t - t0_) / sigma_;
        Real const env   = std::exp(static_cast<Real>(-0.5) * dt_ * dt_);
        Real const value = amp_ * env * std::sin(omega * t);
        if (soft_) {
            u_next(loc_) += dt * dt * value;
        } else {
            u_next(loc_) = value;
        }
    }

    // After 6σ the envelope is < 1e-8 of peak — safe to skip.
    bool active(Real t) const noexcept override {
        return std::abs(t - t0_) < static_cast<Real>(6) * sigma_;
    }

    IVec<D> location() const noexcept { return loc_; }

private:
    IVec<D> loc_;
    Real    amp_;
    Real    freq_;
    Real    t0_;
    Real    sigma_;
    bool    soft_;
};

} // namespace wavelab
