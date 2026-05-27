#pragma once
//
// HarmonicSource<D>: clean sinusoidal point source (overview §4).
//
//     S(t) = A * sin(2π f t + φ)
//
// Injection mode is a *hard* set (u_next at the source cell is overwritten),
// which is the simplest and most predictable for analytical validation
// (plane-wave propagation, dispersion checks). Soft additive injection is
// available via the `soft` ctor flag for energy-balance-style experiments.
//

#include "source.hpp"
#include "core/types.hpp"

#include <cmath>
#include <numbers>

namespace wavelab {

template <int D>
class HarmonicSource : public Source<D> {
public:
    HarmonicSource(IVec<D> loc,
                   Real    amplitude,
                   Real    frequency,
                   Real    phase = Real{0},
                   bool    soft  = false) noexcept
        : loc_(loc), amp_(amplitude), freq_(frequency), phase_(phase), soft_(soft) {}

    void inject(Field<Real, D>& u_next, Real t, Real dt) const override {
        Real const omega = static_cast<Real>(2.0 * std::numbers::pi) * freq_;
        Real const value = amp_ * std::sin(omega * t + phase_);
        if (soft_) {
            // Soft source: additive, scaled by dt² so the discrete equation
            // u^{t+1} = 2u^t - u^{t-1} + λ²∇²u + dt²·S^t holds.
            u_next(loc_) += dt * dt * value;
        } else {
            u_next(loc_) = value;
        }
    }

    IVec<D> location()  const noexcept { return loc_; }
    Real    amplitude() const noexcept { return amp_; }
    Real    frequency() const noexcept { return freq_; }
    Real    phase()     const noexcept { return phase_; }
    bool    is_soft()   const noexcept { return soft_; }

    void set_amplitude(Real a) noexcept { amp_   = a; }
    void set_frequency(Real f) noexcept { freq_  = f; }
    void set_phase(Real p)     noexcept { phase_ = p; }

private:
    IVec<D> loc_;
    Real    amp_;
    Real    freq_;
    Real    phase_;
    bool    soft_;
};

} // namespace wavelab
