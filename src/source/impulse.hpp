#pragma once
//
// Impulse<D>: single-step δ(t - t0) excitation. The discrete Green's
// function of the pocket is what you measure after firing one of these
// at the empty domain (overview §33).
//
// We inject for exactly one time step — the step whose simulation time
// crosses t0 — and stay quiet thereafter.
//

#include "source.hpp"
#include "core/types.hpp"

#include <cmath>

namespace wavelab {

template <int D>
class Impulse : public Source<D> {
public:
    Impulse(IVec<D> loc, Real amplitude, Real t0 = Real{0}, bool soft = false) noexcept
        : loc_(loc), amp_(amplitude), t0_(t0), soft_(soft) {}

    void inject(Field<Real, D>& u_next, Real t, Real dt) const override {
        if (std::abs(t - t0_) < static_cast<Real>(0.5) * dt) {
            if (soft_) {
                u_next(loc_) += dt * dt * amp_;
            } else {
                u_next(loc_) = amp_;
            }
            fired_ = true;
        }
    }

    bool active(Real t) const noexcept override {
        if (fired_) return false;
        return t <= t0_ + dt_window();
    }

    IVec<D> location() const noexcept { return loc_; }

private:
    // Liberal window before t0 to ensure we don't deactivate before firing.
    static constexpr Real dt_window() noexcept { return Real{1e6}; }

    IVec<D>      loc_;
    Real         amp_;
    Real         t0_;
    bool         soft_;
    mutable bool fired_ = false;
};

} // namespace wavelab
