#pragma once
//
// Source<D>: abstract excitation that injects energy into a Field<Real,D>
// at each time step. Implementations decide where (location) and how
// (hard override vs additive soft injection).
//
// Concrete sources live in harmonic.hpp, gaussian_pulse.hpp, impulse.hpp.
//

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"

namespace wavelab {

template <int D>
class Source {
public:
    virtual ~Source() = default;

    // Apply this source's contribution at the given simulation time.
    // `u_next` is the freshly computed u^{t+1} (sources fire AFTER the
    // stencil update, BEFORE boundary conditions are applied).
    virtual void inject(Field<Real, D>& u_next,
                        Real             t,
                        Real             dt) const = 0;

    // Whether this source is still active. Stepper may stop calling
    // inject() once this returns false to save work.
    virtual bool active(Real t) const noexcept { (void)t; return true; }
};

} // namespace wavelab
