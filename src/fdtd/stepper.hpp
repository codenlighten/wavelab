#pragma once
//
// Stepper<D>: abstract interface for any time-stepper that advances a
// scalar field on a D-dimensional grid. Concrete implementations:
//
//   FdtdCpuOmp<D>  — CPU + OpenMP FDTD (Phase 1+)
//   FdtdCuda<D>    — CUDA GPU backend  (Phase 7, optional)
//
// Callers (GUI viewer, scoring code, test fixtures) interact with the
// engine through this interface so backends can be swapped without code
// changes.
//

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"
#include "medium/boundary.hpp"
#include "source/source.hpp"

#include <memory>
#include <vector>

namespace wavelab {

template <int D>
class Stepper {
public:
    virtual ~Stepper() = default;

    // Advance one timestep.
    virtual void step() = 0;

    // Advance N timesteps (default: just loop). Backends may override
    // for batch optimization (e.g. fused kernels).
    virtual void run(Index steps) {
        for (Index i = 0; i < steps; ++i) step();
    }

    // Reset to t = 0 with zeroed fields (sources/boundaries are kept).
    virtual void reset() = 0;

    // --- Query state ---------------------------------------------------

    virtual Grid<D> const&        grid()        const noexcept = 0;
    virtual Field<Real, D> const& current()     const noexcept = 0; // u^t
    virtual Field<Real, D> const& previous()    const noexcept = 0; // u^{t-1}
    virtual Real                  time()        const noexcept = 0;
    virtual Real                  dt()          const noexcept = 0;
    virtual Index                 step_count()  const noexcept = 0;

    // --- Mutate the scene ---------------------------------------------

    virtual void add_source(std::shared_ptr<Source<D>> s) = 0;
    virtual void set_boundary(std::shared_ptr<BoundaryCondition<D>> b) = 0;
    virtual void clear_sources() = 0;
};

} // namespace wavelab
