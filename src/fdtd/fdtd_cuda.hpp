#pragma once
//
// FdtdCuda<D>: GPU FDTD stepper behind the existing Stepper<D> interface.
//
// PIMPL'd so consuming TUs never see CUDA headers (they don't need nvcc).
// The first cut targets D=2 with Dirichlet boundaries; D=3 follows the
// same pattern. Sources are handled host-side via lazy device⇄host
// synchronization — fine for the typical "1 point source" case; multi-
// source / soft-source GPU paths land later if benchmarks demand them.
//
// Memory layout matches FdtdCpuOmp: row-major, three rotating Real
// buffers (u_prev / u_curr / u_next) plus Medium fields (c, α). All
// device-resident; host mirrors of u_curr/u_prev are populated lazily
// when `current()` / `previous()` is called.
//

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"
#include "fdtd/stepper.hpp"
#include "medium/boundary.hpp"
#include "medium/medium.hpp"
#include "source/source.hpp"

#include <memory>
#include <vector>

namespace wavelab {

template <int D>
class FdtdCuda final : public Stepper<D> {
public:
    // Heterogeneous-medium constructor.
    FdtdCuda(Grid<D> grid, Medium<D> medium, Real dt);
    // Uniform-scalar shortcut (mirrors FdtdCpuOmp).
    FdtdCuda(Grid<D> grid, Real c, Real dt, Real gamma = Real{0});

    ~FdtdCuda();
    FdtdCuda(FdtdCuda const&)            = delete;
    FdtdCuda& operator=(FdtdCuda const&) = delete;

    // --- Stepper<D> ----------------------------------------------------
    void step()  override;
    void reset() override;

    Grid<D> const&        grid()       const noexcept override;
    Field<Real, D> const& current()    const noexcept override;
    Field<Real, D> const& previous()   const noexcept override;
    Real                  time()       const noexcept override;
    Real                  dt()         const noexcept override;
    Index                 step_count() const noexcept override;

    void add_source(std::shared_ptr<Source<D>> s)         override;
    void set_boundary(std::shared_ptr<BoundaryCondition<D>> b) override;
    void clear_sources()                                   override;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

extern template class FdtdCuda<2>;

} // namespace wavelab
