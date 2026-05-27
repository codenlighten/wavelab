#pragma once
//
// FdtdCpuOmp<D>: CPU + OpenMP finite-difference time-domain wave stepper.
//
// Discrete update (overview §2, §5, §8) with spatially varying medium:
//
//     u^{t+1}[x] = (2 - α(x)·Δt) · u^t[x]
//                - (1 - α(x)·Δt) · u^{t-1}[x]
//                + (c(x)·Δt/Δx)² · ( Σ neighbors - 2D · u^t[x] )
//
// Constant-medium and heterogeneous-medium cases share the same stencil
// — the medium fields just read scalar c0/α0 from a uniform Medium for
// the constant case. Phase 1's scalar (c, gamma) constructor stays as a
// convenience overload.
//
// Storage rotation: three Field<Real,D> swapped each step
//   u_prev (u^{t-1}) -> read
//   u_curr (u^t)     -> read
//   u_next (u^{t+1}) -> written, then rotates into u_curr
//
// Coverage: D=1 (Phase 1), D=2 (Phase 2); D=3 lands in Phase 5.
//

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"
#include "fdtd/stability.hpp"
#include "fdtd/stepper.hpp"
#include "medium/boundary.hpp"
#include "medium/medium.hpp"
#include "source/source.hpp"

#include <memory>
#include <stdexcept>
#include <vector>

namespace wavelab {

template <int D>
class FdtdCpuOmp final : public Stepper<D> {
public:
    // Heterogeneous-medium constructor (Phase 2+).
    FdtdCpuOmp(Grid<D> grid, Medium<D> medium, Real dt)
        : grid_(grid),
          medium_(std::move(medium)),
          u_prev_(grid, Real{0}),
          u_curr_(grid, Real{0}),
          u_next_(grid, Real{0}),
          dt_(dt) {
        if (medium_.grid().num_cells() != grid_.num_cells()) {
            throw std::invalid_argument(
                "FdtdCpuOmp: medium grid does not match stepper grid");
        }
        Real const cmax = medium_.max_c();
        if (!cfl_satisfied<D>(grid_, cmax, dt_)) {
            throw std::invalid_argument(
                "FdtdCpuOmp: dt violates CFL bound (Courant condition §3)");
        }
        install_default_boundary_();
    }

    // Uniform-scalar convenience (Phase 1 compatibility).
    FdtdCpuOmp(Grid<D> grid, Real c, Real dt, Real gamma = Real{0})
        : FdtdCpuOmp(grid, Medium<D>::uniform(grid, c, gamma), dt) {}

    // --- Stepper<D> interface -----------------------------------------

    void step() override {
        advance_stencil_();
        if (boundary_) boundary_->apply(u_next_, u_curr_, max_c_(), dt_);
        inject_sources_(t_ + dt_);
        rotate_fields_();
        t_ += dt_;
        ++step_;
    }

    void reset() override {
        u_prev_.zero();
        u_curr_.zero();
        u_next_.zero();
        t_     = Real{0};
        step_  = 0;
    }

    Grid<D> const&        grid()       const noexcept override { return grid_; }
    Field<Real, D> const& current()    const noexcept override { return u_curr_; }
    Field<Real, D> const& previous()   const noexcept override { return u_prev_; }
    Real                  time()       const noexcept override { return t_; }
    Real                  dt()         const noexcept override { return dt_; }
    Index                 step_count() const noexcept override { return step_; }

    void add_source(std::shared_ptr<Source<D>> s) override {
        sources_.push_back(std::move(s));
    }
    void set_boundary(std::shared_ptr<BoundaryCondition<D>> b) override {
        boundary_ = std::move(b);
    }
    void clear_sources() override { sources_.clear(); }

    // --- Test / IC / inspection helpers (not in Stepper<D>) -----------

    Field<Real, D>&   current_mut()  noexcept { return u_curr_; }
    Field<Real, D>&   previous_mut() noexcept { return u_prev_; }
    Medium<D> const&  medium()       const noexcept { return medium_; }
    Medium<D>&        medium_mut()   noexcept { return medium_; }
    // Effective wave speed for diagnostic display (max c in the domain).
    Real              wave_speed()   const noexcept { return max_c_(); }

private:
    Real max_c_() const noexcept { return medium_.max_c(); }

    void install_default_boundary_();   // out-of-line per D specialization

    // -----------------------------------------------------------------
    // Discrete stencil — overwrites u_next_ from (u_curr_, u_prev_).
    // -----------------------------------------------------------------
    void advance_stencil_() {
        Real const dx     = grid_.spacing[0];
        Real const inv_dx = Real{1} / dx;

        Real const* up = u_prev_.data();
        Real const* uc = u_curr_.data();
        Real* un       = u_next_.data();
        Real const* cf = medium_.c.data();
        Real const* af = medium_.alpha.data();

        if constexpr (D == 1) {
            Index const n = grid_.shape[0];
            #pragma omp parallel for schedule(static)
            for (Index i = 1; i < n - 1; ++i) {
                Real const c     = cf[i];
                Real const alpha = af[i];
                Real const lam   = c * dt_ * inv_dx;
                Real const lam2  = lam * lam;
                Real const a     = Real{2} - alpha * dt_;
                Real const b     = Real{1} - alpha * dt_;
                Real const lap   = uc[i + 1] + uc[i - 1] - Real{2} * uc[i];
                un[i] = a * uc[i] - b * up[i] + lam2 * lap;
            }
            un[0]      = Real{0};
            un[n - 1]  = Real{0};
        } else if constexpr (D == 2) {
            Index const nx = grid_.shape[0];
            Index const ny = grid_.shape[1];

            #pragma omp parallel for collapse(2) schedule(static)
            for (Index i = 1; i < nx - 1; ++i) {
                for (Index j = 1; j < ny - 1; ++j) {
                    Index const idx  = i * ny + j;
                    Real const c     = cf[idx];
                    Real const alpha = af[idx];
                    Real const lam   = c * dt_ * inv_dx;
                    Real const lam2  = lam * lam;
                    Real const a     = Real{2} - alpha * dt_;
                    Real const b     = Real{1} - alpha * dt_;
                    Real const lap   = uc[idx + ny] + uc[idx - ny]
                                     + uc[idx + 1]  + uc[idx - 1]
                                     - Real{4} * uc[idx];
                    un[idx] = a * uc[idx] - b * up[idx] + lam2 * lap;
                }
            }
            // Boundary rows/cols get zero by default; BC overwrites.
            #pragma omp parallel for schedule(static)
            for (Index j = 0; j < ny; ++j) {
                un[0 * ny + j]        = Real{0};
                un[(nx - 1) * ny + j] = Real{0};
            }
            #pragma omp parallel for schedule(static)
            for (Index i = 0; i < nx; ++i) {
                un[i * ny + 0]        = Real{0};
                un[i * ny + (ny - 1)] = Real{0};
            }
        } else {
            static_assert(D == 1 || D == 2,
                "FdtdCpuOmp<D>: D=3 implemented in Phase 5");
        }
    }

    void inject_sources_(Real t_new) {
        for (auto const& s : sources_) {
            if (s->active(t_new)) s->inject(u_next_, t_new, dt_);
        }
    }

    void rotate_fields_() noexcept {
        u_prev_.swap(u_curr_);
        u_curr_.swap(u_next_);
    }

private:
    Grid<D>        grid_;
    Medium<D>      medium_;
    Field<Real, D> u_prev_;
    Field<Real, D> u_curr_;
    Field<Real, D> u_next_;
    Real           dt_;

    std::vector<std::shared_ptr<Source<D>>>     sources_;
    std::shared_ptr<BoundaryCondition<D>>       boundary_;

    Real  t_    = Real{0};
    Index step_ = 0;
};

// Out-of-line default-BC specializations live in boundary.hpp's section
// below; we just forward-call install_default_boundary_() above.

template <>
inline void FdtdCpuOmp<1>::install_default_boundary_() {
    boundary_ = std::make_shared<Dirichlet1D>();
}

template <>
inline void FdtdCpuOmp<2>::install_default_boundary_() {
    boundary_ = std::make_shared<Dirichlet2D>();
}

} // namespace wavelab
