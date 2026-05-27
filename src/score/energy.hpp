#pragma once
//
// Energy diagnostics (overview §14, §15).
//
//   E(x,t) = ½(∂u/∂t)² + ½c²|∇u|²
//
// Discrete forms:
//   * ∂u/∂t ≈ (u^t - u^{t-1}) / Δt
//   * ∇u    ≈ centered finite differences on u^t (zero at boundary cells)
//
// Two flavors:
//   * total_energy        — physically motivated centered-form sum (§14).
//                            Oscillates ±O(Δt²) per step under leapfrog.
//   * total_energy_conserved — leapfrog-staggered discrete invariant.
//                            Stays flat to float precision; the right
//                            diagnostic for "is the integrator
//                            conservative". Uses the inner product of
//                            ∇u^t and ∇u^{t-1}.
//
// Probe-region overloads sum over a rectangular [lo,hi] subdomain — used
// by reflection/transmission scoring (§17) and scattering loss (§16).
//

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"

#include <algorithm>
#include <cmath>

namespace wavelab {

// ============================================================
// 1D (Phase 1 surface)
// ============================================================

inline Real energy_density_1d(Field<Real, 1> const& u_curr,
                              Field<Real, 1> const& u_prev,
                              Real c, Real dt, Index i) noexcept {
    auto const& g = u_curr.grid();
    Real const dx = g.spacing[0];
    Index const n = g.shape[0];

    Real const du_dt = (u_curr(i) - u_prev(i)) / dt;
    Real du_dx = Real{0};
    if (i > 0 && i < n - 1) {
        du_dx = (u_curr(i + 1) - u_curr(i - 1)) / (Real{2} * dx);
    }
    return Real{0.5} * du_dt * du_dt
         + Real{0.5} * c * c * du_dx * du_dx;
}

inline Real total_energy_1d(Field<Real, 1> const& u_curr,
                            Field<Real, 1> const& u_prev,
                            Real c, Real dt) noexcept {
    auto const& g = u_curr.grid();
    Index const n = g.shape[0];
    Real sum = Real{0};
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (Index i = 0; i < n; ++i) {
        sum += energy_density_1d(u_curr, u_prev, c, dt, i);
    }
    return sum;
}

inline Real total_energy_conserved_1d(Field<Real, 1> const& u_curr,
                                      Field<Real, 1> const& u_prev,
                                      Real c, Real dt) noexcept {
    auto const& g  = u_curr.grid();
    Index const n  = g.shape[0];
    Real const  dx = g.spacing[0];

    Real sum_v = Real{0};
    Real sum_g = Real{0};
    #pragma omp parallel for reduction(+:sum_v) schedule(static)
    for (Index i = 0; i < n; ++i) {
        Real const v = (u_curr(i) - u_prev(i)) / dt;
        sum_v += v * v;
    }
    #pragma omp parallel for reduction(+:sum_g) schedule(static)
    for (Index i = 0; i < n - 1; ++i) {
        Real const gt = (u_curr(i + 1) - u_curr(i)) / dx;
        Real const gp = (u_prev(i + 1) - u_prev(i)) / dx;
        sum_g += gt * gp;
    }
    return Real{0.5} * sum_v + Real{0.5} * c * c * sum_g;
}

inline void energy_density_field_1d(Field<Real, 1>&       out,
                                    Field<Real, 1> const& u_curr,
                                    Field<Real, 1> const& u_prev,
                                    Real c, Real dt) noexcept {
    Index const n = u_curr.grid().shape[0];
    #pragma omp parallel for schedule(static)
    for (Index i = 0; i < n; ++i) {
        out(i) = energy_density_1d(u_curr, u_prev, c, dt, i);
    }
}

// ============================================================
// 2D (Phase 2)
// ============================================================

inline Real energy_density_2d(Field<Real, 2> const& u_curr,
                              Field<Real, 2> const& u_prev,
                              Real c, Real dt,
                              Index i, Index j) noexcept {
    auto const& g = u_curr.grid();
    Real const dx = g.spacing[0];
    Real const dy = g.spacing[1];
    Index const nx = g.shape[0];
    Index const ny = g.shape[1];

    Real const du_dt = (u_curr(i, j) - u_prev(i, j)) / dt;
    Real du_dx = Real{0};
    Real du_dy = Real{0};
    if (i > 0 && i < nx - 1)
        du_dx = (u_curr(i + 1, j) - u_curr(i - 1, j)) / (Real{2} * dx);
    if (j > 0 && j < ny - 1)
        du_dy = (u_curr(i, j + 1) - u_curr(i, j - 1)) / (Real{2} * dy);

    return Real{0.5} * du_dt * du_dt
         + Real{0.5} * c * c * (du_dx * du_dx + du_dy * du_dy);
}

inline Real total_energy_2d(Field<Real, 2> const& u_curr,
                            Field<Real, 2> const& u_prev,
                            Real c, Real dt) noexcept {
    auto const& g = u_curr.grid();
    Index const nx = g.shape[0];
    Index const ny = g.shape[1];
    Real sum = Real{0};
    #pragma omp parallel for reduction(+:sum) collapse(2) schedule(static)
    for (Index i = 0; i < nx; ++i) {
        for (Index j = 0; j < ny; ++j) {
            sum += energy_density_2d(u_curr, u_prev, c, dt, i, j);
        }
    }
    return sum;
}

// Energy over a rectangular probe region [lo, hi] (inclusive).
inline Real region_energy_2d(Field<Real, 2> const& u_curr,
                             Field<Real, 2> const& u_prev,
                             Real c, Real dt,
                             IVec<2> lo, IVec<2> hi) noexcept {
    Real sum = Real{0};
    #pragma omp parallel for reduction(+:sum) collapse(2) schedule(static)
    for (Index i = lo[0]; i <= hi[0]; ++i) {
        for (Index j = lo[1]; j <= hi[1]; ++j) {
            sum += energy_density_2d(u_curr, u_prev, c, dt, i, j);
        }
    }
    return sum;
}

inline void energy_density_field_2d(Field<Real, 2>&       out,
                                    Field<Real, 2> const& u_curr,
                                    Field<Real, 2> const& u_prev,
                                    Real c, Real dt) noexcept {
    auto const& g = u_curr.grid();
    Index const nx = g.shape[0];
    Index const ny = g.shape[1];
    #pragma omp parallel for collapse(2) schedule(static)
    for (Index i = 0; i < nx; ++i) {
        for (Index j = 0; j < ny; ++j) {
            out(i, j) = energy_density_2d(u_curr, u_prev, c, dt, i, j);
        }
    }
}

} // namespace wavelab
