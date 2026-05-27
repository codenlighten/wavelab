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

// ============================================================
// 3D (Phase 5)
// ============================================================

inline Real energy_density_3d(Field<Real, 3> const& u_curr,
                              Field<Real, 3> const& u_prev,
                              Real c, Real dt,
                              Index i, Index j, Index k) noexcept {
    auto const& g = u_curr.grid();
    Real const dx = g.spacing[0];
    Real const dy = g.spacing[1];
    Real const dz = g.spacing[2];
    Index const nx = g.shape[0];
    Index const ny = g.shape[1];
    Index const nz = g.shape[2];

    Real const du_dt = (u_curr(i, j, k) - u_prev(i, j, k)) / dt;
    Real du_dx = Real{0}, du_dy = Real{0}, du_dz = Real{0};
    if (i > 0 && i < nx - 1)
        du_dx = (u_curr(i + 1, j, k) - u_curr(i - 1, j, k)) / (Real{2} * dx);
    if (j > 0 && j < ny - 1)
        du_dy = (u_curr(i, j + 1, k) - u_curr(i, j - 1, k)) / (Real{2} * dy);
    if (k > 0 && k < nz - 1)
        du_dz = (u_curr(i, j, k + 1) - u_curr(i, j, k - 1)) / (Real{2} * dz);

    return Real{0.5} * du_dt * du_dt
         + Real{0.5} * c * c * (du_dx * du_dx + du_dy * du_dy + du_dz * du_dz);
}

inline Real total_energy_3d(Field<Real, 3> const& u_curr,
                            Field<Real, 3> const& u_prev,
                            Real c, Real dt) noexcept {
    auto const& g = u_curr.grid();
    Index const nx = g.shape[0];
    Index const ny = g.shape[1];
    Index const nz = g.shape[2];
    Real sum = Real{0};
    #pragma omp parallel for reduction(+:sum) collapse(3) schedule(static)
    for (Index i = 0; i < nx; ++i) {
        for (Index j = 0; j < ny; ++j) {
            for (Index k = 0; k < nz; ++k) {
                sum += energy_density_3d(u_curr, u_prev, c, dt, i, j, k);
            }
        }
    }
    return sum;
}

inline Real region_energy_3d(Field<Real, 3> const& u_curr,
                             Field<Real, 3> const& u_prev,
                             Real c, Real dt,
                             IVec<3> lo, IVec<3> hi) noexcept {
    Real sum = Real{0};
    #pragma omp parallel for reduction(+:sum) collapse(3) schedule(static)
    for (Index i = lo[0]; i <= hi[0]; ++i) {
        for (Index j = lo[1]; j <= hi[1]; ++j) {
            for (Index k = lo[2]; k <= hi[2]; ++k) {
                sum += energy_density_3d(u_curr, u_prev, c, dt, i, j, k);
            }
        }
    }
    return sum;
}

inline void energy_density_field_3d(Field<Real, 3>&       out,
                                    Field<Real, 3> const& u_curr,
                                    Field<Real, 3> const& u_prev,
                                    Real c, Real dt) noexcept {
    auto const& g = u_curr.grid();
    Index const nx = g.shape[0];
    Index const ny = g.shape[1];
    Index const nz = g.shape[2];
    #pragma omp parallel for collapse(3) schedule(static)
    for (Index i = 0; i < nx; ++i) {
        for (Index j = 0; j < ny; ++j) {
            for (Index k = 0; k < nz; ++k) {
                out(i, j, k) = energy_density_3d(u_curr, u_prev, c, dt, i, j, k);
            }
        }
    }
}

// ============================================================
// Dim-generic dispatchers
// ============================================================
// The _1d / _2d / _3d implementations above ship the actual stencils;
// these templates let callers write `total_energy<D>(...)` instead of
// branching themselves. Used by templated wavecli / Python bindings.

template <int D>
inline Real total_energy(Field<Real, D> const& u_curr,
                         Field<Real, D> const& u_prev,
                         Real c, Real dt) noexcept {
    if constexpr (D == 1) return total_energy_1d(u_curr, u_prev, c, dt);
    else if constexpr (D == 2) return total_energy_2d(u_curr, u_prev, c, dt);
    else if constexpr (D == 3) return total_energy_3d(u_curr, u_prev, c, dt);
}

template <int D>
inline Real region_energy(Field<Real, D> const& u_curr,
                          Field<Real, D> const& u_prev,
                          Real c, Real dt,
                          IVec<D> lo, IVec<D> hi) noexcept {
    if constexpr (D == 2) return region_energy_2d(u_curr, u_prev, c, dt, lo, hi);
    else if constexpr (D == 3) return region_energy_3d(u_curr, u_prev, c, dt, lo, hi);
    else { static_assert(D == 2 || D == 3, "region_energy: D=1 not provided"); return Real{0}; }
}

template <int D>
inline void energy_density_field(Field<Real, D>&       out,
                                 Field<Real, D> const& u_curr,
                                 Field<Real, D> const& u_prev,
                                 Real c, Real dt) noexcept {
    if constexpr (D == 1) energy_density_field_1d(out, u_curr, u_prev, c, dt);
    else if constexpr (D == 2) energy_density_field_2d(out, u_curr, u_prev, c, dt);
    else if constexpr (D == 3) energy_density_field_3d(out, u_curr, u_prev, c, dt);
}

} // namespace wavelab
