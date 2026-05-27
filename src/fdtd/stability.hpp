#pragma once
//
// Courant–Friedrichs–Lewy stability helpers (overview §3).
//
//     Δt ≤ Δx / (c √D)
//
// For heterogeneous media (Phase 2+), the bound must be enforced against
// max(c(x)), since the wave speed varies cell-to-cell. Phase 1 assumes
// uniform c.
//

#include "core/grid.hpp"
#include "core/types.hpp"

#include <cmath>

namespace wavelab {

using namespace wavelab::literals;

template <int D>
constexpr Real cfl_dt_max(Grid<D> const& grid, Real c_max) noexcept {
    Real const h = grid.min_spacing();
    Real const root_d = std::sqrt(static_cast<Real>(D));
    return h / (c_max * root_d);
}

// λ = c·Δt / Δx (Courant number, the dimensionless quantity that appears
// in the discrete stencil). Should satisfy λ ≤ 1/√D for stability.
template <int D>
constexpr Real courant_number(Real c, Real dt, Real dx) noexcept {
    (void)c; (void)dt; (void)dx;
    return c * dt / dx;
}

// True if (c, dt, grid) satisfies the CFL bound with a safety margin.
template <int D>
constexpr bool cfl_satisfied(Grid<D> const& grid, Real c_max, Real dt,
                             Real safety = 0.999_r) noexcept {
    return dt <= safety * cfl_dt_max<D>(grid, c_max);
}

} // namespace wavelab
