#pragma once
//
// Assemble the optical Medium<D> from molecular scalar fields per
// overview §9:
//
//     n(x) = n0 + β_ρ ρ(x) + β_q |Q(x)| + β_h H(x)
//     c(x) = c0 / n(x)
//
// Each contribution is optional — pass nullptr to skip that term.
// Use `apply_index_slab` for clean Fresnel tests where you want to
// specify c(x) directly without going through the molecular pipeline.
//

#include "core/field.hpp"
#include "core/types.hpp"
#include "medium/medium.hpp"

#include <algorithm>
#include <cmath>

namespace wavelab {

struct MediumWeights {
    Real n0       = 1.0_r;
    Real beta_rho = 0.5_r;
    Real beta_q   = 0.0_r;
    Real beta_h   = 0.0_r;
    Real n_max    = 5.0_r;   // hard clip so c(x) doesn't underflow
};

// Build n(x) into `n_out`, then derive c(x) = c0 / n(x) into `medium.c`.
// Both fields must share the same grid.
template <int D>
inline void build_medium_from_fields(Medium<D>&            medium,
                                     Field<Real, D> const* rho,
                                     Field<Real, D> const* charge,
                                     Field<Real, D> const* hydro,
                                     Real                  c0,
                                     MediumWeights const&  w = {}) {
    auto const& g = medium.c.grid();
    Real* cf = medium.c.data();
    Index const N = medium.c.size();

    #pragma omp parallel for schedule(static)
    for (Index k = 0; k < N; ++k) {
        Real n = w.n0;
        if (rho)    n += w.beta_rho * rho->data()[k];
        if (charge) n += w.beta_q   * std::abs(charge->data()[k]);
        if (hydro)  n += w.beta_h   * hydro->data()[k];
        if (n > w.n_max) n = w.n_max;
        if (n < 1.0_r)   n = 1.0_r;   // avoid superluminal artifacts
        cf[k] = c0 / n;
    }
    (void)g;
}

// Drop a rectangular block of refractive index n_inside (with c_inside
// = c0 / n_inside) into the medium. Useful for clean Fresnel reflection
// tests where you want a known abrupt impedance step.
template <int D>
inline void apply_index_slab(Medium<D>& medium, Real x_lo, Real x_hi,
                             Real c_inside) {
    auto const& g = medium.c.grid();
    if constexpr (D == 1) {
        Index const n  = g.shape[0];
        Real const  dx = g.spacing[0];
        for (Index i = 0; i < n; ++i) {
            Real const xc = g.origin[0] + (static_cast<Real>(i) + 0.5_r) * dx;
            if (xc >= x_lo && xc <= x_hi) {
                medium.c(i) = c_inside;
            }
        }
    } else if constexpr (D == 2) {
        Index const nx = g.shape[0];
        Index const ny = g.shape[1];
        Real const  dx = g.spacing[0];
        for (Index i = 0; i < nx; ++i) {
            Real const xc = g.origin[0] + (static_cast<Real>(i) + 0.5_r) * dx;
            if (xc >= x_lo && xc <= x_hi) {
                for (Index j = 0; j < ny; ++j) {
                    medium.c(i, j) = c_inside;
                }
            }
        }
    } else if constexpr (D == 3) {
        Index const nx = g.shape[0];
        Index const ny = g.shape[1];
        Index const nz = g.shape[2];
        Real const  dx = g.spacing[0];
        for (Index i = 0; i < nx; ++i) {
            Real const xc = g.origin[0] + (static_cast<Real>(i) + 0.5_r) * dx;
            if (xc >= x_lo && xc <= x_hi) {
                for (Index j = 0; j < ny; ++j) {
                    for (Index k = 0; k < nz; ++k) {
                        medium.c(i, j, k) = c_inside;
                    }
                }
            }
        }
    }
}

} // namespace wavelab
