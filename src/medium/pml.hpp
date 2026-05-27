#pragma once
//
// Perfectly Matched Layer (graded-damping variant) for 2D / 3D FDTD.
//
// Instead of Berenger's field-splitting PML, we use the simpler
// "sponge-layer" approximation: a region of `cells` width around each
// boundary face where the medium damping α(x) ramps polynomially from
// 0 (interior edge) to α_max (outer wall):
//
//     α(d) = α_max · (d / L_pml)^m       (m = 3 by default)
//
// With Dirichlet at the outer wall, single-pass attenuation A through
// the layer is exp(-α_max·L·dt/(m+1)·c⁻¹). Reflection from the
// Dirichlet wall arrives back with amplitude A², so for ε reflection
// target ε ≤ 0.5%, pick α_max so A² ≤ ε.
//
// Default α_max ~ 4·c_max/dx targets sub-0.5% reflection with 20 cells
// of PML and m = 3 under a CFL = 0.5 timestep. Tune via constructor.
//

#include "core/types.hpp"
#include "medium/medium.hpp"

#include <algorithm>
#include <cmath>

namespace wavelab {

struct PmlSpec {
    Index cells              = 20;        // depth of PML region per face
    Real  alpha_max_factor   = 4.0_r;     // α_max = factor · c_max / dx
    int   polynomial_order   = 3;         // profile exponent m
};

// Stamp a PML region into `medium.alpha`. Idempotent only in the sense
// that calling this with the same spec twice gives the same alpha field;
// it OVERWRITES any prior alpha values within the PML frame.
template <int D>
inline void apply_pml(Medium<D>& medium, PmlSpec const& spec = {}) {
    auto const& g  = medium.alpha.grid();
    Real const  dx = g.spacing[0];          // assume uniform spacing
    Real const  cmax = medium.max_c();
    Real const  alpha_max = spec.alpha_max_factor * cmax / dx;
    Index const L = spec.cells;
    int const   m = spec.polynomial_order;

    auto profile = [=](Index d_into_pml) -> Real {
        // d=L is the inner edge of the PML (no damping); d=0 is outer wall.
        Real const t = static_cast<Real>(L - d_into_pml) / static_cast<Real>(L);
        Real const tm = std::pow(t, static_cast<Real>(m));
        return alpha_max * tm;
    };

    if constexpr (D == 1) {
        Index const n = g.shape[0];
        Real* a = medium.alpha.data();
        for (Index i = 0; i < n; ++i) {
            Index const d_left  = i;
            Index const d_right = n - 1 - i;
            Index const d = std::min(d_left, d_right);
            if (d < L) a[i] = std::max(a[i], profile(d));
        }
    } else if constexpr (D == 2) {
        Index const nx = g.shape[0];
        Index const ny = g.shape[1];
        Real* a = medium.alpha.data();
        for (Index i = 0; i < nx; ++i) {
            for (Index j = 0; j < ny; ++j) {
                Index const dx_left  = i;
                Index const dx_right = nx - 1 - i;
                Index const dy_bot   = j;
                Index const dy_top   = ny - 1 - j;
                Index const d = std::min({dx_left, dx_right, dy_bot, dy_top});
                if (d < L) {
                    Index const idx = i * ny + j;
                    a[idx] = std::max(a[idx], profile(d));
                }
            }
        }
    } else {
        static_assert(D == 1 || D == 2,
            "apply_pml: D=3 lands in Phase 5");
    }
}

} // namespace wavelab
