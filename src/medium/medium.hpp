#pragma once
//
// Medium<D>: bundles the spatially-varying material properties seen by
// the FDTD stepper.
//
//   c(x)     — local wave speed (overview §9, set from refractive index)
//   alpha(x) — local damping (overview §5, §8; used by PML region too)
//
// Future additions (Phase 4+): index_field (n), permittivity, anisotropy.
// Atom-derived fields (rho, charge, hydro) live separately in
// molecule/field_builder; they feed *into* Medium via medium/build_medium.
//

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"

#include <algorithm>

namespace wavelab {

template <int D>
struct Medium {
    Field<Real, D> c;       // wave speed per cell
    Field<Real, D> alpha;   // damping coefficient per cell

    // Construct a uniform medium: c(x) = c0 everywhere, alpha(x) = a0.
    static Medium uniform(Grid<D> const& grid, Real c0, Real a0 = Real{0}) {
        Medium m{
            Field<Real, D>(grid, c0),
            Field<Real, D>(grid, a0)
        };
        return m;
    }

    // Maximum local wave speed — drives the CFL bound.
    Real max_c() const noexcept {
        Real const* p = c.data();
        Index const  n = c.size();
        if (n == 0) return Real{0};
        Real m = p[0];
        for (Index i = 1; i < n; ++i) {
            if (p[i] > m) m = p[i];
        }
        return m;
    }

    Grid<D> const& grid() const noexcept { return c.grid(); }
};

} // namespace wavelab
