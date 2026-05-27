#pragma once
//
// Scattering / reflection / transmission scoring (overview §16, §17).
//
// The pattern: pick two probe regions in the 2D domain — one "incident"
// region between the source and the scatterer, one "transmitted" region
// past the scatterer. Time-integrate the energy passing through each
// region during a controlled experiment.
//
// For pulsed experiments, the typical workflow is:
//   1. Build the scene (empty domain vs scene-with-scatterer).
//   2. Inject a short pulse.
//   3. Track:
//        E_incident   — peak energy in the "incident" region BEFORE the
//                       wavefront reaches the scatterer
//        E_reflected  — peak energy back in the "incident" region AFTER
//                       reflection
//        E_transmitted — peak energy in the "transmitted" region
//
//   4. R = E_reflected / E_incident
//      T = E_transmitted / E_incident
//      A = 1 - R - T              (absorption)
//
// The `ProbeRegion` helper just bundles bounds + reuses region_energy_2d
// from energy.hpp. Coupling to an actual experiment loop lives in the
// test code and (later) the GUI / scene runner.
//

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"
#include "score/energy.hpp"

#include <type_traits>

namespace wavelab {

struct ProbeRegion2D {
    IVec<2> lo{};   // inclusive lower indices
    IVec<2> hi{};   // inclusive upper indices

    Real energy(Field<Real, 2> const& u_curr,
                Field<Real, 2> const& u_prev,
                Real c, Real dt) const noexcept {
        return region_energy_2d(u_curr, u_prev, c, dt, lo, hi);
    }
};

struct ProbeRegion3D {
    IVec<3> lo{};   // inclusive lower indices
    IVec<3> hi{};   // inclusive upper indices

    Real energy(Field<Real, 3> const& u_curr,
                Field<Real, 3> const& u_prev,
                Real c, Real dt) const noexcept {
        return region_energy_3d(u_curr, u_prev, c, dt, lo, hi);
    }
};

// Dim-generic factory: build a probe region from a center cell + half-
// width (in cells), clamped to the grid. Used by wavecli's --region
// flag and by auto-region derivation from a ligand bounding box.
template <int D>
inline auto make_probe_region(Grid<D> const& g,
                              IVec<D> center,
                              Index half_width_cells) {
    using Region = std::conditional_t<D == 2, ProbeRegion2D, ProbeRegion3D>;
    Region r{};
    for (std::size_t d = 0; d < static_cast<std::size_t>(D); ++d) {
        Index lo = center[d] - half_width_cells;
        Index hi = center[d] + half_width_cells;
        if (lo < 0)              lo = 0;
        if (hi >= g.shape[d])    hi = g.shape[d] - 1;
        r.lo[d] = lo;
        r.hi[d] = hi;
    }
    return r;
}

// Ratio: E(post) / E(reference). Used for both R_E (§15: pocket+ligand /
// pocket) and scattering loss (§16: 1 - E_out / E_in).
inline Real energy_ratio(Real e_numerator, Real e_denominator) noexcept {
    if (e_denominator == Real{0}) return Real{0};
    return e_numerator / e_denominator;
}

inline Real scattering_loss(Real e_in, Real e_out) noexcept {
    if (e_in == Real{0}) return Real{0};
    return Real{1} - (e_out / e_in);
}

// Reflection / transmission / absorption packet.
struct RTA {
    Real R;   // reflection coefficient
    Real T;   // transmission coefficient
    Real A;   // absorption (= 1 - R - T)
};

inline RTA compute_rta(Real e_incident,
                       Real e_reflected,
                       Real e_transmitted) noexcept {
    RTA out{};
    if (e_incident == Real{0}) return out;
    out.R = e_reflected   / e_incident;
    out.T = e_transmitted / e_incident;
    out.A = Real{1} - out.R - out.T;
    return out;
}

} // namespace wavelab
