#pragma once
//
// Energy-distribution entropy + concentration scores (overview §23, §24, §25).
//
//   p_x      = E(x) / Σ E(x)
//   H_E      = -Σ p_x log p_x                    Shannon entropy of energy
//   S_focus  = 1 - H_E / log(M)                  normalized focus score
//   C_hotspot = E(Ω_hot) / E(Ω_full)             energy fraction in pocket region
//   O_ligand  = E_empty(Ω_L) / E_empty(Ω_full)   importance of ligand-occupied region
//
// All operate on energy-density fields (built from `energy_density_field_2d`
// etc. in energy.hpp). They are dim-generic in spirit; the implementations
// here are 2D-first to match the active engine path.
//

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"
#include "score/energy.hpp"

#include <algorithm>
#include <cmath>

namespace wavelab {

// ----------------------------------------------------------------------
// Shannon entropy of normalized energy distribution
// ----------------------------------------------------------------------
inline Real energy_entropy(Field<Real, 2> const& e) noexcept {
    Index const n = e.size();
    Real total = 0.0_r;
    Real const* p = e.data();

    #pragma omp parallel for reduction(+:total) schedule(static)
    for (Index k = 0; k < n; ++k) total += p[k];
    if (total <= 0.0_r) return 0.0_r;

    Real H = 0.0_r;
    Real const inv_total = 1.0_r / total;
    #pragma omp parallel for reduction(+:H) schedule(static)
    for (Index k = 0; k < n; ++k) {
        Real const px = p[k] * inv_total;
        if (px > 1e-30_r) {
            H -= px * std::log(px);
        }
    }
    return H;
}

// Normalized "focus" score: 1 when all energy is in one cell (delta),
// 0 when uniformly distributed (entropy = log(M)).
inline Real focus_score(Field<Real, 2> const& e) noexcept {
    Real const H = energy_entropy(e);
    Real const log_M = std::log(static_cast<Real>(e.size()));
    if (log_M <= 0.0_r) return 0.0_r;
    return 1.0_r - H / log_M;
}

// ----------------------------------------------------------------------
// Hotspot concentration (§24): fraction of energy inside Ω_hot ⊂ Ω.
// `mask` is a binary 0/1 field over the same grid as `e`; cells where
// mask > 0 belong to Ω_hot.
// ----------------------------------------------------------------------
inline Real hotspot_concentration(Field<Real, 2> const& e,
                                  Field<Real, 2> const& mask) noexcept {
    Index const n = e.size();
    Real const* pe = e.data();
    Real const* pm = mask.data();

    Real sum_in  = 0.0_r;
    Real sum_all = 0.0_r;
    #pragma omp parallel for reduction(+:sum_in,sum_all) schedule(static)
    for (Index k = 0; k < n; ++k) {
        sum_all += pe[k];
        if (pm[k] > 0.0_r) sum_in += pe[k];
    }
    if (sum_all <= 0.0_r) return 0.0_r;
    return sum_in / sum_all;
}

// Rectangular region overload (Phase 2-style probe boxes).
inline Real hotspot_concentration_region(Field<Real, 2> const& e,
                                         IVec<2> lo, IVec<2> hi) noexcept {
    auto const& g = e.grid();
    Index const nx = g.shape[0];
    Index const ny = g.shape[1];
    Real sum_in  = 0.0_r;
    Real sum_all = 0.0_r;
    #pragma omp parallel for reduction(+:sum_in,sum_all) collapse(2) schedule(static)
    for (Index i = 0; i < nx; ++i) {
        for (Index j = 0; j < ny; ++j) {
            Real const v = e(i, j);
            sum_all += v;
            if (i >= lo[0] && i <= hi[0] && j >= lo[1] && j <= hi[1]) sum_in += v;
        }
    }
    if (sum_all <= 0.0_r) return 0.0_r;
    return sum_in / sum_all;
}

// ----------------------------------------------------------------------
// Ligand occlusion (§25): fraction of EMPTY-pocket energy that resides
// in the region the ligand would occupy. High value = ligand is
// blocking important wave pathways.
// ----------------------------------------------------------------------
inline Real ligand_occlusion(Field<Real, 2> const& e_empty,
                             Field<Real, 2> const& ligand_mask) noexcept {
    return hotspot_concentration(e_empty, ligand_mask);
}

} // namespace wavelab
