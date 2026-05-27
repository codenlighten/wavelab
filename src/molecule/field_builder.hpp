#pragma once
//
// Atoms → scalar Fields by Gaussian splatting.
//
//   ρ(x) = Σ_a       exp(-|x - r_a|² / (2 σ_a²))         (§6 density)
//   Q(x) = Σ_a q_a · exp(-|x - r_a|² / (2 σ_q²))         (§10 charge)
//   H(x) = Σ_a h_a · exp(-|x - r_a|² / (2 σ_h²))         (§11 hydro)
//
// Each splat is windowed to a cube/square of ±4σ around the atom so the
// cost is O(N_atoms · σ³) rather than O(N_atoms · N_cells).
//
// Inputs use the field's own grid; the scene's box_min/box_max is not
// required to align with the grid.
//

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"
#include "molecule/scene.hpp"

#include <algorithm>
#include <cmath>

namespace wavelab {

namespace detail {

// Cell index containing world coordinate `w` on axis `i` of grid `g`.
template <int D>
inline Index world_to_cell(Grid<D> const& g, Real w, int axis) noexcept {
    return static_cast<Index>(std::floor((w - g.origin[static_cast<std::size_t>(axis)])
                                          / g.spacing[static_cast<std::size_t>(axis)]));
}

template <int D>
inline Index clamp_index(Index v, Index lo, Index hi) noexcept {
    return std::max(lo, std::min(v, hi));
}

} // namespace detail

// Add Gaussian splat of `atoms` into the existing `out` field.
// `out` is added to, not overwritten — caller may pre-zero or accumulate.
template <int D>
inline void splat_density(Field<Real, D>&            out,
                          MolecularScene<D> const&   scene,
                          Real                       cutoff_sigmas = 4.0_r) {
    auto const& g  = out.grid();
    Real const  dx = g.spacing[0];

    for (Atom<D> const& a : scene.atoms) {
        Real const sig2_inv_half = 0.5_r / (a.sigma * a.sigma);
        Real const radius        = cutoff_sigmas * a.sigma;
        Index const r_cells      = static_cast<Index>(std::ceil(radius / dx));

        if constexpr (D == 1) {
            Index const c0 = detail::world_to_cell<D>(g, a.pos[0], 0);
            Index const i_lo = detail::clamp_index<D>(c0 - r_cells, 0, g.shape[0] - 1);
            Index const i_hi = detail::clamp_index<D>(c0 + r_cells, 0, g.shape[0] - 1);
            for (Index i = i_lo; i <= i_hi; ++i) {
                Real const xc = g.origin[0] + (static_cast<Real>(i) + 0.5_r) * g.spacing[0];
                Real const dxv = xc - a.pos[0];
                out(i) += std::exp(-(dxv * dxv) * sig2_inv_half);
            }
        } else if constexpr (D == 2) {
            Index const c0 = detail::world_to_cell<D>(g, a.pos[0], 0);
            Index const c1 = detail::world_to_cell<D>(g, a.pos[1], 1);
            Index const i_lo = detail::clamp_index<D>(c0 - r_cells, 0, g.shape[0] - 1);
            Index const i_hi = detail::clamp_index<D>(c0 + r_cells, 0, g.shape[0] - 1);
            Index const j_lo = detail::clamp_index<D>(c1 - r_cells, 0, g.shape[1] - 1);
            Index const j_hi = detail::clamp_index<D>(c1 + r_cells, 0, g.shape[1] - 1);
            for (Index i = i_lo; i <= i_hi; ++i) {
                Real const xc = g.origin[0] + (static_cast<Real>(i) + 0.5_r) * g.spacing[0];
                Real const dxv = xc - a.pos[0];
                for (Index j = j_lo; j <= j_hi; ++j) {
                    Real const yc = g.origin[1] + (static_cast<Real>(j) + 0.5_r) * g.spacing[1];
                    Real const dyv = yc - a.pos[1];
                    out(i, j) += std::exp(-(dxv * dxv + dyv * dyv) * sig2_inv_half);
                }
            }
        } else if constexpr (D == 3) {
            Index const c0 = detail::world_to_cell<D>(g, a.pos[0], 0);
            Index const c1 = detail::world_to_cell<D>(g, a.pos[1], 1);
            Index const c2 = detail::world_to_cell<D>(g, a.pos[2], 2);
            Index const i_lo = detail::clamp_index<D>(c0 - r_cells, 0, g.shape[0] - 1);
            Index const i_hi = detail::clamp_index<D>(c0 + r_cells, 0, g.shape[0] - 1);
            Index const j_lo = detail::clamp_index<D>(c1 - r_cells, 0, g.shape[1] - 1);
            Index const j_hi = detail::clamp_index<D>(c1 + r_cells, 0, g.shape[1] - 1);
            Index const k_lo = detail::clamp_index<D>(c2 - r_cells, 0, g.shape[2] - 1);
            Index const k_hi = detail::clamp_index<D>(c2 + r_cells, 0, g.shape[2] - 1);
            for (Index i = i_lo; i <= i_hi; ++i) {
                Real const xc = g.origin[0] + (static_cast<Real>(i) + 0.5_r) * g.spacing[0];
                Real const dxv = xc - a.pos[0];
                for (Index j = j_lo; j <= j_hi; ++j) {
                    Real const yc = g.origin[1] + (static_cast<Real>(j) + 0.5_r) * g.spacing[1];
                    Real const dyv = yc - a.pos[1];
                    for (Index k = k_lo; k <= k_hi; ++k) {
                        Real const zc = g.origin[2] + (static_cast<Real>(k) + 0.5_r) * g.spacing[2];
                        Real const dzv = zc - a.pos[2];
                        out(i, j, k) += std::exp(
                            -(dxv * dxv + dyv * dyv + dzv * dzv) * sig2_inv_half);
                    }
                }
            }
        }
    }
}

// Same shape as splat_density but weighted by atom.charge / atom.hydro.
template <int D>
inline void splat_charge(Field<Real, D>&          out,
                         MolecularScene<D> const& scene,
                         Real cutoff_sigmas = 4.0_r) {
    splat_weighted_<D>(out, scene, cutoff_sigmas,
                       [](Atom<D> const& a) { return a.charge; });
}

template <int D>
inline void splat_hydro(Field<Real, D>&          out,
                        MolecularScene<D> const& scene,
                        Real cutoff_sigmas = 4.0_r) {
    splat_weighted_<D>(out, scene, cutoff_sigmas,
                       [](Atom<D> const& a) { return a.hydro; });
}

// --- internal helper, weighted version --------------------------------

template <int D, typename F>
inline void splat_weighted_(Field<Real, D>&          out,
                            MolecularScene<D> const& scene,
                            Real                     cutoff_sigmas,
                            F                        weight) {
    auto const& g  = out.grid();
    Real const  dx = g.spacing[0];

    for (Atom<D> const& a : scene.atoms) {
        Real const w = weight(a);
        if (w == Real{0}) continue;
        Real const sig2_inv_half = 0.5_r / (a.sigma * a.sigma);
        Real const radius        = cutoff_sigmas * a.sigma;
        Index const r_cells      = static_cast<Index>(std::ceil(radius / dx));

        if constexpr (D == 2) {
            Index const c0 = detail::world_to_cell<D>(g, a.pos[0], 0);
            Index const c1 = detail::world_to_cell<D>(g, a.pos[1], 1);
            Index const i_lo = detail::clamp_index<D>(c0 - r_cells, 0, g.shape[0] - 1);
            Index const i_hi = detail::clamp_index<D>(c0 + r_cells, 0, g.shape[0] - 1);
            Index const j_lo = detail::clamp_index<D>(c1 - r_cells, 0, g.shape[1] - 1);
            Index const j_hi = detail::clamp_index<D>(c1 + r_cells, 0, g.shape[1] - 1);
            for (Index i = i_lo; i <= i_hi; ++i) {
                Real const xc = g.origin[0] + (static_cast<Real>(i) + 0.5_r) * g.spacing[0];
                Real const dxv = xc - a.pos[0];
                for (Index j = j_lo; j <= j_hi; ++j) {
                    Real const yc = g.origin[1] + (static_cast<Real>(j) + 0.5_r) * g.spacing[1];
                    Real const dyv = yc - a.pos[1];
                    out(i, j) += w * std::exp(-(dxv * dxv + dyv * dyv) * sig2_inv_half);
                }
            }
        } else {
            static_assert(D == 2, "splat_weighted_: D=1/3 add as needed");
        }
    }
}

} // namespace wavelab
