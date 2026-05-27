#pragma once
//
// Synthetic geometry generators — programmatic MolecularScenes for
// validation and demos. All return MolecularScene<2> for Phase 2; the
// 3D variants land in Phase 5.
//
// Patterns:
//   * single_atom          — one Gaussian, sanity check for splatting
//   * dumbbell             — two atoms separated; tests interaction at
//                            a known distance
//   * pocket_arc           — ring of atoms with an opening on one side,
//                            mimics a binding pocket cross-section
//   * slab_atoms           — uniform sheet of atoms forming a wall
//                            (for Fresnel-style reflection tests when
//                            paired with a refractive-index medium)
//
// For pure refractive-index slabs (no atoms), use medium/build_medium
// directly — see apply_index_slab.
//

#include "molecule/atom.hpp"
#include "molecule/scene.hpp"

#include <cmath>
#include <numbers>

namespace wavelab {

inline MolecularScene<2> single_atom_2d(Vec<2> pos, Real sigma,
                                        Vec<2> box_min, Vec<2> box_max) {
    MolecularScene<2> s;
    s.atoms.push_back(Atom<2>{pos, sigma, 0.0_r, 0.0_r, 0});
    s.box_min = box_min;
    s.box_max = box_max;
    return s;
}

inline MolecularScene<2> dumbbell_2d(Vec<2> a, Vec<2> b, Real sigma,
                                     Vec<2> box_min, Vec<2> box_max) {
    MolecularScene<2> s;
    s.atoms.push_back(Atom<2>{a, sigma, 0.0_r, 0.0_r, 0});
    s.atoms.push_back(Atom<2>{b, sigma, 0.0_r, 0.0_r, 0});
    s.box_min = box_min;
    s.box_max = box_max;
    return s;
}

// Ring of `count` atoms at radius R around `center`, with a wedge of
// `opening_rad` radians missing on the side facing `opening_dir_rad`.
// Approximates a 2D cross-section through a binding pocket.
inline MolecularScene<2> pocket_arc_2d(Vec<2> center, Real R, Real sigma,
                                       Index count,
                                       Real opening_dir_rad,
                                       Real opening_rad,
                                       Vec<2> box_min, Vec<2> box_max) {
    MolecularScene<2> s;
    s.box_min = box_min;
    s.box_max = box_max;
    Real const half_open = 0.5_r * opening_rad;
    for (Index k = 0; k < count; ++k) {
        Real const theta = static_cast<Real>(2.0 * std::numbers::pi)
                         * static_cast<Real>(k) / static_cast<Real>(count);
        // Skip atoms whose angle is inside the opening wedge.
        Real d = std::fmod(theta - opening_dir_rad
                           + static_cast<Real>(2.0 * std::numbers::pi),
                           static_cast<Real>(2.0 * std::numbers::pi));
        if (d > static_cast<Real>(std::numbers::pi)) {
            d = static_cast<Real>(2.0 * std::numbers::pi) - d;
        }
        if (d < half_open) continue;
        Vec<2> p{
            center[0] + R * std::cos(theta),
            center[1] + R * std::sin(theta)
        };
        s.atoms.push_back(Atom<2>{p, sigma, 0.0_r, 0.0_r, 0});
    }
    return s;
}

// Densely-packed atoms in a slab x ∈ [x_lo, x_hi], y spanning the box.
// Use for atom-based scattering tests; for clean Fresnel R/T checks
// prefer apply_index_slab on the Medium directly (no Gaussian smearing).
inline MolecularScene<2> slab_atoms_2d(Real x_lo, Real x_hi,
                                       Real spacing, Real sigma,
                                       Vec<2> box_min, Vec<2> box_max) {
    MolecularScene<2> s;
    s.box_min = box_min;
    s.box_max = box_max;
    for (Real x = x_lo; x <= x_hi; x += spacing) {
        for (Real y = box_min[1]; y <= box_max[1]; y += spacing) {
            s.atoms.push_back(Atom<2>{Vec<2>{x, y}, sigma, 0.0_r, 0.0_r, 0});
        }
    }
    return s;
}

} // namespace wavelab
