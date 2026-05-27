#pragma once
//
// Project a MolecularScene<3> onto a 2D slice along the z-axis.
//
// Atoms with z ∈ [z_lo, z_hi] are kept; their (x, y) is the 2D position.
// The 2D box is the 3D box's (x, y) extent (z is collapsed).
//
// This is the bridge that lets the Phase 2 2D engine consume real PDB
// data — until Phase 5 brings up the 3D stepper.
//
// If you'd rather slice along x or y, just permute coordinates in the
// scene before calling — keeps this function shape-simple.
//

#include "core/types.hpp"
#include "molecule/scene.hpp"

namespace wavelab {

inline MolecularScene<2> slice_scene_xy(MolecularScene<3> const& src,
                                        Real z_lo, Real z_hi) {
    MolecularScene<2> out;
    out.box_min = Vec<2>{src.box_min[0], src.box_min[1]};
    out.box_max = Vec<2>{src.box_max[0], src.box_max[1]};
    for (Atom<3> const& a3 : src.atoms) {
        if (a3.pos[2] < z_lo || a3.pos[2] > z_hi) continue;
        Atom<2> a2{};
        a2.pos    = Vec<2>{a3.pos[0], a3.pos[1]};
        a2.sigma  = a3.sigma;
        a2.charge = a3.charge;
        a2.hydro  = a3.hydro;
        a2.element = a3.element;
        out.atoms.push_back(a2);
    }
    return out;
}

// Centered slice of given thickness around z_mid.
inline MolecularScene<2> slice_scene_xy_centered(MolecularScene<3> const& src,
                                                 Real z_mid, Real thickness) {
    Real const half = 0.5_r * thickness;
    return slice_scene_xy(src, z_mid - half, z_mid + half);
}

// Convenience: slice through the geometric midplane of the scene's box.
inline MolecularScene<2> slice_scene_xy_midplane(MolecularScene<3> const& src,
                                                 Real thickness) {
    Real const z_mid = 0.5_r * (src.box_min[2] + src.box_max[2]);
    return slice_scene_xy_centered(src, z_mid, thickness);
}

} // namespace wavelab
