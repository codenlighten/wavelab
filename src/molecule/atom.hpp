#pragma once
//
// Atom<D>: a single point-particle in D-dimensional space with the
// parameters needed to splat itself onto the wavelab field grids.
//
//   pos      — world-coord position (e.g. Ångströms)
//   sigma    — Gaussian radius for density splat (§6); roughly van-der-
//              Waals-derived
//   charge   — partial charge for Q(x) splat (§10)
//   hydro    — hydrophobicity contribution for H(x) splat (§11)
//   element  — atomic number, kept for downstream parameter tables
//
// Phase 3 wires a real element → (sigma, charge, hydro) parameter table.
// Phase 2 uses caller-supplied values on synthetic scenes.
//

#include "core/types.hpp"

namespace wavelab {

template <int D>
struct Atom {
    Vec<D> pos    { };
    Real   sigma  = 1.0_r;
    Real   charge = 0.0_r;
    Real   hydro  = 0.0_r;
    int    element = 0;   // atomic number, 0 = unspecified
};

} // namespace wavelab
