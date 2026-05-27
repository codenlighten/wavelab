#pragma once
//
// Element → (sigma_vdW, default_partial_charge, default_hydrophobicity)
// lookup table.
//
// The σ values are van-der-Waals radii in Ångströms — they control the
// width of the Gaussian splat in density / charge / hydro field builders.
// We use Bondi (1964) consensus values for the most common elements.
//
// Charge and hydrophobicity defaults are chemistry-naive placeholders
// (zero charge, hydro chosen from a simple polar/non-polar split). Real
// applications should override these per atom (AMBER/CHARMM force field
// assignments, charge calculations from QM, etc.). They live here so a
// fresh PDB load lands with usable-but-rough values out of the box.
//

#include "core/types.hpp"
#include "molecule/atom.hpp"

#include <cstring>
#include <string_view>

namespace wavelab {

struct ElementParams {
    int             z;          // atomic number
    char            symbol[3];  // 2-char + NUL
    Real            sigma;      // van-der-Waals radius (Å)
    Real            charge;     // default partial charge
    Real            hydro;      // hydrophobicity contribution
};

// Bondi (1964) van-der-Waals radii for elements we care about most.
//
// Partial charges are *electronegativity-derived pseudo-charges*, NOT
// proper atom-type-aware partials (no Gasteiger / AMBER assignment yet).
// They give the right qualitative ordering — O more negative than N
// more negative than C; H slightly positive on average since most H
// in organic molecules is bonded to C-N-O scaffolds and inherits some
// dipolar character. Real partial charges per molecule are a Phase 8
// concern (RDKit / Gasteiger / QM input).
//
// Hydrophobicity convention: positive = hydrophobic (C, S, halogens),
// negative = hydrophilic (N, O), zero = neutral (H, ions).
inline constexpr ElementParams kElementTable[] = {
    //  Z  sym   σ(Å)    q       h
    {   1, "H ", 1.20_r,  0.1_r,  0.0_r},
    {   6, "C ", 1.70_r,  0.0_r,  1.0_r},
    {   7, "N ", 1.55_r, -0.4_r, -1.0_r},
    {   8, "O ", 1.52_r, -0.5_r, -1.0_r},
    {  15, "P ", 1.80_r,  0.5_r,  0.0_r},
    {  16, "S ", 1.80_r, -0.1_r,  0.5_r},
    {   9, "F ", 1.47_r, -0.6_r, -0.5_r},
    {  17, "CL", 1.75_r, -0.3_r,  0.5_r},
    {  35, "BR", 1.85_r, -0.2_r,  0.5_r},
    {  53, "I ", 1.98_r, -0.1_r,  0.5_r},
    {  11, "NA", 2.27_r,  1.0_r,  0.0_r},
    {  19, "K ", 2.75_r,  1.0_r,  0.0_r},
    {  12, "MG", 1.73_r,  2.0_r,  0.0_r},
    {  20, "CA", 2.31_r,  2.0_r,  0.0_r},
    {  26, "FE", 2.00_r,  2.0_r,  0.0_r},
    {  30, "ZN", 1.39_r,  2.0_r,  0.0_r},
};

// Fallback for unrecognized elements.
inline constexpr ElementParams kElementUnknown{
    0, "X ", 1.70_r, 0.0_r, 0.0_r
};

namespace detail {

// Uppercase a 2-char element string in place (PDB columns can be mixed-case
// for two-letter elements like "Na", "Cl"). Pads to length 2 with space.
inline void normalize_element_symbol(char out[3], std::string_view in) noexcept {
    out[0] = (in.size() > 0) ? static_cast<char>(std::toupper(static_cast<unsigned char>(in[0]))) : ' ';
    out[1] = (in.size() > 1) ? static_cast<char>(std::toupper(static_cast<unsigned char>(in[1]))) : ' ';
    out[2] = '\0';
}

} // namespace detail

inline ElementParams const& element_params_by_symbol(std::string_view sym) noexcept {
    char norm[3];
    detail::normalize_element_symbol(norm, sym);
    for (auto const& e : kElementTable) {
        if (e.symbol[0] == norm[0] && e.symbol[1] == norm[1]) return e;
    }
    return kElementUnknown;
}

inline ElementParams const& element_params_by_z(int z) noexcept {
    for (auto const& e : kElementTable) {
        if (e.z == z) return e;
    }
    return kElementUnknown;
}

// Populate atom.sigma/charge/hydro from element_params (overwrites prior
// values; caller can override after).
template <int D>
inline void assign_atom_params_by_symbol(Atom<D>& a, std::string_view sym) noexcept {
    auto const& p = element_params_by_symbol(sym);
    a.element = p.z;
    a.sigma   = p.sigma;
    a.charge  = p.charge;
    a.hydro   = p.hydro;
}

} // namespace wavelab
