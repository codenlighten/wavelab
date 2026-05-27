#pragma once
//
// Minimal SDF V2000 parser — reads the FIRST molecule from an SDF or
// MOL file and produces a MolecularScene<3>.
//
// Format reference: CTfile (V2000) atom block, columns:
//   * x, y, z   — 10-char fixed-width right-aligned floats
//   * element   — left-justified at column 31 (1-indexed), width 3
//   * additional fields ignored
//
// Bonds, stereo, charge-block extensions, property blocks (M  CHG, etc.)
// are all skipped — the wave engine only consumes atom positions. If
// real charge needs come up later, parsing the M  CHG block lives here.
//

#include "core/types.hpp"
#include "molecule/scene.hpp"

#include <filesystem>
#include <iosfwd>

namespace wavelab {

struct SdfParseOptions {
    Real pad_box = 5.0_r;   // padding around min/max coords (Å)
};

struct SdfParseResult {
    MolecularScene<3> scene;
    Index             atom_count = 0;
    Index             bond_count = 0;
    std::string       molecule_name;
};

SdfParseResult parse_sdf(std::istream& in, SdfParseOptions const& opts = {});
SdfParseResult parse_sdf_file(std::filesystem::path const& path,
                              SdfParseOptions const& opts = {});

} // namespace wavelab
