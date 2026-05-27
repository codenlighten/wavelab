#pragma once
//
// Minimal PDB parser — reads ATOM (and optionally HETATM) records from a
// .pdb file and produces a MolecularScene<3>.
//
// We deliberately do not link gemmi/RDKit/OpenBabel: our consumer only
// needs (position, element, σ, q, h) per atom, and PDB ATOM records are
// fixed-column text that's trivial to parse correctly.
//
// Filtering supported:
//   * alt-loc — keep records where alt_loc is ' ' or matches the chosen
//     preference (default 'A'). Without this, side-chains with both A/B
//     conformers silently double-count, distorting density.
//   * waters / ions — excluded from HETATM by a hard-coded residue
//     blocklist (HOH, WAT, NA, CL, ...). Toggle via include_waters_ions.
//   * residue allowlist — when set, only atoms whose residue name is in
//     the set are kept. Used to extract just a specific ligand (e.g.,
//     SAQ from 1HXB, MK1 from 1HSG).
//

#include "core/types.hpp"
#include "molecule/scene.hpp"

#include <filesystem>
#include <iosfwd>
#include <optional>
#include <set>
#include <string>
#include <string_view>

namespace wavelab {

struct PdbParseOptions {
    bool include_hetatm           = false;   // typically waters / ligands
    bool include_waters_ions      = false;   // include common solvent/buffer HETATM (HOH, GOL, ...)
    bool parse_bfactor_as_charge  = false;   // overview §10: pack charge in B-factor column
    char alt_loc_preference       = 'A';     // when multiple alt-locs, pick this one
    Real pad_box                  = 5.0_r;   // padding around min/max coords for the scene box (Å)
    // When set, only atoms whose 3-letter residue name is in the set
    // are kept. Useful for extracting a single ligand from a co-crystal.
    std::optional<std::set<std::string>> keep_residue_allowlist;
};

struct PdbParseResult {
    MolecularScene<3> scene;
    Index atom_count       = 0;    // accepted ATOM records
    Index hetatm_count     = 0;    // accepted HETATM records
    Index skipped_alt_loc  = 0;
    Index skipped_waters_ions = 0;
    Index skipped_residue  = 0;    // rejected by allowlist
};

// Parse a PDB stream / file. Throws std::runtime_error on I/O failure.
// Malformed individual ATOM lines are skipped silently.
PdbParseResult parse_pdb(std::istream& in, PdbParseOptions const& opts = {});
PdbParseResult parse_pdb_file(std::filesystem::path const& path,
                              PdbParseOptions const& opts = {});

// Convenience: parse a PDB and keep only HETATM records with the given
// residue name (e.g., "SAQ", "MK1", "BEN", "BTN"). Returns a
// MolecularScene<3> containing just that ligand's atoms.
MolecularScene<3> parse_pdb_hetatm_residue(std::filesystem::path const& path,
                                           std::string const& residue);

} // namespace wavelab
