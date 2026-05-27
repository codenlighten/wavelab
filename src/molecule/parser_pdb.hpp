#pragma once
//
// Minimal PDB parser — reads ATOM (and optionally HETATM) records from a
// .pdb file and produces a MolecularScene<3>.
//
// We deliberately do not link gemmi/RDKit/OpenBabel: our consumer only
// needs (position, element, σ, q, h) per atom, and PDB ATOM records are
// fixed-column text that's trivial to parse correctly. Anything beyond
// canonical PDB v3.3 ATOM (multi-model NMR ensembles, mmCIF, ligand
// chemistry, bond orders, alt-loc tiebreaks beyond first) is out of scope
// for Phase 3 and can be added when needed.
//

#include "core/types.hpp"
#include "molecule/scene.hpp"

#include <filesystem>
#include <iosfwd>
#include <string_view>

namespace wavelab {

struct PdbParseOptions {
    bool include_hetatm           = false;   // typically waters / ligands
    bool parse_bfactor_as_charge  = false;   // overview §10: pack charge in B-factor column
    char alt_loc_preference       = 'A';     // when multiple alt-locs, pick this one ( ' ' or 'A' as default)
    Real pad_box                  = 5.0_r;   // padding around min/max coords for the scene box (Å)
};

struct PdbParseResult {
    MolecularScene<3> scene;
    Index atom_count     = 0;
    Index hetatm_count   = 0;
    Index skipped_alt_loc = 0;
};

// Parse a PDB stream / file. Throws std::runtime_error on I/O failure.
// Malformed individual ATOM lines are skipped silently (counted toward
// `skipped_*` in the result if you add tracking).
PdbParseResult parse_pdb(std::istream& in, PdbParseOptions const& opts = {});
PdbParseResult parse_pdb_file(std::filesystem::path const& path,
                              PdbParseOptions const& opts = {});

} // namespace wavelab
