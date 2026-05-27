#include "molecule/parser_pdb.hpp"

#include "molecule/element_params.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace wavelab {

namespace {

// Strip leading/trailing whitespace from a substring view.
std::string_view trim(std::string_view s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))  s.remove_suffix(1);
    return s;
}

// PDB columns we read (1-indexed in the spec, converted to 0-indexed substring).
constexpr std::size_t kColRecord      = 0,  kLenRecord      = 6;
constexpr std::size_t kColAtomName    = 12, kLenAtomName    = 4;
constexpr std::size_t kColAltLoc      = 16, kLenAltLoc      = 1;
constexpr std::size_t kColResName     = 17, kLenResName     = 3;
constexpr std::size_t kColX           = 30, kLenX           = 8;
constexpr std::size_t kColY           = 38, kLenY           = 8;
constexpr std::size_t kColZ           = 46, kLenZ           = 8;
constexpr std::size_t kColBfactor     = 60, kLenBfactor     = 6;
constexpr std::size_t kColElement     = 76, kLenElement     = 2;

// Returns false if the field is empty / out of range / unparseable.
bool parse_real_field(std::string_view line, std::size_t col, std::size_t len, Real& out) {
    if (line.size() < col + len) return false;
    auto sv = trim(line.substr(col, len));
    if (sv.empty()) return false;
    // std::from_chars on float — works for non-special floating-point text.
    // PDB uses fixed-format like "  37.123" or " -1.456".
    float v;
    auto const* first = sv.data();
    auto const* last  = first + sv.size();
    auto res = std::from_chars(first, last, v);
    if (res.ec != std::errc{} || res.ptr != last) {
        // Fallback: strtof for nonstandard layouts.
        try {
            v = std::stof(std::string{sv});
        } catch (...) { return false; }
    }
    out = static_cast<Real>(v);
    return true;
}

// Pull the element symbol — preferring columns 77-78 (PDB v3.3), falling
// back to the first two non-digit characters of the atom name. PDB ATOM
// names follow conventions like " CA ", " N  ", " O1 " where the element
// is the left-justified non-digit portion.
std::string element_symbol_of(std::string_view line) {
    if (line.size() >= kColElement + kLenElement) {
        auto e = trim(line.substr(kColElement, kLenElement));
        if (!e.empty()) return std::string{e};
    }
    if (line.size() >= kColAtomName + kLenAtomName) {
        auto name = line.substr(kColAtomName, kLenAtomName);
        std::string sym;
        for (char c : name) {
            if (std::isspace(static_cast<unsigned char>(c))) {
                if (sym.empty()) continue;
                break;
            }
            if (std::isdigit(static_cast<unsigned char>(c))) break;
            sym.push_back(c);
            if (sym.size() == 2) break;
        }
        return sym;
    }
    return {};
}

// Standard solvent / buffer / counter-ion residue names commonly seen
// in HETATM blocks of crystallographic PDBs. Default behavior excludes
// these so a co-crystal load doesn't silently splat waters into the
// medium.
bool is_water_or_ion(std::string_view res) {
    static const std::set<std::string> kBlocklist = {
        "HOH", "WAT", "DOD", "H2O",                          // water
        "NA", "K", "CL", "BR", "I", "F",                     // monovalent ions
        "MG", "CA", "ZN", "FE", "MN", "CU", "NI", "CO",      // metals
        "SO4", "PO4", "HPO", "NO3",                          // counter-ions
        "GOL", "EDO", "PEG", "DMS", "MPD", "FMT", "ACT",     // cryoprotectants / buffers
    };
    return kBlocklist.contains(std::string{res});
}

} // namespace

PdbParseResult parse_pdb(std::istream& in, PdbParseOptions const& opts) {
    PdbParseResult result;
    auto& scene = result.scene;

    Real xmin =  std::numeric_limits<Real>::max();
    Real ymin =  std::numeric_limits<Real>::max();
    Real zmin =  std::numeric_limits<Real>::max();
    Real xmax = -std::numeric_limits<Real>::max();
    Real ymax = -std::numeric_limits<Real>::max();
    Real zmax = -std::numeric_limits<Real>::max();

    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < kColZ + kLenZ) continue;  // too short to contain coords
        std::string_view sv{line};
        auto record = trim(sv.substr(kColRecord, kLenRecord));

        bool const is_atom   = (record == "ATOM");
        bool const is_hetatm = (record == "HETATM");
        if (!is_atom && (!opts.include_hetatm || !is_hetatm)) continue;

        // alt-loc filter
        char alt = sv.size() > kColAltLoc ? sv[kColAltLoc] : ' ';
        if (alt != ' ' && alt != opts.alt_loc_preference) {
            ++result.skipped_alt_loc;
            continue;
        }

        // Residue-name filtering (waters/ions + optional allowlist).
        std::string_view res_name;
        if (sv.size() >= kColResName + kLenResName) {
            res_name = trim(sv.substr(kColResName, kLenResName));
        }
        if (is_hetatm && !opts.include_waters_ions && is_water_or_ion(res_name)) {
            ++result.skipped_waters_ions;
            continue;
        }
        if (opts.keep_residue_allowlist) {
            if (!opts.keep_residue_allowlist->contains(std::string{res_name})) {
                ++result.skipped_residue;
                continue;
            }
        }

        Real x, y, z;
        if (!parse_real_field(sv, kColX, kLenX, x)) continue;
        if (!parse_real_field(sv, kColY, kLenY, y)) continue;
        if (!parse_real_field(sv, kColZ, kLenZ, z)) continue;

        Atom<3> a{};
        a.pos = Vec<3>{x, y, z};
        auto sym = element_symbol_of(sv);
        assign_atom_params_by_symbol<3>(a, sym);

        if (opts.parse_bfactor_as_charge) {
            Real bf;
            if (parse_real_field(sv, kColBfactor, kLenBfactor, bf)) {
                a.charge = bf;
            }
        }

        scene.atoms.push_back(a);
        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
        if (z < zmin) zmin = z;
        if (z > zmax) zmax = z;

        if (is_atom)   ++result.atom_count;
        if (is_hetatm) ++result.hetatm_count;
    }

    if (!scene.atoms.empty()) {
        Real const p = opts.pad_box;
        scene.box_min = Vec<3>{xmin - p, ymin - p, zmin - p};
        scene.box_max = Vec<3>{xmax + p, ymax + p, zmax + p};
    }
    return result;
}

PdbParseResult parse_pdb_file(std::filesystem::path const& path,
                              PdbParseOptions const& opts) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("parse_pdb_file: cannot open " + path.string());
    }
    return parse_pdb(in, opts);
}

MolecularScene<3> parse_pdb_hetatm_residue(std::filesystem::path const& path,
                                           std::string const& residue) {
    PdbParseOptions opts;
    opts.include_hetatm = true;
    opts.include_waters_ions = true;        // allowlist is the gate; don't double-filter
    opts.keep_residue_allowlist = std::set<std::string>{residue};
    auto r = parse_pdb_file(path, opts);
    return std::move(r.scene);
}

} // namespace wavelab
