#include "molecule/parser_sdf.hpp"

#include "molecule/element_params.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace wavelab {

namespace {

std::string_view trim(std::string_view s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))  s.remove_suffix(1);
    return s;
}

bool parse_real(std::string_view sv, Real& out) {
    sv = trim(sv);
    if (sv.empty()) return false;
    try {
        out = static_cast<Real>(std::stof(std::string{sv}));
        return true;
    } catch (...) { return false; }
}

bool parse_int(std::string_view sv, int& out) {
    sv = trim(sv);
    if (sv.empty()) return false;
    try {
        out = std::stoi(std::string{sv});
        return true;
    } catch (...) { return false; }
}

} // namespace

SdfParseResult parse_sdf(std::istream& in, SdfParseOptions const& opts) {
    SdfParseResult result;
    std::string line;

    // --- Header block (3 lines) ---
    if (!std::getline(in, line)) throw std::runtime_error("parse_sdf: empty stream");
    result.molecule_name = std::string{trim(line)};
    if (!std::getline(in, line)) throw std::runtime_error("parse_sdf: missing line 2");
    if (!std::getline(in, line)) throw std::runtime_error("parse_sdf: missing line 3");

    // --- Counts line ---
    if (!std::getline(in, line)) throw std::runtime_error("parse_sdf: missing counts line");
    if (line.size() < 6) throw std::runtime_error("parse_sdf: counts line too short");
    int n_atoms = 0, n_bonds = 0;
    if (!parse_int(std::string_view{line}.substr(0, 3), n_atoms)) {
        throw std::runtime_error("parse_sdf: bad atom count");
    }
    if (line.size() >= 6) parse_int(std::string_view{line}.substr(3, 3), n_bonds);
    result.atom_count = n_atoms;
    result.bond_count = n_bonds;

    // --- Atom block ---
    auto& scene = result.scene;
    scene.atoms.reserve(static_cast<std::size_t>(n_atoms));

    Real xmin =  std::numeric_limits<Real>::max();
    Real ymin =  std::numeric_limits<Real>::max();
    Real zmin =  std::numeric_limits<Real>::max();
    Real xmax = -std::numeric_limits<Real>::max();
    Real ymax = -std::numeric_limits<Real>::max();
    Real zmax = -std::numeric_limits<Real>::max();

    for (int k = 0; k < n_atoms; ++k) {
        if (!std::getline(in, line)) {
            throw std::runtime_error("parse_sdf: atom block truncated");
        }
        if (line.size() < 34) continue;   // too short to contain element

        Real x, y, z;
        std::string_view sv{line};
        if (!parse_real(sv.substr(0,  10), x)) continue;
        if (!parse_real(sv.substr(10, 10), y)) continue;
        if (!parse_real(sv.substr(20, 10), z)) continue;

        std::string sym{trim(sv.substr(31, 3))};

        Atom<3> a{};
        a.pos = Vec<3>{x, y, z};
        assign_atom_params_by_symbol<3>(a, sym);
        scene.atoms.push_back(a);

        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
        if (z < zmin) zmin = z;
        if (z > zmax) zmax = z;
    }

    // We deliberately skip the bond block, property block, and any
    // additional molecules in the SDF — caller can re-invoke after
    // seeking past the next "$$$$" if they need a multi-record reader.

    if (!scene.atoms.empty()) {
        Real const p = opts.pad_box;
        scene.box_min = Vec<3>{xmin - p, ymin - p, zmin - p};
        scene.box_max = Vec<3>{xmax + p, ymax + p, zmax + p};
    }
    return result;
}

SdfParseResult parse_sdf_file(std::filesystem::path const& path,
                              SdfParseOptions const& opts) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("parse_sdf_file: cannot open " + path.string());
    }
    return parse_sdf(in, opts);
}

} // namespace wavelab
