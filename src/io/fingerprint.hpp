#pragma once
//
// Fingerprint serialization — a versioned JSON blob describing one
// simulated-light measurement (scene + scores + spectral vector + run
// metadata). Hand-rolled JSON writer/reader to avoid pulling in another
// dependency for what's a flat record.
//
//   {
//     "wavelab_fingerprint_version": 1,
//     "scene":   { ... source/grid/medium config ... },
//     "scalars": { name: value, ... },
//     "spectral": [v0, v1, ..., vN-1],
//     "spectral_freqs": [f0, f1, ..., fN-1]
//   }
//

#include "core/types.hpp"

#include <filesystem>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

namespace wavelab {

struct Fingerprint {
    static constexpr int kVersion = 1;

    std::string                                  scene_name;
    std::unordered_map<std::string, Real>        scalars;       // R_E, S_wave, C_φ, H_E, etc.
    std::vector<Real>                            spectral;      // P(ω_k) bins
    std::vector<Real>                            spectral_freqs;
    std::unordered_map<std::string, std::string> meta;          // free-form (engine version, date)
};

// Serialize / deserialize. Throws std::runtime_error on I/O or parse error.
std::string fingerprint_to_json(Fingerprint const& fp);
Fingerprint fingerprint_from_json(std::string const& text);

void write_fingerprint(Fingerprint const& fp, std::filesystem::path const& path);
Fingerprint read_fingerprint(std::filesystem::path const& path);

// --- CSV export (one row per fingerprint) ----------------------------------
//
// Useful for batch analysis without Python — load with pandas, R, etc.
// Header: scene_name + sorted scalar keys + spectral_0..spectral_N-1.
// `append_fingerprint_csv` opens in append mode and writes a header only
// when the file is new/empty.

std::string fingerprint_csv_header(Fingerprint const& fp);
std::string fingerprint_csv_row(Fingerprint const& fp);

void append_fingerprint_csv(Fingerprint const& fp,
                            std::filesystem::path const& path);

// --- Directory loader ------------------------------------------------------
//
// Read every .fp.json under `dir` (non-recursive). Used by wavecli's
// --prototype flag to build a BinderPrototype from a directory.

std::vector<Fingerprint> load_fingerprints_dir(std::filesystem::path const& dir);

} // namespace wavelab
