// wavecli — headless wavelab pipeline runner.
//
//   wavecli --synthetic <kind>    [-o out.fp.json] [options]   (2D only)
//   wavecli --pdb <file.pdb>      [-o out.fp.json] [options]
//   wavecli --sdf <file.sdf>      [-o out.fp.json] [options]
//
// <kind> in: single | dumbbell | pocket | slab     (2D synthetic only)
//
// Common options:
//   --dim {2|3}   simulation dimensionality (default 2)
//   --nx N        grid cells x (default 200)
//   --ny N        grid cells y (default 200)
//   --nz N        grid cells z (default 200; used when --dim 3)
//   --dx F        cell spacing in Å (default 0.5; uniform across axes)
//   --freq F      source frequency (default 2.0)
//   --steps N     FDTD steps (default 400)
//   --probe I,J[,K]   probe cell for spectrum
//                     (default: PAST scatterer if --add-sdf, else nx/4,ny/2)
//   --slice-z Z   z-pos for PDB midplane slice (2D only)
//   --slice-t T   slice thickness in Å (2D only, default 4.0)
//   --pulse       broadband Gaussian-pulse source (recommended for spectral
//                 fingerprinting); default is harmonic
//   --pulse-sigma S  Gaussian pulse sigma_t (default 0.6)
//   --beta-rho B  density-channel refractive-index weight (default 0.3)
//   --beta-q B    polar/charge channel weight (default 0 = off)
//   --beta-h B    hydrophobic channel weight (default 0 = off)
//   --add-sdf F   load additional ligand atoms from SDF and merge (§15)
//   --add-pdb F   load additional atoms from PDB and merge
//   --place-at X,Y,Z  center the added atoms at this world coord
//                     (default: primary scene centroid)
//   --include-hetatm  include HETATM records when parsing --pdb
//   --prototype <dir> compute binder_correlation against avg of .fp.json
//                     fingerprints in this directory
//   --csv <path>      append a CSV row of this fingerprint to <path>
//
// Output: a .fp.json fingerprint with scalar scores and the spectral
// bin vector. Determinism: same args + same input → identical scalars
// to within float-roundoff (OpenMP reduction ordering, see Phase 4 test).

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "fdtd/stability.hpp"
#include "io/fingerprint.hpp"
#include "medium/boundary.hpp"
#include "medium/build_medium.hpp"
#include "medium/medium.hpp"
#include "medium/pml.hpp"
#include "molecule/field_builder.hpp"
#include "molecule/parser_pdb.hpp"
#include "molecule/parser_sdf.hpp"
#include "molecule/scene.hpp"
#include "molecule/slice.hpp"
#include "molecule/synthetic.hpp"
#include "score/calibration.hpp"
#include "score/energy.hpp"
#include "score/entropy.hpp"
#include "score/scattering.hpp"
#include "score/similarity.hpp"
#include "score/spectral.hpp"
#include "source/gaussian_pulse.hpp"
#include "source/harmonic.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

struct Args {
    enum class Mode { Synthetic, Pdb, Sdf, None };
    Mode        mode      = Mode::None;
    std::string subject;                   // kind or path
    std::string output    = "out.fp.json";
    std::string prototype_dir;             // optional: dir of .fp.json
    std::string csv_path;                  // optional: append CSV row

    int         dim       = 2;             // 2 or 3
    Index       nx        = 200;
    Index       ny        = 200;
    Index       nz        = 200;
    Real        dx        = 0.5_r;
    Real        freq      = 2.0_r;
    Index       steps     = 400;
    Index       probe_i   = -1;
    Index       probe_j   = -1;
    Index       probe_k   = -1;            // 3D only
    Real        slice_z   = 0.0_r;
    bool        slice_z_set = false;
    Real        slice_t   = 4.0_r;
    bool        pulse_source = false;
    Real        pulse_sigma_t = 0.6_r;
    Real        beta_rho  = 0.3_r;
    Real        beta_q    = 0.0_r;
    Real        beta_h    = 0.0_r;

    // §15 pocket+ligand support.
    std::string add_sdf_path;
    std::string add_pdb_path;
    bool        place_at_set = false;
    Vec<3>      place_at{0.0_r, 0.0_r, 0.0_r};
    bool        include_hetatm = false;

    // §15 regional R_E (sum energy in a box around the ligand placement
    // instead of over the whole domain — discriminates ligand
    // perturbation from pocket-dominated total energy).
    bool        region_set    = false;
    Vec<3>      region_center{0.0_r, 0.0_r, 0.0_r};
    Real        region_half   = 5.0_r;     // half-width in Å
};

void print_usage() {
    std::puts("wavecli — wavelab headless pipeline");
    std::puts("Usage:");
    std::puts("  wavecli --synthetic <single|dumbbell|pocket|slab> [options]   (2D only)");
    std::puts("  wavecli --pdb <file.pdb> [options]");
    std::puts("  wavecli --sdf <file.sdf> [options]");
    std::puts("");
    std::puts("Common options:");
    std::puts("  -o <path>        output fingerprint path (default out.fp.json)");
    std::puts("  --dim {2|3}      simulation dimensionality (default 2)");
    std::puts("  --nx N           grid cells x (default 200)");
    std::puts("  --ny N           grid cells y (default 200)");
    std::puts("  --nz N           grid cells z (default 200; used when --dim 3)");
    std::puts("  --dx F           spacing (Å, uniform; default 0.5)");
    std::puts("  --freq F         source frequency (default 2.0)");
    std::puts("  --steps N        FDTD steps (default 400)");
    std::puts("  --probe I,J[,K]  probe cell for spectrum");
    std::puts("  --slice-z Z      z-pos for PDB slice (2D only)");
    std::puts("  --slice-t T      slice thickness (Å, 2D only, default 4.0)");
    std::puts("  --pulse          broadband Gaussian-pulse source");
    std::puts("  --pulse-sigma S  Gaussian pulse sigma_t (default 0.6)");
    std::puts("  --beta-rho B     density-channel weight (default 0.3)");
    std::puts("  --beta-q B       polar/charge channel weight (default 0)");
    std::puts("  --beta-h B       hydrophobic channel weight (default 0)");
    std::puts("  --add-sdf F      load + merge additional ligand atoms (§15)");
    std::puts("  --add-pdb F      load + merge additional atoms from PDB");
    std::puts("  --place-at X,Y,Z center added atoms at world coord");
    std::puts("  --include-hetatm include HETATM records when parsing --pdb");
    std::puts("  --region X,Y,Z,R world-coord center + half-width (Å) for");
    std::puts("                   regional_energy (§15 regional R_E).");
    std::puts("                   Auto-derived from placement when --add-sdf is");
    std::puts("                   given (defaults to half-width 5 Å around place_at).");
    std::puts("  --prototype <dir> compute binder_correlation vs prototype");
    std::puts("  --csv <path>     append CSV row to <path>");
}

// --- Parsing ---------------------------------------------------------------

bool parse_probe(char const* v, Args& a) {
    // Accepts "I,J" or "I,J,K". Sets probe_i/j/(k).
    char const* c1 = std::strchr(v, ',');
    if (!c1) return false;
    a.probe_i = std::atoll(std::string(v, c1).c_str());
    char const* c2 = std::strchr(c1 + 1, ',');
    if (!c2) {
        a.probe_j = std::atoll(c1 + 1);
        a.probe_k = -1;
    } else {
        a.probe_j = std::atoll(std::string(c1 + 1, c2).c_str());
        a.probe_k = std::atoll(c2 + 1);
    }
    return true;
}

std::optional<Args> parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string_view s{argv[i]};
        auto need = [&](char const* name) -> char const* {
            if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", name); std::exit(2); }
            return argv[++i];
        };
        if (s == "-h" || s == "--help") { print_usage(); return std::nullopt; }
        else if (s == "--synthetic") { a.mode = Args::Mode::Synthetic; a.subject = need("--synthetic"); }
        else if (s == "--pdb")       { a.mode = Args::Mode::Pdb;       a.subject = need("--pdb"); }
        else if (s == "--sdf")       { a.mode = Args::Mode::Sdf;       a.subject = need("--sdf"); }
        else if (s == "-o")          { a.output = need("-o"); }
        else if (s == "--dim")       { a.dim = std::atoi(need("--dim")); }
        else if (s == "--nx")        { a.nx = std::atoll(need("--nx")); }
        else if (s == "--ny")        { a.ny = std::atoll(need("--ny")); }
        else if (s == "--nz")        { a.nz = std::atoll(need("--nz")); }
        else if (s == "--dx")        { a.dx = static_cast<Real>(std::atof(need("--dx"))); }
        else if (s == "--freq")      { a.freq = static_cast<Real>(std::atof(need("--freq"))); }
        else if (s == "--steps")     { a.steps = std::atoll(need("--steps")); }
        else if (s == "--probe")     {
            if (!parse_probe(need("--probe"), a)) {
                std::fprintf(stderr, "--probe expects I,J or I,J,K\n");
                std::exit(2);
            }
        }
        else if (s == "--slice-z")   { a.slice_z = static_cast<Real>(std::atof(need("--slice-z"))); a.slice_z_set = true; }
        else if (s == "--slice-t")   { a.slice_t = static_cast<Real>(std::atof(need("--slice-t"))); }
        else if (s == "--prototype") { a.prototype_dir = need("--prototype"); }
        else if (s == "--csv")       { a.csv_path = need("--csv"); }
        else if (s == "--pulse")     { a.pulse_source = true; }
        else if (s == "--pulse-sigma") { a.pulse_sigma_t = static_cast<Real>(std::atof(need("--pulse-sigma"))); }
        else if (s == "--beta-rho")  { a.beta_rho = static_cast<Real>(std::atof(need("--beta-rho"))); }
        else if (s == "--beta-q")    { a.beta_q   = static_cast<Real>(std::atof(need("--beta-q"))); }
        else if (s == "--beta-h")    { a.beta_h   = static_cast<Real>(std::atof(need("--beta-h"))); }
        else if (s == "--add-sdf")   { a.add_sdf_path = need("--add-sdf"); }
        else if (s == "--add-pdb")   { a.add_pdb_path = need("--add-pdb"); }
        else if (s == "--include-hetatm") { a.include_hetatm = true; }
        else if (s == "--place-at")  {
            char const* v = need("--place-at");
            Real x, y, z;
            if (std::sscanf(v, "%g,%g,%g", &x, &y, &z) != 3) {
                std::fprintf(stderr, "--place-at expects X,Y,Z (got %s)\n", v);
                std::exit(2);
            }
            a.place_at = Vec<3>{x, y, z};
            a.place_at_set = true;
        }
        else if (s == "--region")    {
            char const* v = need("--region");
            Real x, y, z, r;
            if (std::sscanf(v, "%g,%g,%g,%g", &x, &y, &z, &r) != 4) {
                std::fprintf(stderr, "--region expects X,Y,Z,R (got %s)\n", v);
                std::exit(2);
            }
            a.region_center = Vec<3>{x, y, z};
            a.region_half   = r;
            a.region_set    = true;
        }
        else { std::fprintf(stderr, "unknown arg: %.*s\n", static_cast<int>(s.size()), s.data()); std::exit(2); }
    }
    if (a.mode == Args::Mode::None) {
        print_usage();
        std::fprintf(stderr, "\nerror: --synthetic, --pdb, or --sdf required\n");
        return std::nullopt;
    }
    if (a.dim != 2 && a.dim != 3) {
        std::fprintf(stderr, "error: --dim must be 2 or 3 (got %d)\n", a.dim);
        return std::nullopt;
    }
    return a;
}

// --- Synthetic (2D-only) ---------------------------------------------------

MolecularScene<2> build_synthetic_2d(Args const& a) {
    Real const box_x = static_cast<Real>(a.nx) * a.dx;
    Real const box_y = static_cast<Real>(a.ny) * a.dx;
    Vec<2> bmin{0.0_r, 0.0_r};
    Vec<2> bmax{box_x, box_y};
    Vec<2> center{0.5_r * box_x, 0.5_r * box_y};
    Real const sigma = 2.0_r;

    if (a.subject == "single")   return single_atom_2d(center, sigma, bmin, bmax);
    if (a.subject == "dumbbell") return dumbbell_2d(
        Vec<2>{center[0] - 3.0_r, center[1]},
        Vec<2>{center[0] + 3.0_r, center[1]},
        sigma, bmin, bmax);
    if (a.subject == "pocket")   return pocket_arc_2d(
        center, /*R=*/5.0_r, sigma, /*count=*/24,
        /*opening_dir=*/0.0_r, /*opening_rad=*/1.0_r, bmin, bmax);
    if (a.subject == "slab")     return slab_atoms_2d(
        center[0] - 1.0_r, center[0] + 1.0_r,
        /*spacing=*/1.5_r, sigma, bmin, bmax);
    std::fprintf(stderr, "unknown synthetic kind: %s\n", a.subject.c_str());
    std::exit(2);
}

// --- 3D scene loading: --pdb / --sdf primary + optional --add-sdf/--add-pdb ---

MolecularScene<3> build_scene_3d(Args const& a) {
    MolecularScene<3> scene3d;
    if (a.mode == Args::Mode::Pdb) {
        PdbParseOptions opts;
        opts.include_hetatm = a.include_hetatm;
        auto r3 = parse_pdb_file(a.subject, opts);
        if (r3.scene.atoms.empty()) {
            std::fprintf(stderr, "no atoms parsed from %s\n", a.subject.c_str());
            std::exit(2);
        }
        scene3d = std::move(r3.scene);
    } else if (a.mode == Args::Mode::Sdf) {
        auto r3 = parse_sdf_file(a.subject);
        if (r3.scene.atoms.empty()) {
            std::fprintf(stderr, "no atoms parsed from %s\n", a.subject.c_str());
            std::exit(2);
        }
        scene3d = std::move(r3.scene);
    }

    // Append translated secondary atom set(s) (§15).
    auto append_translated = [&](MolecularScene<3>& extra) {
        if (extra.atoms.empty()) return;
        Vec<3> ec{0, 0, 0};
        for (auto const& at : extra.atoms)
            for (std::size_t d = 0; d < 3; ++d) ec[d] += at.pos[d];
        Real const inv = Real{1} / static_cast<Real>(extra.atoms.size());
        for (std::size_t d = 0; d < 3; ++d) ec[d] *= inv;

        Vec<3> target = a.place_at;
        if (!a.place_at_set) {
            Vec<3> pc{0, 0, 0};
            for (auto const& at : scene3d.atoms)
                for (std::size_t d = 0; d < 3; ++d) pc[d] += at.pos[d];
            Real const pinv = scene3d.atoms.empty()
                ? Real{0}
                : Real{1} / static_cast<Real>(scene3d.atoms.size());
            for (std::size_t d = 0; d < 3; ++d) target[d] = pc[d] * pinv;
        }
        Vec<3> shift{target[0] - ec[0], target[1] - ec[1], target[2] - ec[2]};
        for (auto& at : extra.atoms)
            for (std::size_t d = 0; d < 3; ++d) at.pos[d] += shift[d];
        scene3d.atoms.insert(scene3d.atoms.end(),
                             extra.atoms.begin(), extra.atoms.end());
    };
    if (!a.add_sdf_path.empty()) {
        auto r = parse_sdf_file(a.add_sdf_path);
        append_translated(r.scene);
    }
    if (!a.add_pdb_path.empty()) {
        PdbParseOptions opts;
        opts.include_hetatm = a.include_hetatm;
        auto r = parse_pdb_file(a.add_pdb_path, opts);
        append_translated(r.scene);
    }

    // Re-derive box around all atoms with a 5Å pad.
    if (!scene3d.atoms.empty()) {
        Vec<3> mn = scene3d.atoms.front().pos;
        Vec<3> mx = mn;
        for (auto const& at : scene3d.atoms) {
            for (std::size_t d = 0; d < 3; ++d) {
                if (at.pos[d] < mn[d]) mn[d] = at.pos[d];
                if (at.pos[d] > mx[d]) mx[d] = at.pos[d];
            }
        }
        Real const pad = 5.0_r;
        for (std::size_t d = 0; d < 3; ++d) {
            scene3d.box_min[d] = mn[d] - pad;
            scene3d.box_max[d] = mx[d] + pad;
        }
    }
    return scene3d;
}

// --- D-specific helpers ---------------------------------------------------

template <int D>
Grid<D> make_grid(Args const& a, Vec<D> origin) {
    if constexpr (D == 2) {
        return Grid<2>{IVec<2>{a.nx, a.ny}, Vec<2>{a.dx, a.dx}, origin};
    } else {
        return Grid<3>{IVec<3>{a.nx, a.ny, a.nz},
                       Vec<3>{a.dx, a.dx, a.dx}, origin};
    }
}

template <int D>
IVec<D> default_source_loc(Args const& a) {
    // 2D: source on left edge, mid-y (centered protein is far enough from
    // left edge that mid-y is fine).
    // 3D: source in a grid CORNER, well outside protein bbox. Centered
    // 3D defaults like (25, ny/2, nz/2) land near the protein center,
    // which puts the source inside the regional probe and washes out
    // ligand-induced perturbation.
    if constexpr (D == 2) return IVec<2>{Index{25}, a.ny / 2};
    else                  return IVec<3>{Index{25}, Index{25}, Index{25}};
}

template <int D>
IVec<D> default_probe_loc(Args const& a) {
    // If user added a ligand, probe at scene mid (closer to placement).
    // Otherwise probe past the scatterer (3*nx/4).
    bool const has_ligand = !a.add_sdf_path.empty() || !a.add_pdb_path.empty();
    Index const px = has_ligand ? a.nx / 2 : 3 * a.nx / 4;
    if constexpr (D == 2) return IVec<2>{px, a.ny / 2};
    else                  return IVec<3>{px, a.ny / 2, a.nz / 2};
}

template <int D>
IVec<D> resolve_probe(Args const& a) {
    if (a.probe_i < 0 && a.probe_j < 0) return default_probe_loc<D>(a);
    if constexpr (D == 2) {
        return IVec<2>{
            a.probe_i < 0 ? a.nx / 2 : a.probe_i,
            a.probe_j < 0 ? a.ny / 2 : a.probe_j
        };
    } else {
        return IVec<3>{
            a.probe_i < 0 ? a.nx / 2 : a.probe_i,
            a.probe_j < 0 ? a.ny / 2 : a.probe_j,
            a.probe_k < 0 ? a.nz / 2 : a.probe_k
        };
    }
}

// --- Pipeline body, templated on D ---------------------------------------

template <int D>
int run_pipeline(Args const& a, MolecularScene<D> scene) {
    if (scene.atoms.empty()) {
        std::fprintf(stderr, "scene has no atoms\n");
        return 2;
    }

    Grid<D> grid = make_grid<D>(a, scene.box_min);
    Real const c0 = 1.0_r;

    // Splat density + optional charge / hydrophobicity fields, assemble n(x).
    Field<Real, D> rho(grid, 0.0_r);
    splat_density(rho, scene);

    std::unique_ptr<Field<Real, D>> charge_field;
    std::unique_ptr<Field<Real, D>> hydro_field;
    if (a.beta_q != Real{0}) {
        charge_field = std::make_unique<Field<Real, D>>(grid, 0.0_r);
        splat_charge(*charge_field, scene);
    }
    if (a.beta_h != Real{0}) {
        hydro_field = std::make_unique<Field<Real, D>>(grid, 0.0_r);
        splat_hydro(*hydro_field, scene);
    }

    auto medium = Medium<D>::uniform(grid, c0);
    MediumWeights w;
    w.beta_rho = a.beta_rho;
    w.beta_q   = a.beta_q;
    w.beta_h   = a.beta_h;
    build_medium_from_fields<D>(medium, &rho,
                                charge_field.get(), hydro_field.get(),
                                c0, w);
    apply_pml(medium, PmlSpec{/*cells=*/20, /*alpha_max_factor=*/2.0_r,
                              /*polynomial_order=*/3});

    Real const dt = 0.4_r * cfl_dt_max<D>(grid, c0);
    FdtdCpuOmp<D> sim(grid, std::move(medium), dt);
    sim.set_boundary(make_dirichlet<D>());

    IVec<D> src = default_source_loc<D>(a);
    if (a.pulse_source) {
        sim.add_source(std::make_shared<GaussianPulse<D>>(
            src, /*amp=*/5.0_r, a.freq,
            /*t0=*/3.0_r * a.pulse_sigma_t, a.pulse_sigma_t));
    } else {
        sim.add_source(std::make_shared<HarmonicSource<D>>(
            src, /*amp=*/1.0_r, a.freq));
    }

    IVec<D> probe = resolve_probe<D>(a);
    std::vector<Real> series;
    series.reserve(static_cast<std::size_t>(a.steps));
    for (Index s = 0; s < a.steps; ++s) {
        sim.step();
        probe_record<D>(sim.current(), probe, series);
    }

    Fingerprint fp;
    fp.scene_name = a.subject;
    fp.scalars["total_energy"] = total_energy<D>(
        sim.current(), sim.previous(), sim.wave_speed(), dt);

    // Regional energy (§15 regional R_E support). Active when either
    // --region is supplied OR the user added a ligand (in which case
    // the region defaults to a box around the placement point).
    bool has_region = a.region_set
                   || !a.add_sdf_path.empty()
                   || !a.add_pdb_path.empty();
    if constexpr (D == 2 || D == 3) {
        if (has_region) {
            Vec<3> rc = a.region_set ? a.region_center : a.place_at;
            // If not user-set, also fall back to scene centroid when
            // no place_at was given.
            if (!a.region_set && !a.place_at_set && !scene.atoms.empty()) {
                Vec<D> pc{};
                for (auto const& at : scene.atoms)
                    for (std::size_t d = 0; d < static_cast<std::size_t>(D); ++d)
                        pc[d] += at.pos[d];
                Real const inv = Real{1} / static_cast<Real>(scene.atoms.size());
                for (std::size_t d = 0; d < static_cast<std::size_t>(D); ++d) {
                    rc[d] = pc[d] * inv;
                }
            }
            // World → cell index (round, then clamp via make_probe_region).
            IVec<D> center{};
            for (std::size_t d = 0; d < static_cast<std::size_t>(D); ++d) {
                center[d] = static_cast<Index>(std::round(
                    (rc[d] - grid.origin[d]) / grid.spacing[d]));
            }
            Index const half_cells = static_cast<Index>(std::ceil(a.region_half / grid.spacing[0]));
            auto region = make_probe_region<D>(grid, center, half_cells);
            fp.scalars["regional_energy"] = region.energy(
                sim.current(), sim.previous(), sim.wave_speed(), dt);
            fp.meta["region_center_x"] = std::to_string(static_cast<double>(rc[0]));
            fp.meta["region_center_y"] = std::to_string(static_cast<double>(rc[1]));
            if constexpr (D == 3) fp.meta["region_center_z"] = std::to_string(static_cast<double>(rc[2]));
            fp.meta["region_half_A"]   = std::to_string(static_cast<double>(a.region_half));
            fp.meta["region_half_cells"] = std::to_string(half_cells);
        }
    }

    Field<Real, D> e_field(grid, 0.0_r);
    energy_density_field<D>(e_field, sim.current(), sim.previous(),
                            sim.wave_speed(), dt);
    fp.scalars["entropy"]     = energy_entropy<D>(e_field);
    fp.scalars["focus_score"] = focus_score<D>(e_field);

    auto power = compute_power_spectrum(series);
    auto freqs = compute_frequency_axis(series.size(), dt);
    auto bins  = spectral_fingerprint_logbins(power, freqs, /*nbands=*/16,
                                              /*f_min=*/0.1_r, /*f_max=*/10.0_r);
    fp.spectral       = bins;
    fp.spectral_freqs = std::vector<Real>(bins.size(), 0.0_r);
    if (!bins.empty()) {
        Real const log_lo = std::log(0.1_r);
        Real const log_hi = std::log(10.0_r);
        Real const dlog   = (log_hi - log_lo) / static_cast<Real>(bins.size());
        for (std::size_t i = 0; i < bins.size(); ++i) {
            fp.spectral_freqs[i] = std::exp(log_lo + (static_cast<Real>(i) + 0.5_r) * dlog);
        }
    }

    fp.meta["wavelab_version"] = "0.0.1";
    fp.meta["dim"]             = std::to_string(D);
    fp.meta["grid_nx"]         = std::to_string(a.nx);
    fp.meta["grid_ny"]         = std::to_string(a.ny);
    if constexpr (D == 3) fp.meta["grid_nz"] = std::to_string(a.nz);
    fp.meta["dx"]              = std::to_string(static_cast<double>(a.dx));
    fp.meta["freq"]            = std::to_string(static_cast<double>(a.freq));
    fp.meta["steps"]           = std::to_string(a.steps);
    fp.meta["atoms"]           = std::to_string(scene.atoms.size());
    fp.meta["mode"]            = (a.mode == Args::Mode::Pdb) ? "pdb"
                               : (a.mode == Args::Mode::Sdf) ? "sdf"
                               : "synthetic";

    if (!a.prototype_dir.empty()) {
        try {
            auto const binders = load_fingerprints_dir(a.prototype_dir);
            if (binders.empty()) {
                std::fprintf(stderr,
                    "warning: --prototype dir %s contained no .fp.json files\n",
                    a.prototype_dir.c_str());
            } else {
                auto proto = build_prototype(binders);
                Real const r = correlate_with_prototype(fp, proto);
                fp.scalars["binder_correlation"] = r;
                fp.meta["prototype_count"] = std::to_string(proto.sample_count);
            }
        } catch (std::exception const& ex) {
            std::fprintf(stderr, "prototype loading failed: %s\n", ex.what());
        }
    }

    write_fingerprint(fp, a.output);
    std::printf("wrote fingerprint: %s (D=%d, %zu scalars, %zu spectral bins)\n",
                a.output.c_str(), D, fp.scalars.size(), fp.spectral.size());

    if (!a.csv_path.empty()) {
        try {
            append_fingerprint_csv(fp, a.csv_path);
            std::printf("appended CSV row: %s\n", a.csv_path.c_str());
        } catch (std::exception const& ex) {
            std::fprintf(stderr, "csv append failed: %s\n", ex.what());
        }
    }
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    auto opt = parse_args(argc, argv);
    if (!opt) return 0;
    Args const a = *opt;

    if (a.mode == Args::Mode::Synthetic) {
        if (a.dim != 2) {
            std::fprintf(stderr, "error: --synthetic currently 2D-only (use --dim 2)\n");
            return 2;
        }
        return run_pipeline<2>(a, build_synthetic_2d(a));
    }

    auto scene3d = build_scene_3d(a);
    if (scene3d.atoms.empty()) {
        std::fprintf(stderr, "scene has no atoms after loading\n");
        return 2;
    }

    if (a.dim == 2) {
        Real const z_mid = a.slice_z_set
            ? a.slice_z
            : 0.5_r * (scene3d.box_min[2] + scene3d.box_max[2]);
        return run_pipeline<2>(a, slice_scene_xy_centered(scene3d, z_mid, a.slice_t));
    }
    return run_pipeline<3>(a, std::move(scene3d));
}
