// wavecli — headless wavelab pipeline runner.
//
//   wavecli --synthetic <kind>    [-o out.fp.json] [options]
//   wavecli --pdb <file.pdb>      [-o out.fp.json] [options]
//
// <kind> in: single | dumbbell | pocket | slab
//
// Options:
//   --nx N       grid cells in x (default 200)
//   --ny N       grid cells in y (default 200)
//   --dx F       cell spacing in Å (default 0.5)
//   --freq F     source frequency (default 2.0)
//   --steps N    FDTD steps (default 400)
//   --probe I,J  probe cell for spectral fingerprint (default = nx/4, ny/2)
//   --slice-z Z  z-position for PDB midplane slice (default = midplane)
//   --slice-t T  slice thickness in Å (default 4.0)
//
// Output: a .fp.json fingerprint with scalar scores and the spectral
// bin vector. Run twice with the same args → byte-identical output
// (determinism check used by Phase 4 acceptance test).

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
#include "score/energy.hpp"
#include "score/entropy.hpp"
#include "score/scattering.hpp"
#include "score/similarity.hpp"
#include "score/calibration.hpp"
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
    Index       nx        = 200;
    Index       ny        = 200;
    Real        dx        = 0.5_r;
    Real        freq      = 2.0_r;
    Index       steps     = 400;
    Index       probe_i   = -1;
    Index       probe_j   = -1;
    Real        slice_z   = 0.0_r;
    bool        slice_z_set = false;
    Real        slice_t   = 4.0_r;
    bool        pulse_source = false;     // GaussianPulse instead of HarmonicSource
    Real        pulse_sigma_t = 0.6_r;
    Real        beta_rho  = 0.3_r;        // refractive-index weight
};

void print_usage() {
    std::puts("wavecli — wavelab headless pipeline");
    std::puts("Usage:");
    std::puts("  wavecli --synthetic <single|dumbbell|pocket|slab> [options]");
    std::puts("  wavecli --pdb <file.pdb> [options]");
    std::puts("  wavecli --sdf <file.sdf> [options]");
    std::puts("");
    std::puts("Options:");
    std::puts("  -o <path>     output fingerprint path (default out.fp.json)");
    std::puts("  --nx N        grid cells x (default 200)");
    std::puts("  --ny N        grid cells y (default 200)");
    std::puts("  --dx F        spacing (Å, default 0.5)");
    std::puts("  --freq F      source frequency (default 2.0)");
    std::puts("  --steps N     FDTD steps (default 400)");
    std::puts("  --probe I,J   probe cell for spectrum");
    std::puts("  --slice-z Z   z-pos for PDB slice");
    std::puts("  --slice-t T   slice thickness (Å, default 4.0)");
    std::puts("  --pulse       Gaussian-pulse source (broadband, recommended for");
    std::puts("                spectral fingerprinting); default is harmonic");
    std::puts("  --pulse-sigma S  Gaussian pulse sigma_t (default 0.6)");
    std::puts("  --beta-rho B  refractive-index weight per density unit (default 0.3)");
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
        else if (s == "--nx")        { a.nx = std::atoll(need("--nx")); }
        else if (s == "--ny")        { a.ny = std::atoll(need("--ny")); }
        else if (s == "--dx")        { a.dx = static_cast<Real>(std::atof(need("--dx"))); }
        else if (s == "--freq")      { a.freq = static_cast<Real>(std::atof(need("--freq"))); }
        else if (s == "--steps")     { a.steps = std::atoll(need("--steps")); }
        else if (s == "--probe")     {
            char const* v = need("--probe");
            char const* comma = std::strchr(v, ',');
            if (!comma) { std::fprintf(stderr, "--probe needs I,J\n"); std::exit(2); }
            a.probe_i = std::atoll(std::string(v, comma).c_str());
            a.probe_j = std::atoll(comma + 1);
        }
        else if (s == "--slice-z")   { a.slice_z = static_cast<Real>(std::atof(need("--slice-z"))); a.slice_z_set = true; }
        else if (s == "--slice-t")   { a.slice_t = static_cast<Real>(std::atof(need("--slice-t"))); }
        else if (s == "--prototype") { a.prototype_dir = need("--prototype"); }
        else if (s == "--csv")       { a.csv_path = need("--csv"); }
        else if (s == "--pulse")     { a.pulse_source = true; }
        else if (s == "--pulse-sigma") { a.pulse_sigma_t = static_cast<Real>(std::atof(need("--pulse-sigma"))); }
        else if (s == "--beta-rho")  { a.beta_rho = static_cast<Real>(std::atof(need("--beta-rho"))); }
        else { std::fprintf(stderr, "unknown arg: %.*s\n", static_cast<int>(s.size()), s.data()); std::exit(2); }
    }
    if (a.mode == Args::Mode::None) {
        print_usage();
        std::fprintf(stderr, "\nerror: --synthetic, --pdb, or --sdf required\n");
        return std::nullopt;
    }
    // probe default now set later (after we know source mode), see main().
    return a;
}

MolecularScene<2> build_scene(Args const& a) {
    if (a.mode == Args::Mode::Synthetic) {
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

    // 3D scene loaders. Both PDB and SDF produce MolecularScene<3>; we
    // slice through the midplane (or user z) to feed the 2D engine.
    MolecularScene<3> scene3d;
    if (a.mode == Args::Mode::Pdb) {
        auto r3 = parse_pdb_file(a.subject);
        if (r3.scene.atoms.empty()) {
            std::fprintf(stderr, "no atoms parsed from %s\n", a.subject.c_str());
            std::exit(2);
        }
        scene3d = std::move(r3.scene);
    } else {
        auto r3 = parse_sdf_file(a.subject);
        if (r3.scene.atoms.empty()) {
            std::fprintf(stderr, "no atoms parsed from %s\n", a.subject.c_str());
            std::exit(2);
        }
        scene3d = std::move(r3.scene);
    }
    Real const z_mid = a.slice_z_set
        ? a.slice_z
        : 0.5_r * (scene3d.box_min[2] + scene3d.box_max[2]);
    return slice_scene_xy_centered(scene3d, z_mid, a.slice_t);
}

} // namespace

int main(int argc, char** argv) {
    auto opt = parse_args(argc, argv);
    if (!opt) return 0;
    Args const a = *opt;

    auto scene = build_scene(a);
    if (scene.atoms.empty()) {
        std::fprintf(stderr, "scene has no atoms after slicing\n");
        return 2;
    }

    // Grid: just the user-specified dimensions; origin = scene box_min.
    Grid<2> grid{IVec<2>{a.nx, a.ny}, Vec<2>{a.dx, a.dx}, scene.box_min};
    Real const c0 = 1.0_r;

    // Density splat → refractive index → wave speed.
    Field<Real, 2> rho(grid, 0.0_r);
    splat_density(rho, scene);

    auto medium = Medium<2>::uniform(grid, c0);
    MediumWeights w;
    w.beta_rho = a.beta_rho;
    build_medium_from_fields<2>(medium, &rho, nullptr, nullptr, c0, w);
    apply_pml(medium, PmlSpec{/*cells=*/20, /*alpha_max_factor=*/2.0_r,
                              /*polynomial_order=*/3});

    Real const dt = 0.4_r * cfl_dt_max<2>(grid, c0);
    FdtdCpuOmp<2> sim(grid, std::move(medium), dt);
    sim.set_boundary(std::make_shared<Dirichlet2D>());

    // Source: harmonic (default) or broadband Gaussian pulse. The pulse
    // mode gives a richer spectrum and is what discrimination/scoring
    // workflows usually want; harmonic stays for steady-state probing.
    IVec<2> src{25, a.ny / 2};
    if (a.pulse_source) {
        sim.add_source(std::make_shared<GaussianPulse<2>>(
            src, /*amp=*/5.0_r, a.freq,
            /*t0=*/3.0_r * a.pulse_sigma_t, a.pulse_sigma_t));
    } else {
        sim.add_source(std::make_shared<HarmonicSource<2>>(
            src, /*amp=*/1.0_r, a.freq));
    }

    // Spectral probe — by default sample PAST the scatterer (3*nx/4)
    // so the recorded signal carries the scattered field; user can
    // override with --probe.
    Index probe_i = a.probe_i < 0 ? 3 * a.nx / 4 : a.probe_i;
    Index probe_j = a.probe_j < 0 ? a.ny / 2     : a.probe_j;
    IVec<2> probe{probe_i, probe_j};
    std::vector<Real> series;
    series.reserve(static_cast<std::size_t>(a.steps));
    for (Index s = 0; s < a.steps; ++s) {
        sim.step();
        probe_record<2>(sim.current(), probe, series);
    }

    // Build fingerprint.
    Fingerprint fp;
    fp.scene_name = a.subject;
    fp.scalars["total_energy"] = total_energy_2d(
        sim.current(), sim.previous(), sim.wave_speed(), dt);

    // Energy field for entropy/focus
    Field<Real, 2> e_field(grid, 0.0_r);
    energy_density_field_2d(e_field, sim.current(), sim.previous(),
                            sim.wave_speed(), dt);
    fp.scalars["entropy"]      = energy_entropy(e_field);
    fp.scalars["focus_score"]  = focus_score(e_field);

    // Spectral
    auto power = compute_power_spectrum(series);
    auto freqs = compute_frequency_axis(series.size(), dt);
    auto bins  = spectral_fingerprint_logbins(power, freqs, /*nbands=*/16,
                                              /*f_min=*/0.1_r, /*f_max=*/10.0_r);
    fp.spectral       = bins;
    fp.spectral_freqs = std::vector<Real>(bins.size(), 0.0_r);
    // Mid-bin frequency (log-mid)
    if (!bins.empty()) {
        Real const log_lo = std::log(0.1_r);
        Real const log_hi = std::log(10.0_r);
        Real const dlog   = (log_hi - log_lo) / static_cast<Real>(bins.size());
        for (std::size_t i = 0; i < bins.size(); ++i) {
            fp.spectral_freqs[i] = std::exp(log_lo + (static_cast<Real>(i) + 0.5_r) * dlog);
        }
    }

    fp.meta["wavelab_version"] = "0.0.1";
    fp.meta["grid_nx"]         = std::to_string(a.nx);
    fp.meta["grid_ny"]         = std::to_string(a.ny);
    fp.meta["dx"]              = std::to_string(static_cast<double>(a.dx));
    fp.meta["freq"]            = std::to_string(static_cast<double>(a.freq));
    fp.meta["steps"]           = std::to_string(a.steps);
    fp.meta["atoms"]           = std::to_string(scene.atoms.size());
    fp.meta["mode"]            = (a.mode == Args::Mode::Pdb) ? "pdb" : "synthetic";

    // Optional: correlate against a binder prototype directory.
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
    std::printf("wrote fingerprint: %s (%zu scalars, %zu spectral bins)\n",
                a.output.c_str(), fp.scalars.size(), fp.spectral.size());

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
