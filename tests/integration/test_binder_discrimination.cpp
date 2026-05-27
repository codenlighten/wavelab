#include <doctest/doctest.h>

#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "fdtd/stability.hpp"
#include "io/fingerprint.hpp"
#include "medium/boundary.hpp"
#include "medium/build_medium.hpp"
#include "medium/medium.hpp"
#include "medium/pml.hpp"
#include "molecule/field_builder.hpp"
#include "molecule/scene.hpp"
#include "molecule/synthetic.hpp"
#include "score/calibration.hpp"
#include "score/energy.hpp"
#include "score/entropy.hpp"
#include "score/spectral.hpp"
#include "source/gaussian_pulse.hpp"
#include "source/harmonic.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

Fingerprint simulate_scene(MolecularScene<2> const& scene,
                           Real center_freq, Index steps,
                           std::string const& name) {
    constexpr Index nx = 80;
    constexpr Index ny = 80;
    constexpr Real  dx = 0.5_r;
    constexpr Real  c0 = 1.0_r;
    Grid<2> grid{IVec<2>{nx, ny}, Vec<2>{dx, dx}, scene.box_min};

    Field<Real, 2> rho(grid, 0.0_r);
    splat_density(rho, scene);

    auto medium = Medium<2>::uniform(grid, c0);
    MediumWeights w;
    w.beta_rho = 0.3_r;
    build_medium_from_fields<2>(medium, &rho, nullptr, nullptr, c0, w);
    apply_pml(medium, PmlSpec{15, 2.0_r, 3});

    Real const dt = 0.4_r * cfl_dt_max<2>(grid, c0);
    FdtdCpuOmp<2> sim(grid, std::move(medium), dt);
    sim.set_boundary(std::make_shared<Dirichlet2D>());
    // Broadband Gaussian pulse — the spectral fingerprint then reflects
    // how the scene reshapes a wide-frequency excitation.
    sim.add_source(std::make_shared<GaussianPulse<2>>(
        IVec<2>{20, ny / 2}, /*amp=*/5.0_r, center_freq,
        /*t0=*/3.0_r, /*sigma_t=*/0.6_r));

    // Probe past the scatterer so the recorded signal encodes scattering.
    IVec<2> probe{3 * nx / 4, ny / 2};
    std::vector<Real> series;
    series.reserve(static_cast<std::size_t>(steps));
    for (Index s = 0; s < steps; ++s) {
        sim.step();
        probe_record<2>(sim.current(), probe, series);
    }

    Fingerprint fp;
    fp.scene_name = name;
    auto power = compute_power_spectrum(series);
    auto freqs = compute_frequency_axis(series.size(), dt);
    fp.spectral = spectral_fingerprint_logbins(power, freqs, /*nbands=*/16,
                                               /*f_min=*/0.1_r, /*f_max=*/10.0_r);
    fp.spectral_freqs = freqs;
    return fp;
}

// Binder family: dumbbells of varied separation (similar wave signatures).
// Non-binder: a single isolated atom — qualitatively different scattering.
TEST_CASE("Binder prototype discriminates binder from non-binder") {
    constexpr Real freq = 2.0_r;
    constexpr Index steps = 200;

    Real const box = 40.0_r;
    Vec<2> bmin{0.0_r, 0.0_r};
    Vec<2> bmax{box, box};
    Vec<2> center{0.5_r * box, 0.5_r * box};

    // Build K binder fingerprints — dumbbells at varied separation.
    std::vector<Fingerprint> binders;
    Real const seps[] = {2.5_r, 3.0_r, 3.5_r, 4.0_r, 4.5_r};
    for (Real sep : seps) {
        auto scene = dumbbell_2d(
            Vec<2>{center[0] - sep, center[1]},
            Vec<2>{center[0] + sep, center[1]},
            /*sigma=*/2.0_r, bmin, bmax);
        binders.push_back(simulate_scene(scene, freq, steps, "binder"));
    }

    auto proto = build_prototype(binders);
    REQUIRE(proto.sample_count == 5);
    REQUIRE(!proto.spectral_mean.empty());

    // Candidate 1: dumbbell-like (a new separation in the binder family).
    auto binder_like_scene = dumbbell_2d(
        Vec<2>{center[0] - 3.25_r, center[1]},
        Vec<2>{center[0] + 3.25_r, center[1]},
        2.0_r, bmin, bmax);
    auto binder_like_fp = simulate_scene(binder_like_scene, freq, steps, "candidate_binder");

    // Candidate 2: non-binder shape — a single isolated atom.
    auto non_binder_scene = single_atom_2d(center, 2.0_r, bmin, bmax);
    auto non_binder_fp = simulate_scene(non_binder_scene, freq, steps, "candidate_nonbinder");

    Real const r_binder    = correlate_with_prototype(binder_like_fp, proto);
    Real const r_nonbinder = correlate_with_prototype(non_binder_fp, proto);

    CAPTURE(r_binder);
    CAPTURE(r_nonbinder);
    CHECK(r_binder > r_nonbinder);
    CHECK(r_binder > 0.9_r);          // binder-family candidate should look very similar
}

TEST_CASE("CSV round-trip via append_fingerprint_csv") {
    Fingerprint fp;
    fp.scene_name = "csv_test";
    fp.scalars["a"] = 1.0_r;
    fp.scalars["b"] = 2.5_r;
    fp.spectral = {0.1_r, 0.2_r, 0.3_r};

    auto path = std::filesystem::temp_directory_path() / "wavelab_csv_test.csv";
    std::filesystem::remove(path);
    append_fingerprint_csv(fp, path);
    append_fingerprint_csv(fp, path);

    std::ifstream in(path);
    REQUIRE(in.good());
    std::string line;
    int lines = 0;
    while (std::getline(in, line)) ++lines;
    CHECK(lines == 3);  // header + 2 data rows

    std::filesystem::remove(path);
}

} // namespace
