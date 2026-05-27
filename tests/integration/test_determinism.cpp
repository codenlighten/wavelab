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
#include "molecule/synthetic.hpp"
#include "score/energy.hpp"
#include "score/entropy.hpp"
#include "score/spectral.hpp"
#include "source/harmonic.hpp"

#include <memory>
#include <vector>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

// Reproducible pipeline: build scene, sim, and fingerprint with fixed
// inputs. Returns the serialized JSON. Two invocations with the same
// parameters must produce byte-identical output (no clock-derived noise,
// no unordered iteration in serializer).
std::string run_pipeline() {
    constexpr Index nx = 80;
    constexpr Index ny = 80;
    constexpr Real  dx = 0.5_r;
    constexpr Real  c0 = 1.0_r;
    constexpr Real  freq = 2.0_r;
    constexpr Index steps = 150;

    Real const box_x = static_cast<Real>(nx) * dx;
    Real const box_y = static_cast<Real>(ny) * dx;
    Vec<2> bmin{0.0_r, 0.0_r};
    Vec<2> bmax{box_x, box_y};
    Vec<2> center{0.5_r * box_x, 0.5_r * box_y};
    auto scene = dumbbell_2d(
        Vec<2>{center[0] - 3.0_r, center[1]},
        Vec<2>{center[0] + 3.0_r, center[1]},
        /*sigma=*/2.0_r, bmin, bmax);

    Grid<2> grid{IVec<2>{nx, ny}, Vec<2>{dx, dx}, bmin};

    Field<Real, 2> rho(grid, 0.0_r);
    splat_density(rho, scene);

    auto medium = Medium<2>::uniform(grid, c0);
    MediumWeights w;
    w.beta_rho = 0.3_r;
    build_medium_from_fields<2>(medium, &rho, nullptr, nullptr, c0, w);
    apply_pml(medium, PmlSpec{/*cells=*/15, /*alpha_max_factor=*/2.0_r,
                              /*polynomial_order=*/3});

    Real const dt = 0.4_r * cfl_dt_max<2>(grid, c0);
    FdtdCpuOmp<2> sim(grid, std::move(medium), dt);
    sim.set_boundary(std::make_shared<Dirichlet2D>());
    sim.add_source(std::make_shared<HarmonicSource<2>>(
        IVec<2>{20, ny / 2}, /*amp=*/1.0_r, freq));

    IVec<2> probe{nx / 4, ny / 2};
    std::vector<Real> series;
    series.reserve(static_cast<std::size_t>(steps));
    for (Index s = 0; s < steps; ++s) {
        sim.step();
        probe_record<2>(sim.current(), probe, series);
    }

    Fingerprint fp;
    fp.scene_name = "determinism";
    fp.scalars["total_energy"] = total_energy_2d(
        sim.current(), sim.previous(), sim.wave_speed(), dt);

    Field<Real, 2> e_field(grid, 0.0_r);
    energy_density_field_2d(e_field, sim.current(), sim.previous(),
                            sim.wave_speed(), dt);
    fp.scalars["entropy"]     = energy_entropy(e_field);
    fp.scalars["focus_score"] = focus_score(e_field);

    auto power = compute_power_spectrum(series);
    auto freqs = compute_frequency_axis(series.size(), dt);
    fp.spectral = spectral_fingerprint_logbins(power, freqs, /*nbands=*/8,
                                               /*f_min=*/0.1_r, /*f_max=*/10.0_r);
    fp.spectral_freqs = std::vector<Real>(fp.spectral.size(), 0.0_r);
    fp.meta["fixed"] = "yes";

    return fingerprint_to_json(fp);
}

TEST_CASE("wavecli pipeline is deterministic across runs") {
    // The engine is deterministic at the float-arithmetic level, but
    // OpenMP `reduction(+:sum)` does not promise a fixed summation
    // order across threads. So two runs produce values agreeing to
    // ~1e-6, not bit-identity. The fingerprint comparison reflects this.
    std::string const a = run_pipeline();
    std::string const b = run_pipeline();
    CHECK(a.size() > 200u);

    Fingerprint const fa = fingerprint_from_json(a);
    Fingerprint const fb = fingerprint_from_json(b);

    REQUIRE(fa.scene_name == fb.scene_name);
    REQUIRE(fa.scalars.size() == fb.scalars.size());
    for (auto const& [k, v] : fa.scalars) {
        REQUIRE(fb.scalars.count(k) == 1u);
        CAPTURE(k);
        CHECK(fb.scalars.at(k) == doctest::Approx(static_cast<double>(v)).epsilon(1e-5));
    }
    REQUIRE(fa.spectral.size() == fb.spectral.size());
    for (std::size_t i = 0; i < fa.spectral.size(); ++i) {
        CAPTURE(i);
        CHECK(fb.spectral[i] == doctest::Approx(static_cast<double>(fa.spectral[i])).epsilon(1e-5));
    }
    CHECK(fa.meta == fb.meta);
}

} // namespace
