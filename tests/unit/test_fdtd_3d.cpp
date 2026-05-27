#include <doctest/doctest.h>

#include "core/memory_budget.hpp"
#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "fdtd/stability.hpp"
#include "medium/boundary.hpp"
#include "medium/medium.hpp"
#include "medium/pml.hpp"
#include "score/energy.hpp"
#include "score/entropy.hpp"
#include "source/gaussian_pulse.hpp"
#include "source/harmonic.hpp"

#include <cmath>
#include <memory>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

TEST_CASE("FdtdCpuOmp<3>: uniform medium runs without exploding") {
    auto grid = Grid<3>::uniform(IVec<3>{32, 32, 32}, 0.5_r);
    Real const c0 = 1.0_r;
    Real const dt = 0.4_r * cfl_dt_max<3>(grid, c0);

    FdtdCpuOmp<3> sim(grid, c0, dt);
    sim.add_source(std::make_shared<HarmonicSource<3>>(
        IVec<3>{16, 16, 16}, 1.0_r, 2.0_r));
    sim.run(80);

    Real const e = total_energy_3d(sim.current(), sim.previous(), c0, dt);
    CHECK(std::isfinite(e));
    CHECK(e > 0.0_r);
    CHECK(e < 1e8_r);
}

// 3D Green's function check (§33): a point impulse / Gaussian pulse at
// the origin produces a spherically-symmetric outgoing wave with
// amplitude envelope ∝ 1/r in 3D. We verify spherical symmetry by
// comparing peak amplitudes at three equidistant points along the x,
// y, and z axes.
TEST_CASE("FdtdCpuOmp<3>: point source is spherically symmetric") {
    constexpr Index N = 64;
    auto grid = Grid<3>::uniform(IVec<3>{N, N, N}, 0.5_r);
    Real const c0 = 1.0_r;
    Real const dt = 0.4_r * cfl_dt_max<3>(grid, c0);

    auto medium = Medium<3>::uniform(grid, c0);
    apply_pml(medium, PmlSpec{/*cells=*/10, /*alpha_max_factor=*/2.0_r,
                              /*polynomial_order=*/3});

    FdtdCpuOmp<3> sim(grid, std::move(medium), dt);
    sim.set_boundary(std::make_shared<Dirichlet3D>());

    Index const cx = N / 2, cy = N / 2, cz = N / 2;
    sim.add_source(std::make_shared<GaussianPulse<3>>(
        IVec<3>{cx, cy, cz}, /*amp=*/100.0_r,
        /*freq=*/2.0_r, /*t0=*/1.5_r, /*sigma_t=*/0.4_r));

    // Track max |u| at three equidistant probes (r = 8 cells = 4 units)
    Index const r = 8;
    Real peak_x = 0.0_r, peak_y = 0.0_r, peak_z = 0.0_r;

    // Run long enough for the pulse front to reach the probes
    // (t_arrival = r·dx/c = 4) and for the peak to pass.
    Real const t_end = 10.0_r;
    auto const steps = static_cast<Index>(std::round(t_end / dt));
    for (Index s = 0; s < steps; ++s) {
        sim.step();
        peak_x = std::max(peak_x, std::abs(sim.current()(cx + r, cy, cz)));
        peak_y = std::max(peak_y, std::abs(sim.current()(cx, cy + r, cz)));
        peak_z = std::max(peak_z, std::abs(sim.current()(cx, cy, cz + r)));
    }

    CAPTURE(peak_x);
    CAPTURE(peak_y);
    CAPTURE(peak_z);
    REQUIRE(peak_x > 0.0_r);
    REQUIRE(peak_y > 0.0_r);
    REQUIRE(peak_z > 0.0_r);

    // Symmetry: all three should agree within numerical error (the grid
    // isn't isotropic at this resolution, so allow ~25% spread).
    Real const peak_mean = (peak_x + peak_y + peak_z) / 3.0_r;
    CHECK(std::abs(peak_x - peak_mean) / peak_mean < 0.25_r);
    CHECK(std::abs(peak_y - peak_mean) / peak_mean < 0.25_r);
    CHECK(std::abs(peak_z - peak_mean) / peak_mean < 0.25_r);
}

// Verify the templatized scoring functions compile and execute for D=3.
TEST_CASE("Scoring functions work for D=3") {
    auto grid = Grid<3>::uniform(IVec<3>{16, 16, 16}, 1.0_r);
    Field<Real, 3> e(grid, 1.0_r);
    Field<Real, 3> mask(grid, 0.0_r);

    Real const H     = energy_entropy<3>(e);
    Real const log_M = std::log(static_cast<Real>(e.size()));
    CHECK(H == doctest::Approx(static_cast<double>(log_M)).epsilon(1e-5));

    CHECK(focus_score<3>(e) == doctest::Approx(0.0).epsilon(1e-5));

    // Fill quarter mask, expect 1/8 of cells = 0.125 hotspot fraction
    for (Index i = 0; i < 8; ++i)
        for (Index j = 0; j < 8; ++j)
            for (Index k = 0; k < 8; ++k)
                mask(i, j, k) = 1.0_r;
    CHECK(hotspot_concentration<3>(e, mask) == doctest::Approx(0.125));
}

TEST_CASE("Memory budget estimator: catches oversized grids") {
    auto grid_ok = Grid<3>::uniform(IVec<3>{32, 32, 32}, 1.0_r);
    CHECK_NOTHROW(enforce_memory_budget<3>(grid_ok, /*cap=*/100ull * 1024 * 1024));

    auto grid_huge = Grid<3>::uniform(IVec<3>{1024, 1024, 1024}, 1.0_r);
    CHECK_THROWS_AS(enforce_memory_budget<3>(grid_huge, /*cap=*/1ull * 1024 * 1024 * 1024),
                    std::runtime_error);
}

} // namespace
