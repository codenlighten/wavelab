#include <doctest/doctest.h>

#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "medium/boundary.hpp"
#include "score/energy.hpp"

#include <cmath>
#include <memory>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

// Initialize u_curr and u_prev with the same Gaussian bump. Zero initial
// velocity (du/dt = 0 at t=0) means the bump splits into two identical
// halves traveling left and right (d'Alembert).
void seed_gaussian_split(FdtdCpuOmp<1>& sim, Real x0, Real sigma, Real amp) {
    auto const& g = sim.grid();
    Index const n = g.shape[0];
    Real const  dx = g.spacing[0];

    auto& curr = sim.current_mut();
    auto& prev = sim.previous_mut();
    for (Index i = 0; i < n; ++i) {
        Real const x  = (static_cast<Real>(i) + 0.5_r) * dx;
        Real const dx_ = (x - x0) / sigma;
        Real const v  = amp * std::exp(-0.5_r * dx_ * dx_);
        curr(i) = v;
        prev(i) = v;  // u_t(0) = 0
    }
}

// Lossless reflective box: total energy drifts < 1% over many steps.
TEST_CASE("FdtdCpuOmp<1>: energy conservation under Dirichlet, lossless") {
    constexpr Index nx = 500;
    constexpr Real  dx = 0.02_r;   // length 10
    constexpr Real  c  = 1.0_r;
    constexpr Real  dt = 0.5_r * dx / c;

    auto grid = Grid<1>::uniform(IVec<1>{nx}, dx);
    FdtdCpuOmp<1> sim(grid, c, dt);  // Dirichlet by default
    seed_gaussian_split(sim, /*x0=*/5.0_r, /*sigma=*/0.5_r, /*amp=*/1.0_r);

    // Use the leapfrog-conserved discrete energy. The centered-form
    // total_energy_1d (§14) oscillates by O(dt²) per step under leapfrog
    // — that's "breathing", not real drift. The conserved diagnostic
    // stays flat to float precision.
    sim.run(5);
    Real const e0 = total_energy_conserved_1d(sim.current(), sim.previous(), c, dt);
    REQUIRE(e0 > 0.0_r);

    Real e_min = e0, e_max = e0;
    constexpr Index total_steps = 10000;
    constexpr Index sample_each = 100;
    for (Index k = 0; k < total_steps / sample_each; ++k) {
        sim.run(sample_each);
        Real const e = total_energy_conserved_1d(sim.current(), sim.previous(), c, dt);
        if (e < e_min) e_min = e;
        if (e > e_max) e_max = e;
    }

    Real const drift = (e_max - e_min) / std::abs(e0);
    CAPTURE(e0);
    CAPTURE(e_min);
    CAPTURE(e_max);
    CAPTURE(drift);
    CHECK(drift < 0.01_r);   // < 1% per plan gate
}

// Mur-1 absorbing boundary: Gaussian pulse should radiate away to under
// 5% residual energy after the wave fronts have time to exit.
TEST_CASE("FdtdCpuOmp<1>: Mur-1 absorbs a Gaussian pulse") {
    constexpr Index nx = 1000;
    constexpr Real  dx = 0.02_r;     // length 20
    constexpr Real  c  = 1.0_r;
    constexpr Real  dt = 0.5_r * dx / c;

    auto grid = Grid<1>::uniform(IVec<1>{nx}, dx);
    FdtdCpuOmp<1> sim(grid, c, dt);
    sim.set_boundary(std::make_shared<Mur1D>());
    seed_gaussian_split(sim, /*x0=*/10.0_r, /*sigma=*/0.5_r, /*amp=*/1.0_r);

    sim.run(5);
    Real const e0 = total_energy_conserved_1d(sim.current(), sim.previous(), c, dt);
    REQUIRE(e0 > 0.0_r);

    // Pulse needs to travel L/2 = 10 units to each boundary at speed c=1,
    // so 10 simulated seconds suffice. Add headroom.
    Real const t_target = 15.0_r;
    auto const steps    = static_cast<Index>(std::round(t_target / dt));
    sim.run(steps);

    Real const e_residual = total_energy_conserved_1d(sim.current(), sim.previous(), c, dt);
    Real const fraction   = std::abs(e_residual / e0);
    CAPTURE(e0);
    CAPTURE(e_residual);
    CAPTURE(fraction);
    CHECK(fraction < 0.05_r);  // < 5% per plan gate
}

// Damping (§5) actually damps the field over time.
TEST_CASE("FdtdCpuOmp<1>: damping reduces energy over time") {
    constexpr Index nx = 500;
    constexpr Real  dx = 0.02_r;
    constexpr Real  c  = 1.0_r;
    constexpr Real  dt = 0.5_r * dx / c;
    constexpr Real  gamma = 0.5_r;  // moderate damping

    auto grid = Grid<1>::uniform(IVec<1>{nx}, dx);
    FdtdCpuOmp<1> sim(grid, c, dt, gamma);
    seed_gaussian_split(sim, 5.0_r, 0.5_r, 1.0_r);
    sim.run(5);

    Real const e0 = total_energy_conserved_1d(sim.current(), sim.previous(), c, dt);
    sim.run(2000);
    Real const e_after = total_energy_conserved_1d(sim.current(), sim.previous(), c, dt);

    CAPTURE(e0);
    CAPTURE(e_after);
    CHECK(e_after < 0.5_r * e0);
}

} // namespace
