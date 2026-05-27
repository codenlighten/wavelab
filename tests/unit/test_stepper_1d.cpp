#include <doctest/doctest.h>

#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "fdtd/stability.hpp"
#include "medium/boundary.hpp"
#include "source/harmonic.hpp"

#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

// Sample a 1D plane wave traveling rightward from a hard-source boundary
// at x=0 and check it matches sin(2π f (t - x/c)) within tolerance.
TEST_CASE("FdtdCpuOmp<1>: plane wave matches analytical solution") {
    constexpr Index nx = 500;
    constexpr Real  dx = 0.02_r;   // domain length 10
    constexpr Real  c  = 1.0_r;
    constexpr Real  dt = 0.5_r * dx / c;  // λ = 0.5, well below CFL=1
    constexpr Real  freq = 0.2_r;          // wavelength = 5.0 = 250 dx
    constexpr Real  amp  = 1.0_r;

    auto grid = Grid<1>::uniform(IVec<1>{nx}, dx);
    FdtdCpuOmp<1> sim(grid, c, dt);

    sim.add_source(std::make_shared<HarmonicSource<1>>(
        IVec<1>{0}, amp, freq));
    sim.set_boundary(std::make_shared<Mur1D>());  // absorb on right

    // Run for 3 full wave periods so steady state has propagated past
    // the sample point.
    Real const period = 1.0_r / freq;          // 5.0
    Real const t_end  = 3.0_r * period;        // 15.0
    auto const steps  = static_cast<Index>(std::round(t_end / dt));
    sim.run(steps);

    // Sample at x = 1 wavelength from the source.
    Index const sample_i = 250;
    Real const  x        = static_cast<Real>(sample_i) * dx;  // 5.0
    Real const  t        = sim.time();
    Real const  omega    = static_cast<Real>(2.0 * std::numbers::pi) * freq;
    Real const  expected = amp * std::sin(omega * (t - x / c));
    Real const  actual   = sim.current()(sample_i);

    // With 250 cells/wavelength and CFL=0.5, numerical dispersion is
    // sub-percent — but float accumulates ~1e-3 over thousands of steps.
    CHECK(actual == doctest::Approx(expected).epsilon(0.05));

    // Sanity: at x = 0 we should still see the hard-driven source value.
    Real const expected_at_source = amp * std::sin(omega * t);
    CHECK(sim.current()(Index{0}) == doctest::Approx(expected_at_source).epsilon(1e-4));
}

TEST_CASE("FdtdCpuOmp<1>: CFL violation throws") {
    auto grid = Grid<1>::uniform(IVec<1>{100}, 0.01_r);
    Real const c        = 1.0_r;
    Real const dt_bad   = 0.5_r;  // way too large (CFL says dt <= 0.01/1)

    CHECK_THROWS_AS(FdtdCpuOmp<1>(grid, c, dt_bad), std::invalid_argument);
}

TEST_CASE("FdtdCpuOmp<1>: CFL boundary acceptance") {
    auto grid = Grid<1>::uniform(IVec<1>{100}, 0.01_r);
    Real const c    = 1.0_r;
    Real const dt_max = cfl_dt_max<1>(grid, c);  // 0.01 in 1D
    Real const dt_ok  = 0.9_r * dt_max;

    CHECK_NOTHROW(FdtdCpuOmp<1>(grid, c, dt_ok));
}

TEST_CASE("FdtdCpuOmp<1>: step count and time advance") {
    auto grid = Grid<1>::uniform(IVec<1>{64}, 0.1_r);
    Real const c  = 1.0_r;
    Real const dt = 0.05_r;
    FdtdCpuOmp<1> sim(grid, c, dt);

    CHECK(sim.step_count() == 0);
    CHECK(sim.time() == doctest::Approx(0.0));

    sim.run(10);

    CHECK(sim.step_count() == 10);
    CHECK(sim.time() == doctest::Approx(0.5));  // 10 * 0.05
}

} // namespace
