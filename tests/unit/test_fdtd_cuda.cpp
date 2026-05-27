// CPU↔GPU parity for FdtdCuda<2>: identical IC, identical source,
// identical (Dirichlet) boundary on both backends → identical-within-
// roundoff field after N steps.
//
// This file is only compiled when WAVELAB_BUILD_CUDA is ON.

#include <doctest/doctest.h>

#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "fdtd/fdtd_cuda.hpp"
#include "fdtd/stability.hpp"
#include "medium/boundary.hpp"
#include "medium/medium.hpp"
#include "source/harmonic.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

TEST_CASE("FdtdCuda<2> agrees with FdtdCpuOmp<2> after 100 steps") {
    constexpr Index nx = 96;
    constexpr Index ny = 96;
    auto grid = Grid<2>::uniform(IVec<2>{nx, ny}, 0.25_r);
    Real const c0 = 1.0_r;
    Real const dt = 0.4_r * cfl_dt_max<2>(grid, c0);

    auto medium_cpu = Medium<2>::uniform(grid, c0);
    auto medium_gpu = Medium<2>::uniform(grid, c0);

    FdtdCpuOmp<2> cpu(grid, std::move(medium_cpu), dt);
    FdtdCuda<2>   gpu(grid, std::move(medium_gpu), dt);

    cpu.set_boundary(std::make_shared<Dirichlet2D>());
    gpu.set_boundary(std::make_shared<Dirichlet2D>());

    auto src_cpu = std::make_shared<HarmonicSource<2>>(
        IVec<2>{nx / 2, ny / 2}, /*amp=*/1.0_r, /*freq=*/2.0_r);
    auto src_gpu = std::make_shared<HarmonicSource<2>>(
        IVec<2>{nx / 2, ny / 2}, 1.0_r, 2.0_r);

    cpu.add_source(src_cpu);
    gpu.add_source(src_gpu);

    cpu.run(100);
    gpu.run(100);

    // Compare fields cell-by-cell.
    auto const& uc = cpu.current();
    auto const& ug = gpu.current();
    REQUIRE(uc.size() == ug.size());

    Real max_abs = 0.0_r;
    Real max_diff = 0.0_r;
    for (Index k = 0; k < uc.size(); ++k) {
        max_abs  = std::max(max_abs,  std::abs(uc.data()[k]));
        max_diff = std::max(max_diff, std::abs(uc.data()[k] - ug.data()[k]));
    }
    REQUIRE(max_abs > 0.0_r);

    Real const rel = max_diff / max_abs;
    CAPTURE(max_abs);
    CAPTURE(max_diff);
    CAPTURE(rel);
    // Float roundoff + non-associative reductions in source download/upload
    // path; both should agree to a few parts per thousand.
    CHECK(rel < 0.005_r);
}

TEST_CASE("FdtdCuda<2>: CFL violation throws") {
    auto grid = Grid<2>::uniform(IVec<2>{32, 32}, 0.1_r);
    CHECK_THROWS_AS(FdtdCuda<2>(grid, /*c=*/1.0_r, /*dt=*/0.5_r),
                    std::invalid_argument);
}

TEST_CASE("FdtdCuda<2>: time and step_count advance correctly") {
    auto grid = Grid<2>::uniform(IVec<2>{32, 32}, 0.1_r);
    Real const c0 = 1.0_r;
    Real const dt = 0.4_r * cfl_dt_max<2>(grid, c0);
    FdtdCuda<2> sim(grid, c0, dt);

    CHECK(sim.step_count() == 0);
    sim.run(50);
    CHECK(sim.step_count() == 50);
    CHECK(sim.time() == doctest::Approx(static_cast<double>(50 * dt)).epsilon(1e-5));
}

#include "fdtd/fdtd_cuda.hpp"

TEST_CASE("FdtdCuda<3> agrees with FdtdCpuOmp<3> after 50 steps") {
    constexpr Index N = 32;
    auto grid = Grid<3>::uniform(IVec<3>{N, N, N}, 0.25_r);
    Real const c0 = 1.0_r;
    Real const dt = 0.4_r * cfl_dt_max<3>(grid, c0);

    auto medium_cpu = Medium<3>::uniform(grid, c0);
    auto medium_gpu = Medium<3>::uniform(grid, c0);
    FdtdCpuOmp<3> cpu(grid, std::move(medium_cpu), dt);
    FdtdCuda<3>   gpu(grid, std::move(medium_gpu), dt);

    cpu.set_boundary(std::make_shared<Dirichlet3D>());
    gpu.set_boundary(std::make_shared<Dirichlet3D>());

    cpu.add_source(std::make_shared<HarmonicSource<3>>(
        IVec<3>{N / 2, N / 2, N / 2}, 1.0_r, 2.0_r));
    gpu.add_source(std::make_shared<HarmonicSource<3>>(
        IVec<3>{N / 2, N / 2, N / 2}, 1.0_r, 2.0_r));

    cpu.run(50);
    gpu.run(50);

    auto const& uc = cpu.current();
    auto const& ug = gpu.current();
    REQUIRE(uc.size() == ug.size());

    Real max_abs = 0.0_r;
    Real max_diff = 0.0_r;
    for (Index k = 0; k < uc.size(); ++k) {
        max_abs  = std::max(max_abs,  std::abs(uc.data()[k]));
        max_diff = std::max(max_diff, std::abs(uc.data()[k] - ug.data()[k]));
    }
    REQUIRE(max_abs > 0.0_r);
    Real const rel = max_diff / max_abs;
    CAPTURE(max_abs);
    CAPTURE(max_diff);
    CAPTURE(rel);
    CHECK(rel < 0.005_r);
}

TEST_CASE("FdtdCuda<3>: CFL violation throws") {
    auto grid = Grid<3>::uniform(IVec<3>{16, 16, 16}, 0.1_r);
    CHECK_THROWS_AS(FdtdCuda<3>(grid, /*c=*/1.0_r, /*dt=*/0.5_r),
                    std::invalid_argument);
}

} // namespace
