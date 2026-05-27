#include <doctest/doctest.h>

#include "core/memory_budget.hpp"
#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "fdtd/stability.hpp"
#include "medium/boundary.hpp"
#include "medium/build_medium.hpp"
#include "medium/medium.hpp"
#include "medium/pml.hpp"
#include "molecule/field_builder.hpp"
#include "molecule/parser_pdb.hpp"
#include "molecule/scene.hpp"
#include "score/energy.hpp"
#include "score/entropy.hpp"
#include "source/gaussian_pulse.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <memory>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

std::filesystem::path locate_1ubq() {
    if (char const* env = std::getenv("WAVELAB_TEST_DATA")) {
        std::filesystem::path p = std::filesystem::path{env} / "1ubq.pdb";
        if (std::filesystem::exists(p)) return p;
    }
    std::filesystem::path candidates[] = {
        std::filesystem::current_path() / "data" / "pdb" / "1ubq.pdb",
        std::filesystem::current_path() / ".." / "data" / "pdb" / "1ubq.pdb",
        "/home/greg/Documents/CLAUDE_PROJECTS/light-simulation/data/pdb/1ubq.pdb",
    };
    for (auto const& p : candidates) {
        if (std::filesystem::exists(p)) return p;
    }
    return {};
}

// Full 3D 1ubq run on a grid that fits the whole protein. We don't try
// to hit the 128³ from the plan — 64³ keeps the unit-test loop fast
// while still exercising the full pipeline (load, 3D splat, 3D medium,
// 3D PML, 3D FDTD, 3D scoring).
TEST_CASE("1ubq full 3D volume runs to completion with finite scores") {
    auto path = locate_1ubq();
    if (path.empty()) {
        MESSAGE("Skipping: data/pdb/1ubq.pdb not found");
        return;
    }
    auto r = parse_pdb_file(path);
    REQUIRE(r.scene.atoms.size() > 0u);

    // Size the grid to span the protein + padding, with isotropic
    // spacing chosen to land at a 64³ box.
    Real const ext_x = r.scene.box_max[0] - r.scene.box_min[0];
    Real const ext_y = r.scene.box_max[1] - r.scene.box_min[1];
    Real const ext_z = r.scene.box_max[2] - r.scene.box_min[2];
    Real const max_ext = std::max({ext_x, ext_y, ext_z});
    Index const N = 64;
    Real const dx = max_ext / static_cast<Real>(N);

    Grid<3> grid{IVec<3>{N, N, N}, Vec<3>{dx, dx, dx}, r.scene.box_min};

    // Sanity: memory budget. 64³ × 5 fields × 4B = ~5 MB.
    auto est = estimate_engine_memory<3>(grid);
    CAPTURE(format_bytes(est.total_bytes));
    REQUIRE_NOTHROW(enforce_memory_budget<3>(grid, 200ull * 1024 * 1024));

    // Splat density into a 3D field.
    Field<Real, 3> rho(grid, 0.0_r);
    splat_density(rho, r.scene);

    Real const c0 = 1.0_r;
    auto medium = Medium<3>::uniform(grid, c0);
    MediumWeights w;
    w.beta_rho = 0.3_r;
    build_medium_from_fields<3>(medium, &rho, nullptr, nullptr, c0, w);
    apply_pml(medium, PmlSpec{/*cells=*/10, /*alpha_max_factor=*/2.0_r,
                              /*polynomial_order=*/3});

    Real const dt = 0.4_r * cfl_dt_max<3>(grid, c0);
    FdtdCpuOmp<3> sim(grid, std::move(medium), dt);
    sim.set_boundary(std::make_shared<Dirichlet3D>());

    // Pulse source near the grid center.
    sim.add_source(std::make_shared<GaussianPulse<3>>(
        IVec<3>{N / 2, N / 2, N / 2}, /*amp=*/5.0_r,
        /*freq=*/2.0_r, /*t0=*/2.0_r, /*sigma_t=*/0.5_r));

    // Run for a moderate number of steps — enough to feel the molecule.
    sim.run(100);

    Real const e_tot = total_energy_3d(sim.current(), sim.previous(),
                                       sim.wave_speed(), dt);

    Field<Real, 3> e_field(grid, 0.0_r);
    energy_density_field_3d(e_field, sim.current(), sim.previous(),
                            sim.wave_speed(), dt);
    Real const H = energy_entropy<3>(e_field);

    CAPTURE(r.atom_count);
    CAPTURE(e_tot);
    CAPTURE(H);
    CHECK(std::isfinite(e_tot));
    CHECK(e_tot > 0.0_r);
    CHECK(std::isfinite(H));
    CHECK(H > 0.0_r);
}

} // namespace
