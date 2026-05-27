#include <doctest/doctest.h>

#include "core/field.hpp"
#include "core/grid.hpp"
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
#include "molecule/slice.hpp"
#include "score/energy.hpp"
#include "score/similarity.hpp"
#include "source/gaussian_pulse.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <numbers>
#include <string>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

// Locate the PDB sample. Allows override via WAVELAB_TEST_DATA env var.
std::filesystem::path locate_1ubq() {
    if (char const* env = std::getenv("WAVELAB_TEST_DATA")) {
        std::filesystem::path p = std::filesystem::path{env} / "1ubq.pdb";
        if (std::filesystem::exists(p)) return p;
    }
    // Project layout: tests/ sibling to data/
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

TEST_CASE("PDB parser: 1ubq loads with 602 ATOM records") {
    auto path = locate_1ubq();
    if (path.empty()) {
        MESSAGE("Skipping 1ubq test: data/pdb/1ubq.pdb not found");
        return;
    }
    auto const r = parse_pdb_file(path);
    CAPTURE(r.atom_count);
    CAPTURE(r.hetatm_count);
    CHECK(r.atom_count == 602);
    CHECK(r.scene.atoms.size() == 602u);
    // Hetatm excluded by default — should be zero in the scene.
    CHECK(r.hetatm_count == 0);
    // Box should encompass actual coords with padding.
    CHECK(r.scene.box_max[0] > r.scene.box_min[0]);
    CHECK(r.scene.box_max[1] > r.scene.box_min[1]);
    CHECK(r.scene.box_max[2] > r.scene.box_min[2]);
}

TEST_CASE("PDB parser: 1ubq with HETATM includes waters") {
    auto path = locate_1ubq();
    if (path.empty()) return;
    PdbParseOptions opts;
    opts.include_hetatm = true;
    auto const r = parse_pdb_file(path, opts);
    CAPTURE(r.atom_count);
    CAPTURE(r.hetatm_count);
    CHECK(r.atom_count == 602);
    CHECK(r.hetatm_count == 58);
    CHECK(r.scene.atoms.size() == 660u);
}

// Density splatting integral check: for D-dimensional Gaussians of width σ,
// ∫ exp(-|x-r|² / 2σ²) dV = (2π σ²)^(D/2). Sum over atoms gives total
// integral; the discrete sum × cell volume should match within a few %
// (Gaussian cut at 4σ truncates a small tail, plus discretization).
TEST_CASE("Density splat volume integral ≈ N·(2πσ²)^(D/2)") {
    auto path = locate_1ubq();
    if (path.empty()) return;
    auto const r = parse_pdb_file(path);

    // 2D slice through midplane, 4 Å thick.
    auto scene2d = slice_scene_xy_midplane(r.scene, 4.0_r);
    REQUIRE(scene2d.atoms.size() > 0u);

    // Set up grid spanning the scene's xy extent at 0.5 Å spacing.
    Real const h = 0.5_r;
    Real const ext_x = scene2d.box_max[0] - scene2d.box_min[0];
    Real const ext_y = scene2d.box_max[1] - scene2d.box_min[1];
    Index const nx = static_cast<Index>(std::ceil(ext_x / h)) + 1;
    Index const ny = static_cast<Index>(std::ceil(ext_y / h)) + 1;

    Grid<2> grid{IVec<2>{nx, ny}, Vec<2>{h, h}, scene2d.box_min};
    Field<Real, 2> rho(grid, 0.0_r);
    splat_density(rho, scene2d);

    Real sum = 0.0_r;
    for (Index k = 0; k < rho.size(); ++k) sum += rho.data()[k];
    Real const cell_area = h * h;
    Real const integral  = sum * cell_area;

    // Expected — sum over atoms of (2π σ²) for D=2.
    Real expected = 0.0_r;
    Real const pi = static_cast<Real>(std::numbers::pi);
    for (auto const& a : scene2d.atoms) {
        expected += 2.0_r * pi * a.sigma * a.sigma;
    }

    Real const rel_err = std::abs(integral - expected) / expected;
    CAPTURE(integral);
    CAPTURE(expected);
    CAPTURE(rel_err);
    CHECK(rel_err < 0.05_r);   // 5% per plan gate
}

// End-to-end: load 1ubq, slice, splat density into a refractive-index
// field, run a 2D simulation, verify scattering is non-trivial relative
// to an empty-medium baseline.
TEST_CASE("1ubq slice produces non-trivial scattering vs empty baseline") {
    auto path = locate_1ubq();
    if (path.empty()) return;
    auto const r = parse_pdb_file(path);
    auto scene2d = slice_scene_xy_midplane(r.scene, 4.0_r);
    REQUIRE(scene2d.atoms.size() > 0u);

    // Grid covering scene's xy extent at 0.5 Å spacing.
    Real const h = 0.5_r;
    Real const ext_x = scene2d.box_max[0] - scene2d.box_min[0];
    Real const ext_y = scene2d.box_max[1] - scene2d.box_min[1];
    Index const nx = static_cast<Index>(std::ceil(ext_x / h)) + 1;
    Index const ny = static_cast<Index>(std::ceil(ext_y / h)) + 1;

    Grid<2> grid{IVec<2>{nx, ny}, Vec<2>{h, h}, scene2d.box_min};
    Real const c0 = 1.0_r;

    // --- empty-medium baseline ---
    auto medium_empty = Medium<2>::uniform(grid, c0);
    apply_pml(medium_empty, PmlSpec{/*cells=*/15, /*alpha_max_factor=*/2.0_r, /*polynomial_order=*/3});
    Real const dt = 0.4_r * cfl_dt_max<2>(grid, c0);
    FdtdCpuOmp<2> sim_empty(grid, std::move(medium_empty), dt);
    sim_empty.set_boundary(std::make_shared<Dirichlet2D>());

    // Pulse near domain center.
    IVec<2> src_loc{nx / 2, ny / 2};
    sim_empty.add_source(std::make_shared<GaussianPulse<2>>(
        src_loc, /*amp=*/5.0_r, /*freq=*/3.0_r, /*t0=*/2.0_r, /*sigma_t=*/0.5_r));
    sim_empty.run(150);

    Real const e_empty = total_energy_2d(sim_empty.current(),
                                         sim_empty.previous(), c0, dt);

    // --- molecule-loaded medium ---
    Field<Real, 2> rho(grid, 0.0_r);
    splat_density(rho, scene2d);

    auto medium_mol = Medium<2>::uniform(grid, c0);
    MediumWeights w;
    w.beta_rho = 0.3_r;
    build_medium_from_fields<2>(medium_mol, &rho, nullptr, nullptr, c0, w);
    apply_pml(medium_mol, PmlSpec{/*cells=*/15, /*alpha_max_factor=*/2.0_r, /*polynomial_order=*/3});

    FdtdCpuOmp<2> sim_mol(grid, std::move(medium_mol), dt);
    sim_mol.set_boundary(std::make_shared<Dirichlet2D>());
    sim_mol.add_source(std::make_shared<GaussianPulse<2>>(
        src_loc, 5.0_r, 3.0_r, 2.0_r, 0.5_r));
    sim_mol.run(150);

    Real const e_mol = total_energy_2d(sim_mol.current(),
                                       sim_mol.previous(),
                                       sim_mol.wave_speed(), dt);

    // Scattering causes the field configurations to differ — measure via
    // wave_difference normalized by the empty-baseline norm.
    Real const d = wave_difference<2>(sim_empty.current(), sim_mol.current());
    Real norm_empty = 0.0_r;
    for (Index k = 0; k < sim_empty.current().size(); ++k) {
        Real const v = sim_empty.current().data()[k];
        norm_empty += v * v;
    }
    norm_empty = std::sqrt(norm_empty);
    Real const rel_diff = d / std::max(norm_empty, 1e-6_r);

    CAPTURE(e_empty);
    CAPTURE(e_mol);
    CAPTURE(d);
    CAPTURE(norm_empty);
    CAPTURE(rel_diff);
    // Non-trivial scattering: the molecule-loaded field differs from the
    // empty-baseline field by at least 10% of the baseline's L2 norm.
    CHECK(rel_diff > 0.1_r);
    // Sanity: still finite, not blown up.
    CHECK(std::isfinite(e_mol));
    CHECK(e_mol > 0.0_r);
}

} // namespace
