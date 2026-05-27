// Python bindings for wavelab.
//
// Exposes a small surface (no full FDTD object hierarchy) — most users
// just want to run a pipeline and inspect the resulting fingerprint
// scalars + spectral vector. Power users can fall back to the C++ API.
//
//   import wavelab
//   fp2 = wavelab.run_pdb("1ubq.pdb", nx=80, ny=80, steps=200)         # 2D
//   fp3 = wavelab.run_pdb("1ubq.pdb", dim=3, nx=64, ny=64, nz=64)      # 3D
//   print(fp2["scalars"]["total_energy"], fp2["spectral"][:4])
//
// `run_synthetic`, `build_prototype`, and `correlate` mirror wavecli's
// CLI but return Python dicts directly. `run_synthetic` is 2D only.

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
#include "molecule/scene.hpp"
#include "molecule/slice.hpp"
#include "molecule/synthetic.hpp"
#include "score/calibration.hpp"
#include "score/energy.hpp"
#include "score/entropy.hpp"
#include "score/spectral.hpp"
#include "source/harmonic.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace py = pybind11;
using namespace wavelab;
using namespace wavelab::literals;

namespace {

struct PipelineParams {
    int   dim    = 2;          // 2 or 3
    Index nx     = 100;
    Index ny     = 100;
    Index nz     = 100;
    Real  dx     = 0.5_r;
    Real  freq   = 2.0_r;
    Index steps  = 200;
    Index pml_cells = 15;
    Real  c0     = 1.0_r;
    Real  beta_rho = 0.3_r;
    Real  slice_thickness = 4.0_r;
};

// Pipeline body, templated on D. Mirrors wavecli's run_pipeline<D>.
template <int D>
Fingerprint run_one_(MolecularScene<D> const& scene,
                     PipelineParams const& p,
                     std::string const& scene_name) {
    if (scene.atoms.empty()) {
        throw std::runtime_error("run_pipeline: scene has no atoms");
    }
    Grid<D> grid = [&]() {
        if constexpr (D == 2) {
            return Grid<2>{IVec<2>{p.nx, p.ny}, Vec<2>{p.dx, p.dx},
                           scene.box_min};
        } else {
            return Grid<3>{IVec<3>{p.nx, p.ny, p.nz},
                           Vec<3>{p.dx, p.dx, p.dx}, scene.box_min};
        }
    }();

    Field<Real, D> rho(grid, 0.0_r);
    splat_density(rho, scene);

    auto medium = Medium<D>::uniform(grid, p.c0);
    MediumWeights w;
    w.beta_rho = p.beta_rho;
    build_medium_from_fields<D>(medium, &rho, nullptr, nullptr, p.c0, w);
    apply_pml(medium, PmlSpec{p.pml_cells, /*alpha_max_factor=*/2.0_r,
                              /*polynomial_order=*/3});

    Real const dt = 0.4_r * cfl_dt_max<D>(grid, p.c0);
    FdtdCpuOmp<D> sim(grid, std::move(medium), dt);
    sim.set_boundary(make_dirichlet<D>());

    IVec<D> src = [&]() {
        if constexpr (D == 2)
            return IVec<2>{p.pml_cells + 5, p.ny / 2};
        else
            return IVec<3>{p.pml_cells + 5, p.ny / 2, p.nz / 2};
    }();
    sim.add_source(std::make_shared<HarmonicSource<D>>(src, /*amp=*/1.0_r, p.freq));

    IVec<D> probe = [&]() {
        if constexpr (D == 2) return IVec<2>{p.nx / 4, p.ny / 2};
        else                  return IVec<3>{p.nx / 4, p.ny / 2, p.nz / 2};
    }();
    std::vector<Real> series;
    series.reserve(static_cast<std::size_t>(p.steps));
    for (Index s = 0; s < p.steps; ++s) {
        sim.step();
        probe_record<D>(sim.current(), probe, series);
    }

    Fingerprint fp;
    fp.scene_name = scene_name;
    fp.scalars["total_energy"] = total_energy<D>(sim.current(), sim.previous(),
                                                  sim.wave_speed(), dt);
    Field<Real, D> e_field(grid, 0.0_r);
    energy_density_field<D>(e_field, sim.current(), sim.previous(),
                            sim.wave_speed(), dt);
    fp.scalars["entropy"]     = energy_entropy<D>(e_field);
    fp.scalars["focus_score"] = focus_score<D>(e_field);

    auto power = compute_power_spectrum(series);
    auto freqs = compute_frequency_axis(series.size(), dt);
    fp.spectral = spectral_fingerprint_logbins(power, freqs, /*nbands=*/16,
                                               /*f_min=*/0.1_r, /*f_max=*/10.0_r);
    fp.spectral_freqs = std::vector<Real>(fp.spectral.size(), 0.0_r);
    if (!fp.spectral.empty()) {
        Real const log_lo = std::log(0.1_r);
        Real const log_hi = std::log(10.0_r);
        Real const dlog   = (log_hi - log_lo) / static_cast<Real>(fp.spectral.size());
        for (std::size_t i = 0; i < fp.spectral.size(); ++i) {
            fp.spectral_freqs[i] = std::exp(log_lo + (static_cast<Real>(i) + 0.5_r) * dlog);
        }
    }
    fp.meta["dim"]   = std::to_string(D);
    fp.meta["atoms"] = std::to_string(scene.atoms.size());
    return fp;
}

py::dict fingerprint_to_dict_(Fingerprint const& fp) {
    py::dict d;
    d["scene_name"]     = fp.scene_name;
    d["scalars"]        = fp.scalars;
    d["spectral"]       = fp.spectral;
    d["spectral_freqs"] = fp.spectral_freqs;
    d["meta"]           = fp.meta;
    return d;
}

Fingerprint dict_to_fingerprint_(py::dict const& d) {
    Fingerprint fp;
    if (d.contains("scene_name"))    fp.scene_name     = d["scene_name"].cast<std::string>();
    if (d.contains("scalars"))       fp.scalars        = d["scalars"].cast<std::unordered_map<std::string, Real>>();
    if (d.contains("spectral"))      fp.spectral       = d["spectral"].cast<std::vector<Real>>();
    if (d.contains("spectral_freqs")) fp.spectral_freqs = d["spectral_freqs"].cast<std::vector<Real>>();
    if (d.contains("meta"))          fp.meta           = d["meta"].cast<std::unordered_map<std::string, std::string>>();
    return fp;
}

PipelineParams params_from_kwargs(py::kwargs const& k) {
    PipelineParams p;
    if (k.contains("dim"))           p.dim   = k["dim"].cast<int>();
    if (k.contains("nx"))            p.nx    = k["nx"].cast<Index>();
    if (k.contains("ny"))            p.ny    = k["ny"].cast<Index>();
    if (k.contains("nz"))            p.nz    = k["nz"].cast<Index>();
    if (k.contains("dx"))            p.dx    = k["dx"].cast<Real>();
    if (k.contains("freq"))          p.freq  = k["freq"].cast<Real>();
    if (k.contains("steps"))         p.steps = k["steps"].cast<Index>();
    if (k.contains("pml_cells"))     p.pml_cells = k["pml_cells"].cast<Index>();
    if (k.contains("c0"))            p.c0    = k["c0"].cast<Real>();
    if (k.contains("beta_rho"))      p.beta_rho = k["beta_rho"].cast<Real>();
    if (k.contains("slice_thickness")) p.slice_thickness = k["slice_thickness"].cast<Real>();
    if (p.dim != 2 && p.dim != 3) {
        throw std::invalid_argument("dim must be 2 or 3");
    }
    return p;
}

py::dict run_pdb_(std::string const& path, py::kwargs k) {
    auto const p = params_from_kwargs(k);
    auto r3 = parse_pdb_file(path);
    if (r3.scene.atoms.empty()) {
        throw std::runtime_error("run_pdb: no atoms parsed from " + path);
    }
    if (p.dim == 2) {
        Real const z_mid = 0.5_r * (r3.scene.box_min[2] + r3.scene.box_max[2]);
        auto scene = slice_scene_xy_centered(r3.scene, z_mid, p.slice_thickness);
        return fingerprint_to_dict_(run_one_<2>(scene, p, path));
    }
    return fingerprint_to_dict_(run_one_<3>(r3.scene, p, path));
}

py::dict run_synthetic_(std::string const& kind, py::kwargs k) {
    auto const p = params_from_kwargs(k);
    if (p.dim != 2) {
        throw std::invalid_argument("run_synthetic: D=3 synthetics not implemented");
    }
    Real const box_x = static_cast<Real>(p.nx) * p.dx;
    Real const box_y = static_cast<Real>(p.ny) * p.dx;
    Vec<2> bmin{0.0_r, 0.0_r};
    Vec<2> bmax{box_x, box_y};
    Vec<2> center{0.5_r * box_x, 0.5_r * box_y};

    MolecularScene<2> scene;
    if (kind == "single")        scene = single_atom_2d(center, 2.0_r, bmin, bmax);
    else if (kind == "dumbbell") scene = dumbbell_2d(
        Vec<2>{center[0] - 3.0_r, center[1]},
        Vec<2>{center[0] + 3.0_r, center[1]},
        2.0_r, bmin, bmax);
    else if (kind == "pocket")   scene = pocket_arc_2d(
        center, /*R=*/5.0_r, 2.0_r, /*count=*/24,
        /*opening_dir=*/0.0_r, /*opening_rad=*/1.0_r, bmin, bmax);
    else if (kind == "slab")     scene = slab_atoms_2d(
        center[0] - 1.0_r, center[0] + 1.0_r, 1.5_r, 2.0_r, bmin, bmax);
    else throw std::runtime_error("unknown synthetic kind: " + kind);

    return fingerprint_to_dict_(run_one_<2>(scene, p, kind));
}

py::dict build_prototype_(py::list const& binders_list) {
    std::vector<Fingerprint> binders;
    binders.reserve(binders_list.size());
    for (auto h : binders_list) {
        binders.push_back(dict_to_fingerprint_(py::cast<py::dict>(h)));
    }
    auto proto = build_prototype(binders);
    py::dict d;
    d["scalar_mean"]    = proto.scalar_mean;
    d["spectral_mean"]  = proto.spectral_mean;
    d["spectral_freqs"] = proto.spectral_freqs;
    d["sample_count"]   = proto.sample_count;
    return d;
}

Real correlate_(py::dict const& candidate, py::dict const& proto_d) {
    Fingerprint c = dict_to_fingerprint_(candidate);
    BinderPrototype p;
    if (proto_d.contains("scalar_mean"))    p.scalar_mean    = proto_d["scalar_mean"].cast<std::unordered_map<std::string, Real>>();
    if (proto_d.contains("spectral_mean"))  p.spectral_mean  = proto_d["spectral_mean"].cast<std::vector<Real>>();
    if (proto_d.contains("spectral_freqs")) p.spectral_freqs = proto_d["spectral_freqs"].cast<std::vector<Real>>();
    if (proto_d.contains("sample_count"))   p.sample_count   = proto_d["sample_count"].cast<Index>();
    return correlate_with_prototype(c, p);
}

} // namespace

PYBIND11_MODULE(wavelab, m) {
    m.doc() = "wavelab — wave-based molecular geometry engine";

    m.def("run_pdb",       &run_pdb_,        py::arg("path"),
          "Load a PDB, run a 2D midplane-slice or full 3D simulation "
          "(via kwarg dim=2|3), return a fingerprint dict.");
    m.def("run_synthetic", &run_synthetic_,  py::arg("kind"),
          "Run on a synthetic 2D scene (single|dumbbell|pocket|slab); "
          "returns a fingerprint dict.");

    m.def("build_prototype", &build_prototype_, py::arg("fingerprints"),
          "Aggregate K fingerprint dicts into a prototype dict (§34).");
    m.def("correlate",       &correlate_,       py::arg("candidate"), py::arg("prototype"),
          "Cosine similarity between a candidate fingerprint and a prototype.");
}
