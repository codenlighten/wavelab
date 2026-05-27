#pragma once
//
// View1D — interactive ImGui + raylib view for a 1D FDTD simulation.
//
// Owns the FdtdCpuOmp<1> stepper and the active HarmonicSource. Provides:
//   * draw_controls() — ImGui panel (play/pause, reset, freq/amp/damping
//                      sliders, boundary selector, step count, energy).
//   * draw_field(area) — raylib line plot of u(x, t) inside the given
//                        screen-space rectangle.
//
// Phase 1 runs the stepper inline (N steps per frame). Engine⊥viewer
// threading lands when 2D simulations make per-frame stepping too slow.
//

#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "medium/boundary.hpp"
#include "source/harmonic.hpp"

#include <raylib.h>

#include <memory>

namespace wavelab {

class View1D {
public:
    struct Config {
        Index nx          = 800;
        Real  dx          = 0.01_r;
        Real  c           = 1.0_r;
        Real  cfl_safety  = 0.5_r;
        Real  source_freq = 2.0_r;
        Real  source_amp  = 1.0_r;
        Real  damping     = 0.0_r;
        Index source_cell = 0;
    };

    enum class BoundaryKind { Dirichlet, Neumann, Mur };

    explicit View1D(Config cfg = {});

    // Advance the simulation by `steps_per_frame_` if playing.
    void update();

    // ImGui controls panel.
    void draw_controls();

    // Field plot inside the given screen-space rectangle.
    void draw_field(Rectangle area) const;

    FdtdCpuOmp<1>&       sim()       noexcept { return *sim_; }
    FdtdCpuOmp<1> const& sim() const noexcept { return *sim_; }

private:
    void rebuild_sim_();
    void install_boundary_();

    Config                              cfg_;
    std::unique_ptr<FdtdCpuOmp<1>>      sim_;
    std::shared_ptr<HarmonicSource<1>>  source_;
    BoundaryKind                        boundary_kind_ = BoundaryKind::Mur;

    bool  playing_         = true;
    int   steps_per_frame_ = 4;
    Real  ymin_            = -1.5_r;
    Real  ymax_            =  1.5_r;
};

} // namespace wavelab
