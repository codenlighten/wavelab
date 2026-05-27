#pragma once
//
// View2D — interactive ImGui + raylib view for a 2D FDTD simulation.
//
// Owns the FdtdCpuOmp<2> stepper, a HarmonicSource, and the current
// Medium. Renders the field as a diverging-colormap heatmap into a
// raylib Texture. Supports overlay toggles for the refractive-index
// field and the energy-density field.
//
// Like View1D, the stepper runs inline (N steps per frame). Engine⊥viz
// threading lands when 3D makes per-frame stepping too slow (Phase 5).
//

#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "medium/medium.hpp"
#include "source/harmonic.hpp"

#include <raylib.h>

#include <memory>
#include <vector>

namespace wavelab {

class View2D {
public:
    struct Config {
        Index nx          = 256;
        Index ny          = 256;
        Real  dx          = 0.1_r;
        Real  c0          = 1.0_r;
        Real  cfl_safety  = 0.4_r;
        Real  source_freq = 2.0_r;
        Real  source_amp  = 1.0_r;
        Index pml_cells   = 24;
        bool  insert_slab = false;
        Real  slab_x_lo   = 14.0_r;
        Real  slab_x_hi   = 18.0_r;
        Real  slab_n      = 1.6_r;
    };

    enum class OverlayMode {
        Field,           // u(x,y,t)
        RefractiveIndex, // n = c0/c per cell
        EnergyDensity,   // ½u_t² + ½c²|∇u|²
    };

    explicit View2D(Config cfg = {});
    ~View2D();

    View2D(View2D const&)            = delete;
    View2D& operator=(View2D const&) = delete;

    void update();
    void draw_controls();
    void draw_field(Rectangle area);

    FdtdCpuOmp<2>&       sim()       noexcept { return *sim_; }
    FdtdCpuOmp<2> const& sim() const noexcept { return *sim_; }

private:
    void rebuild_sim_();
    void recolor_();

    // Paint a single value into RGBA8. Used for all overlay modes; the
    // value range and divergence behavior change per mode.
    static void colormap_diverging_(Real v, Real vmax,
                                    unsigned char& r, unsigned char& g,
                                    unsigned char& b);
    static void colormap_sequential_(Real v, Real vmax,
                                     unsigned char& r, unsigned char& g,
                                     unsigned char& b);

    Config                              cfg_;
    std::unique_ptr<FdtdCpuOmp<2>>      sim_;
    std::shared_ptr<HarmonicSource<2>>  source_;

    bool  playing_         = true;
    int   steps_per_frame_ = 4;
    Real  vmax_field_      = 1.0_r;
    Real  vmax_energy_     = 0.0_r;     // auto-scale
    OverlayMode mode_      = OverlayMode::Field;

    // Texture-backed RGBA buffer for raylib upload.
    Image       image_{};
    Texture2D   tex_{};
    bool        gpu_ready_ = false;
    std::vector<unsigned char> rgba_;
};

} // namespace wavelab
