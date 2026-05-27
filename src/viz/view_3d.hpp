#pragma once
//
// View3D — orthographic three-plane slicer for a 3D FDTD simulation.
//
// Renders three orthogonal 2D slices (xy, xz, yz) of the chosen overlay
// (field u, refractive index n, or energy density). Each slice is drawn
// as a raylib texture using the same diverging/sequential colormaps as
// View2D. Slice indices are ImGui-controlled sliders.
//
// Volumetric raymarcher is the Phase 5 stretch goal — deferred.
//

#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "medium/medium.hpp"
#include "source/harmonic.hpp"

#include <raylib.h>

#include <memory>
#include <vector>

namespace wavelab {

class View3D {
public:
    struct Config {
        Index nx          = 64;
        Index ny          = 64;
        Index nz          = 64;
        Real  dx          = 0.5_r;
        Real  c0          = 1.0_r;
        Real  cfl_safety  = 0.4_r;
        Real  source_freq = 2.0_r;
        Real  source_amp  = 1.0_r;
        Index pml_cells   = 12;
    };

    enum class OverlayMode {
        Field,
        RefractiveIndex,
        EnergyDensity,
    };

    explicit View3D(Config cfg = {});
    ~View3D();

    View3D(View3D const&)            = delete;
    View3D& operator=(View3D const&) = delete;

    void update();
    void draw_controls();
    void draw_field(Rectangle area);

    FdtdCpuOmp<3>&       sim()       noexcept { return *sim_; }
    FdtdCpuOmp<3> const& sim() const noexcept { return *sim_; }

private:
    void rebuild_sim_();
    void recolor_();
    void resize_textures_();

    static void colormap_diverging_(Real v, Real vmax,
                                    unsigned char& r, unsigned char& g, unsigned char& b);
    static void colormap_sequential_(Real v, Real vmax,
                                     unsigned char& r, unsigned char& g, unsigned char& b);

    Config                              cfg_;
    std::unique_ptr<FdtdCpuOmp<3>>      sim_;
    std::shared_ptr<HarmonicSource<3>>  source_;

    bool  playing_         = true;
    int   steps_per_frame_ = 2;
    Real  vmax_            = 1.0_r;
    OverlayMode mode_      = OverlayMode::Field;

    int slice_x_ = 0;
    int slice_y_ = 0;
    int slice_z_ = 0;

    // Three slice textures (xy, xz, yz)
    Image     img_xy_{}, img_xz_{}, img_yz_{};
    Texture2D tex_xy_{}, tex_xz_{}, tex_yz_{};
    bool      gpu_ready_ = false;
    std::vector<unsigned char> rgba_xy_, rgba_xz_, rgba_yz_;
};

} // namespace wavelab
