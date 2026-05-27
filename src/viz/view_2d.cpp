#include "viz/view_2d.hpp"

#include "fdtd/stability.hpp"
#include "medium/boundary.hpp"
#include "medium/build_medium.hpp"
#include "medium/pml.hpp"
#include "score/energy.hpp"

#include <imgui.h>
#include <raylib.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace wavelab {

View2D::View2D(Config cfg) : cfg_(cfg) {
    rebuild_sim_();
}

View2D::~View2D() {
    if (gpu_ready_) {
        UnloadTexture(tex_);
        UnloadImage(image_);
    }
}

void View2D::rebuild_sim_() {
    auto grid = Grid<2>::uniform(IVec<2>{cfg_.nx, cfg_.ny}, cfg_.dx);
    auto medium = Medium<2>::uniform(grid, cfg_.c0);

    if (cfg_.insert_slab) {
        Real const c_inside = cfg_.c0 / cfg_.slab_n;
        apply_index_slab(medium, cfg_.slab_x_lo, cfg_.slab_x_hi, c_inside);
    }

    PmlSpec pml_spec;
    pml_spec.cells = cfg_.pml_cells;
    pml_spec.alpha_max_factor = 2.0_r;
    apply_pml(medium, pml_spec);

    Real const dt = cfg_.cfl_safety * cfl_dt_max<2>(grid, cfg_.c0);
    sim_ = std::make_unique<FdtdCpuOmp<2>>(grid, std::move(medium), dt);
    sim_->set_boundary(std::make_shared<Dirichlet2D>());

    // Point source near the left edge, on-axis.
    Index const src_i = cfg_.pml_cells + 4;
    Index const src_j = cfg_.ny / 2;
    source_ = std::make_shared<HarmonicSource<2>>(
        IVec<2>{src_i, src_j}, cfg_.source_amp, cfg_.source_freq);
    sim_->add_source(source_);

    // Allocate / resize the GPU-backed image.
    if (gpu_ready_) {
        UnloadTexture(tex_);
        UnloadImage(image_);
        gpu_ready_ = false;
    }
    rgba_.assign(static_cast<std::size_t>(cfg_.nx * cfg_.ny * 4), 0);
    image_.data    = rgba_.data();
    image_.width   = static_cast<int>(cfg_.nx);
    image_.height  = static_cast<int>(cfg_.ny);
    image_.mipmaps = 1;
    image_.format  = PIXELFORMAT_UNCOMPRESSED_R8G8B8A8;
    tex_           = LoadTextureFromImage(image_);
    SetTextureFilter(tex_, TEXTURE_FILTER_BILINEAR);
    gpu_ready_     = true;
}

void View2D::update() {
    if (playing_) sim_->run(static_cast<Index>(steps_per_frame_));
    recolor_();
}

void View2D::recolor_() {
    Index const nx = cfg_.nx;
    Index const ny = cfg_.ny;
    unsigned char* px = rgba_.data();

    if (mode_ == OverlayMode::Field) {
        auto const& u = sim_->current();
        for (Index i = 0; i < nx; ++i) {
            for (Index j = 0; j < ny; ++j) {
                unsigned char r, g, b;
                colormap_diverging_(u(i, j), vmax_field_, r, g, b);
                Index const k = (i * ny + j) * 4;
                px[k + 0] = r; px[k + 1] = g; px[k + 2] = b; px[k + 3] = 255;
            }
        }
    } else if (mode_ == OverlayMode::RefractiveIndex) {
        auto const& cf = sim_->medium().c;
        Real const c0 = cfg_.c0;
        Real n_max = 1.0_r;
        for (Index k = 0; k < cf.size(); ++k) {
            Real const n = c0 / cf.data()[k];
            if (n > n_max) n_max = n;
        }
        for (Index i = 0; i < nx; ++i) {
            for (Index j = 0; j < ny; ++j) {
                Real const n = c0 / cf(i, j);
                unsigned char r, g, b;
                colormap_sequential_(n - 1.0_r, std::max(0.001_r, n_max - 1.0_r),
                                      r, g, b);
                Index const k = (i * ny + j) * 4;
                px[k + 0] = r; px[k + 1] = g; px[k + 2] = b; px[k + 3] = 255;
            }
        }
    } else { // EnergyDensity
        Real const cmax = sim_->wave_speed();
        // Cheap moving-max auto-scale
        Real vmax = vmax_energy_ * 0.95_r;
        for (Index i = 1; i < nx - 1; ++i) {
            for (Index j = 1; j < ny - 1; ++j) {
                Real const e = energy_density_2d(sim_->current(),
                                                 sim_->previous(),
                                                 cmax, sim_->dt(), i, j);
                if (e > vmax) vmax = e;
            }
        }
        vmax_energy_ = std::max(vmax, 1e-6_r);

        for (Index i = 0; i < nx; ++i) {
            for (Index j = 0; j < ny; ++j) {
                Real e = 0.0_r;
                if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1) {
                    e = energy_density_2d(sim_->current(), sim_->previous(),
                                          cmax, sim_->dt(), i, j);
                }
                unsigned char r, g, b;
                colormap_sequential_(e, vmax_energy_, r, g, b);
                Index const k = (i * ny + j) * 4;
                px[k + 0] = r; px[k + 1] = g; px[k + 2] = b; px[k + 3] = 255;
            }
        }
    }

    UpdateTexture(tex_, rgba_.data());
}

void View2D::draw_controls() {
    ImGui::Begin("View2D — controls");

    ImGui::Text("Phase 2: 2D wave engine");
    ImGui::Separator();

    if (ImGui::Button(playing_ ? "Pause" : "Play")) playing_ = !playing_;
    ImGui::SameLine();
    if (ImGui::Button("Reset")) sim_->reset();
    ImGui::SameLine();
    if (ImGui::Button("Rebuild")) rebuild_sim_();
    ImGui::SliderInt("Steps / frame", &steps_per_frame_, 1, 32);

    ImGui::Separator();
    ImGui::Text("Time:  %.4f", static_cast<double>(sim_->time()));
    ImGui::Text("Step:  %lld", static_cast<long long>(sim_->step_count()));
    ImGui::Text("dt:    %.5f", static_cast<double>(sim_->dt()));
    ImGui::Text("Grid:  %lld x %lld", static_cast<long long>(cfg_.nx),
                                       static_cast<long long>(cfg_.ny));

    Real const e = total_energy_2d(sim_->current(), sim_->previous(),
                                   sim_->wave_speed(), sim_->dt());
    ImGui::Text("Total energy: %.4f", static_cast<double>(e));

    ImGui::Separator();
    ImGui::Text("Source");
    float freq = static_cast<float>(source_->frequency());
    float amp  = static_cast<float>(source_->amplitude());
    if (ImGui::SliderFloat("Frequency", &freq, 0.1f, 8.0f, "%.2f")) {
        source_->set_frequency(static_cast<Real>(freq));
    }
    if (ImGui::SliderFloat("Amplitude", &amp, 0.0f, 2.0f, "%.2f")) {
        source_->set_amplitude(static_cast<Real>(amp));
    }

    ImGui::Separator();
    ImGui::Text("Medium");
    bool slab = cfg_.insert_slab;
    if (ImGui::Checkbox("Insert refractive slab", &slab)) {
        cfg_.insert_slab = slab;
        rebuild_sim_();
    }
    if (cfg_.insert_slab) {
        float n_val = static_cast<float>(cfg_.slab_n);
        if (ImGui::SliderFloat("Slab n", &n_val, 1.1f, 4.0f, "%.2f")) {
            cfg_.slab_n = static_cast<Real>(n_val);
            rebuild_sim_();
        }
    }

    ImGui::Separator();
    ImGui::Text("Overlay");
    int m = static_cast<int>(mode_);
    bool changed = false;
    changed |= ImGui::RadioButton("Field u",       &m, 0); ImGui::SameLine();
    changed |= ImGui::RadioButton("Index n",       &m, 1); ImGui::SameLine();
    changed |= ImGui::RadioButton("Energy",        &m, 2);
    if (changed) mode_ = static_cast<OverlayMode>(m);

    if (mode_ == OverlayMode::Field) {
        float v = static_cast<float>(vmax_field_);
        if (ImGui::SliderFloat("Field scale", &v, 0.01f, 5.0f, "%.2f")) {
            vmax_field_ = static_cast<Real>(v);
        }
    }

    ImGui::End();
}

void View2D::draw_field(Rectangle area) {
    if (!gpu_ready_) return;

    // Preserve aspect ratio of the simulation grid.
    float const grid_aspect = static_cast<float>(cfg_.nx)
                            / static_cast<float>(cfg_.ny);
    float w = area.width;
    float h = area.height;
    if (w / h > grid_aspect) {
        w = h * grid_aspect;
    } else {
        h = w / grid_aspect;
    }
    Rectangle dst{area.x + 0.5f * (area.width - w),
                  area.y + 0.5f * (area.height - h),
                  w, h};

    // Source rectangle = full texture; texture is (nx, ny) which is stored
    // as (i along x, j along y). Raylib texture coords assume (width=nx,
    // height=ny), so we draw rotated/swapped as desired. For this engine
    // i is the first index (x), j is second (y). We render so the i-axis
    // goes RIGHTWARD on screen — that means texture rows correspond to i
    // (x), columns to j (y). Standard image layout has rows = y, cols = x,
    // so we draw with width=nx and height=ny but with no rotation — the
    // resulting image has x running DOWN. To make x → right and y → up,
    // we'd want a transpose; for Phase 2 the canonical view just shows
    // the field in matrix order, which is fine for inspection.
    Rectangle src{0, 0, static_cast<float>(cfg_.nx),
                       static_cast<float>(cfg_.ny)};
    DrawTexturePro(tex_, src, dst, Vector2{0, 0}, 0.0f, WHITE);
    DrawRectangleLinesEx(dst, 1.0f, Color{60, 64, 78, 255});
}

// --- colormaps -----------------------------------------------------------

void View2D::colormap_diverging_(Real v, Real vmax,
                                 unsigned char& r,
                                 unsigned char& g,
                                 unsigned char& b) {
    if (vmax <= Real{0}) { r = g = b = 0; return; }
    Real t = std::clamp(v / vmax, -1.0_r, 1.0_r);
    if (t >= Real{0}) {
        // white → red
        Real const w = Real{1} - t;
        r = static_cast<unsigned char>(255);
        g = static_cast<unsigned char>(255 * w);
        b = static_cast<unsigned char>(255 * w);
    } else {
        // white → blue
        Real const w = Real{1} + t;
        r = static_cast<unsigned char>(255 * w);
        g = static_cast<unsigned char>(255 * w);
        b = static_cast<unsigned char>(255);
    }
}

void View2D::colormap_sequential_(Real v, Real vmax,
                                  unsigned char& r,
                                  unsigned char& g,
                                  unsigned char& b) {
    if (vmax <= Real{0}) { r = g = b = 0; return; }
    Real const t = std::clamp(v / vmax, 0.0_r, 1.0_r);
    // Viridis-ish: dark purple → cyan → yellow
    if (t < 0.5_r) {
        Real const s = t * 2.0_r;
        r = static_cast<unsigned char>(60  + s * (40  - 60));
        g = static_cast<unsigned char>(20  + s * (180 - 20));
        b = static_cast<unsigned char>(120 + s * (180 - 120));
    } else {
        Real const s = (t - 0.5_r) * 2.0_r;
        r = static_cast<unsigned char>(40  + s * (250 - 40));
        g = static_cast<unsigned char>(180 + s * (230 - 180));
        b = static_cast<unsigned char>(180 + s * (40  - 180));
    }
}

} // namespace wavelab
