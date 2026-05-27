#include "viz/view_3d.hpp"

#include "core/memory_budget.hpp"
#include "fdtd/stability.hpp"
#include "medium/boundary.hpp"
#include "medium/medium.hpp"
#include "medium/pml.hpp"
#include "score/energy.hpp"

#include <imgui.h>
#include <raylib.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace wavelab {

View3D::View3D(Config cfg) : cfg_(cfg) {
    rebuild_sim_();
}

View3D::~View3D() {
    if (gpu_ready_) {
        UnloadTexture(tex_xy_); UnloadImage(img_xy_);
        UnloadTexture(tex_xz_); UnloadImage(img_xz_);
        UnloadTexture(tex_yz_); UnloadImage(img_yz_);
    }
}

void View3D::rebuild_sim_() {
    auto grid = Grid<3>::uniform(IVec<3>{cfg_.nx, cfg_.ny, cfg_.nz}, cfg_.dx);
    // Guardrail: prevent accidental gigabyte allocations.
    enforce_memory_budget<3>(grid, /*cap=*/2ull * 1024 * 1024 * 1024);

    auto medium = Medium<3>::uniform(grid, cfg_.c0);
    PmlSpec pml{cfg_.pml_cells, /*alpha_max_factor=*/2.0_r, /*polynomial_order=*/3};
    apply_pml(medium, pml);

    Real const dt = cfg_.cfl_safety * cfl_dt_max<3>(grid, cfg_.c0);
    sim_ = std::make_unique<FdtdCpuOmp<3>>(grid, std::move(medium), dt);
    sim_->set_boundary(std::make_shared<Dirichlet3D>());

    source_ = std::make_shared<HarmonicSource<3>>(
        IVec<3>{cfg_.nx / 2, cfg_.ny / 2, cfg_.nz / 2},
        cfg_.source_amp, cfg_.source_freq);
    sim_->add_source(source_);

    slice_x_ = static_cast<int>(cfg_.nx / 2);
    slice_y_ = static_cast<int>(cfg_.ny / 2);
    slice_z_ = static_cast<int>(cfg_.nz / 2);

    resize_textures_();
}

void View3D::resize_textures_() {
    if (gpu_ready_) {
        UnloadTexture(tex_xy_); UnloadImage(img_xy_);
        UnloadTexture(tex_xz_); UnloadImage(img_xz_);
        UnloadTexture(tex_yz_); UnloadImage(img_yz_);
        gpu_ready_ = false;
    }
    auto make_img = [&](Image& im, Texture2D& tx, std::vector<unsigned char>& buf,
                        int w, int h) {
        buf.assign(static_cast<std::size_t>(w * h * 4), 0);
        im.data    = buf.data();
        im.width   = w;
        im.height  = h;
        im.mipmaps = 1;
        im.format  = PIXELFORMAT_UNCOMPRESSED_R8G8B8A8;
        tx = LoadTextureFromImage(im);
        SetTextureFilter(tx, TEXTURE_FILTER_BILINEAR);
    };
    make_img(img_xy_, tex_xy_, rgba_xy_, static_cast<int>(cfg_.nx), static_cast<int>(cfg_.ny));
    make_img(img_xz_, tex_xz_, rgba_xz_, static_cast<int>(cfg_.nx), static_cast<int>(cfg_.nz));
    make_img(img_yz_, tex_yz_, rgba_yz_, static_cast<int>(cfg_.ny), static_cast<int>(cfg_.nz));
    gpu_ready_ = true;
}

void View3D::update() {
    if (playing_) sim_->run(static_cast<Index>(steps_per_frame_));
    recolor_();
}

void View3D::recolor_() {
    if (!gpu_ready_) return;
    auto const& u = sim_->current();
    Index const nx = cfg_.nx, ny = cfg_.ny, nz = cfg_.nz;

    auto value_at = [&](Index i, Index j, Index k) -> Real {
        if (mode_ == OverlayMode::Field) return u(i, j, k);
        if (mode_ == OverlayMode::RefractiveIndex) {
            return cfg_.c0 / sim_->medium().c(i, j, k);
        }
        // EnergyDensity
        Real const cmax = sim_->wave_speed();
        return energy_density_3d(sim_->current(), sim_->previous(),
                                 cmax, sim_->dt(), i, j, k);
    };
    auto paint = [&](unsigned char* px, Real v, Real vmax) {
        unsigned char r, g, b;
        if (mode_ == OverlayMode::Field) {
            colormap_diverging_(v, vmax, r, g, b);
        } else {
            colormap_sequential_(v - (mode_ == OverlayMode::RefractiveIndex ? 1.0_r : 0.0_r),
                                 vmax, r, g, b);
        }
        px[0] = r; px[1] = g; px[2] = b; px[3] = 255;
    };

    Real vmax = (mode_ == OverlayMode::Field) ? vmax_ : 1.0_r;
    if (mode_ != OverlayMode::Field) {
        // Quick max scan over the 3 slices for auto-scale.
        Real m = 0.0_r;
        for (Index i = 0; i < nx; ++i)
            for (Index j = 0; j < ny; ++j) m = std::max(m, std::abs(value_at(i, j, slice_z_)));
        for (Index i = 0; i < nx; ++i)
            for (Index k = 0; k < nz; ++k) m = std::max(m, std::abs(value_at(i, slice_y_, k)));
        for (Index j = 0; j < ny; ++j)
            for (Index k = 0; k < nz; ++k) m = std::max(m, std::abs(value_at(slice_x_, j, k)));
        vmax = std::max(m, 1e-6_r);
    }

    for (Index i = 0; i < nx; ++i)
        for (Index j = 0; j < ny; ++j) {
            std::size_t off = static_cast<std::size_t>(i * ny + j) * 4u;
            paint(rgba_xy_.data() + off, value_at(i, j, slice_z_), vmax);
        }
    for (Index i = 0; i < nx; ++i)
        for (Index k = 0; k < nz; ++k) {
            std::size_t off = static_cast<std::size_t>(i * nz + k) * 4u;
            paint(rgba_xz_.data() + off, value_at(i, slice_y_, k), vmax);
        }
    for (Index j = 0; j < ny; ++j)
        for (Index k = 0; k < nz; ++k) {
            std::size_t off = static_cast<std::size_t>(j * nz + k) * 4u;
            paint(rgba_yz_.data() + off, value_at(slice_x_, j, k), vmax);
        }

    UpdateTexture(tex_xy_, rgba_xy_.data());
    UpdateTexture(tex_xz_, rgba_xz_.data());
    UpdateTexture(tex_yz_, rgba_yz_.data());
}

void View3D::draw_controls() {
    ImGui::Begin("View3D — controls");
    ImGui::Text("Phase 5: 3D wave engine (slicer)");
    ImGui::Separator();

    if (ImGui::Button(playing_ ? "Pause" : "Play")) playing_ = !playing_;
    ImGui::SameLine();
    if (ImGui::Button("Reset")) sim_->reset();
    ImGui::SameLine();
    if (ImGui::Button("Rebuild")) rebuild_sim_();
    ImGui::SliderInt("Steps / frame", &steps_per_frame_, 1, 16);

    ImGui::Separator();
    ImGui::Text("Time: %.4f", static_cast<double>(sim_->time()));
    ImGui::Text("Step: %lld", static_cast<long long>(sim_->step_count()));
    ImGui::Text("Grid: %lld x %lld x %lld",
                static_cast<long long>(cfg_.nx),
                static_cast<long long>(cfg_.ny),
                static_cast<long long>(cfg_.nz));
    Real const e = total_energy_3d(sim_->current(), sim_->previous(),
                                   sim_->wave_speed(), sim_->dt());
    ImGui::Text("Total energy: %.4f", static_cast<double>(e));

    ImGui::Separator();
    int mx = static_cast<int>(cfg_.nx) - 1;
    int my = static_cast<int>(cfg_.ny) - 1;
    int mz = static_cast<int>(cfg_.nz) - 1;
    ImGui::SliderInt("Slice X (yz plane)", &slice_x_, 0, mx);
    ImGui::SliderInt("Slice Y (xz plane)", &slice_y_, 0, my);
    ImGui::SliderInt("Slice Z (xy plane)", &slice_z_, 0, mz);

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
    int m = static_cast<int>(mode_);
    bool changed = false;
    changed |= ImGui::RadioButton("Field u",   &m, 0); ImGui::SameLine();
    changed |= ImGui::RadioButton("Index n",   &m, 1); ImGui::SameLine();
    changed |= ImGui::RadioButton("Energy",    &m, 2);
    if (changed) mode_ = static_cast<OverlayMode>(m);

    if (mode_ == OverlayMode::Field) {
        float v = static_cast<float>(vmax_);
        if (ImGui::SliderFloat("Field scale", &v, 0.01f, 5.0f, "%.2f")) {
            vmax_ = static_cast<Real>(v);
        }
    }
    ImGui::End();
}

void View3D::draw_field(Rectangle area) {
    if (!gpu_ready_) return;

    // Lay the three slices in a 2×2 grid (xy top-left, xz top-right,
    // yz bottom-left). Each panel gets ~half the area.
    float pad = 8.0f;
    float w = (area.width  - pad) * 0.5f;
    float h = (area.height - pad) * 0.5f;

    auto draw_slice = [&](Texture2D const& tex, Rectangle dst, char const* label) {
        Rectangle src{0, 0, static_cast<float>(tex.width),
                            static_cast<float>(tex.height)};
        DrawTexturePro(tex, src, dst, Vector2{0, 0}, 0.0f, WHITE);
        DrawRectangleLinesEx(dst, 1.0f, Color{60, 64, 78, 255});
        DrawText(label, static_cast<int>(dst.x + 6),
                       static_cast<int>(dst.y + 6), 14, RAYWHITE);
    };
    draw_slice(tex_xy_, Rectangle{area.x,            area.y,            w, h}, "XY (z slice)");
    draw_slice(tex_xz_, Rectangle{area.x + w + pad,  area.y,            w, h}, "XZ (y slice)");
    draw_slice(tex_yz_, Rectangle{area.x,            area.y + h + pad,  w, h}, "YZ (x slice)");
}

void View3D::colormap_diverging_(Real v, Real vmax,
                                 unsigned char& r, unsigned char& g, unsigned char& b) {
    if (vmax <= Real{0}) { r = g = b = 0; return; }
    Real t = std::clamp(v / vmax, -1.0_r, 1.0_r);
    if (t >= Real{0}) {
        Real const w = Real{1} - t;
        r = 255;
        g = static_cast<unsigned char>(255 * w);
        b = static_cast<unsigned char>(255 * w);
    } else {
        Real const w = Real{1} + t;
        r = static_cast<unsigned char>(255 * w);
        g = static_cast<unsigned char>(255 * w);
        b = 255;
    }
}

void View3D::colormap_sequential_(Real v, Real vmax,
                                  unsigned char& r, unsigned char& g, unsigned char& b) {
    if (vmax <= Real{0}) { r = g = b = 0; return; }
    Real const t = std::clamp(v / vmax, 0.0_r, 1.0_r);
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
