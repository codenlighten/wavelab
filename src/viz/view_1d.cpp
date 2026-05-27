#include "viz/view_1d.hpp"

#include "fdtd/stability.hpp"
#include "score/energy.hpp"

#include <imgui.h>
#include <raylib.h>

#include <cmath>
#include <cstdio>

namespace wavelab {

View1D::View1D(Config cfg) : cfg_(cfg) {
    rebuild_sim_();
}

void View1D::rebuild_sim_() {
    auto grid = Grid<1>::uniform(IVec<1>{cfg_.nx}, cfg_.dx);
    Real const dt = cfg_.cfl_safety * cfl_dt_max<1>(grid, cfg_.c);
    sim_ = std::make_unique<FdtdCpuOmp<1>>(grid, cfg_.c, dt, cfg_.damping);

    source_ = std::make_shared<HarmonicSource<1>>(
        IVec<1>{cfg_.source_cell}, cfg_.source_amp, cfg_.source_freq);
    sim_->add_source(source_);

    install_boundary_();
}

void View1D::install_boundary_() {
    switch (boundary_kind_) {
        case BoundaryKind::Dirichlet:
            sim_->set_boundary(std::make_shared<Dirichlet1D>()); break;
        case BoundaryKind::Neumann:
            sim_->set_boundary(std::make_shared<Neumann1D>()); break;
        case BoundaryKind::Mur:
            sim_->set_boundary(std::make_shared<Mur1D>()); break;
    }
}

void View1D::update() {
    if (!playing_) return;
    sim_->run(static_cast<Index>(steps_per_frame_));
}

void View1D::draw_controls() {
    ImGui::Begin("View1D — controls");

    ImGui::Text("Phase 1: 1D wave engine");
    ImGui::Separator();

    // Playback
    if (ImGui::Button(playing_ ? "Pause" : "Play")) playing_ = !playing_;
    ImGui::SameLine();
    if (ImGui::Button("Reset")) {
        sim_->reset();
    }
    ImGui::SameLine();
    if (ImGui::Button("Rebuild")) {
        rebuild_sim_();
    }
    ImGui::SliderInt("Steps / frame", &steps_per_frame_, 1, 64);

    ImGui::Separator();
    ImGui::Text("Time:  %.4f", static_cast<double>(sim_->time()));
    ImGui::Text("Step:  %lld", static_cast<long long>(sim_->step_count()));
    ImGui::Text("dt:    %.5f", static_cast<double>(sim_->dt()));
    ImGui::Text("CFL λ: %.3f", static_cast<double>(courant_number<1>(
        sim_->wave_speed(), sim_->dt(), sim_->grid().spacing[0])));

    Real const e = total_energy_conserved_1d(
        sim_->current(), sim_->previous(), sim_->wave_speed(), sim_->dt());
    ImGui::Text("Energy (conserved): %.4f", static_cast<double>(e));

    // Source params — live editable
    ImGui::Separator();
    ImGui::Text("Source (HarmonicSource at i=%lld)",
                static_cast<long long>(cfg_.source_cell));
    float freq = static_cast<float>(source_->frequency());
    float amp  = static_cast<float>(source_->amplitude());
    if (ImGui::SliderFloat("Frequency", &freq, 0.1f, 20.0f, "%.2f Hz")) {
        source_->set_frequency(static_cast<Real>(freq));
    }
    if (ImGui::SliderFloat("Amplitude", &amp, 0.0f, 2.0f, "%.2f")) {
        source_->set_amplitude(static_cast<Real>(amp));
    }

    // Boundary
    ImGui::Separator();
    int b = static_cast<int>(boundary_kind_);
    bool changed = false;
    changed |= ImGui::RadioButton("Dirichlet", &b, 0); ImGui::SameLine();
    changed |= ImGui::RadioButton("Neumann",   &b, 1); ImGui::SameLine();
    changed |= ImGui::RadioButton("Mur (abs)", &b, 2);
    if (changed) {
        boundary_kind_ = static_cast<BoundaryKind>(b);
        install_boundary_();
    }

    // View axis
    ImGui::Separator();
    float yr[2] = {static_cast<float>(ymin_), static_cast<float>(ymax_)};
    if (ImGui::SliderFloat2("y range", yr, -5.0f, 5.0f, "%.2f")) {
        ymin_ = static_cast<Real>(yr[0]);
        ymax_ = static_cast<Real>(yr[1]);
    }

    ImGui::End();
}

void View1D::draw_field(Rectangle area) const {
    auto const& u = sim_->current();
    Index const  n = u.grid().shape[0];
    if (n < 2) return;

    // Background
    DrawRectangleRec(area, Color{14, 16, 22, 255});
    DrawRectangleLinesEx(area, 1.0f, Color{60, 64, 78, 255});

    // Zero axis
    Real const yspan = ymax_ - ymin_;
    if (yspan <= Real{0}) return;
    float const y0 = area.y + area.height
                   - static_cast<float>((Real{0} - ymin_) / yspan)
                   * area.height;
    DrawLine(static_cast<int>(area.x),
             static_cast<int>(y0),
             static_cast<int>(area.x + area.width),
             static_cast<int>(y0),
             Color{60, 64, 78, 255});

    // Field polyline
    auto y_screen = [&](Real v) {
        Real const t = (v - ymin_) / yspan;
        return area.y + area.height - static_cast<float>(t) * area.height;
    };
    auto x_screen = [&](Index i) {
        return area.x + (static_cast<float>(i) / static_cast<float>(n - 1)) * area.width;
    };

    Color const line = Color{114, 188, 255, 255};
    float prev_x = x_screen(0);
    float prev_y = y_screen(u(Index{0}));
    for (Index i = 1; i < n; ++i) {
        float const cx = x_screen(i);
        float const cy = y_screen(u(i));
        DrawLine(static_cast<int>(prev_x), static_cast<int>(prev_y),
                 static_cast<int>(cx),     static_cast<int>(cy),
                 line);
        prev_x = cx;
        prev_y = cy;
    }
}

} // namespace wavelab
