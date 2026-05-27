// wavelab — interactive GUI front-end for the wave-based molecular
// geometry engine.
//
// Top-level menu toggles between View1D (Phase 1 — 1D wave engine) and
// View2D (Phase 2 — 2D engine with heterogeneous medium + PML).
//

#include "viz/view_1d.hpp"
#include "viz/view_2d.hpp"

#include <imgui.h>
#include <raylib.h>
#include <rlImGui.h>

#include <memory>

namespace {

constexpr int kInitialWidth  = 1400;
constexpr int kInitialHeight = 900;

enum class ViewMode { Plot1D, Heatmap2D };

} // namespace

int main() {
    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_VSYNC_HINT | FLAG_MSAA_4X_HINT);
    InitWindow(kInitialWidth, kInitialHeight, "wavelab");
    SetTargetFPS(60);
    rlImGuiSetup(true);

    ViewMode mode = ViewMode::Heatmap2D;
    std::unique_ptr<wavelab::View1D> view1d;
    std::unique_ptr<wavelab::View2D> view2d = std::make_unique<wavelab::View2D>();

    while (!WindowShouldClose()) {
        if (mode == ViewMode::Plot1D) {
            if (!view1d) view1d = std::make_unique<wavelab::View1D>();
            view1d->update();
        } else {
            if (!view2d) view2d = std::make_unique<wavelab::View2D>();
            view2d->update();
        }

        BeginDrawing();
        ClearBackground(Color{20, 22, 28, 255});

        float const margin = 16.0f;
        float const ctrl_w = 360.0f;
        Rectangle field_area{
            ctrl_w + margin,
            margin + 32.0f, // leave space for menu bar
            static_cast<float>(GetScreenWidth())  - ctrl_w - 2.0f * margin,
            static_cast<float>(GetScreenHeight()) - 2.0f * margin - 32.0f
        };
        if (mode == ViewMode::Plot1D) view1d->draw_field(field_area);
        else                          view2d->draw_field(field_area);

        rlImGuiBegin();

        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("View")) {
                if (ImGui::MenuItem("1D (Phase 1)", nullptr, mode == ViewMode::Plot1D)) {
                    mode = ViewMode::Plot1D;
                }
                if (ImGui::MenuItem("2D (Phase 2)", nullptr, mode == ViewMode::Heatmap2D)) {
                    mode = ViewMode::Heatmap2D;
                }
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }

        if (mode == ViewMode::Plot1D) view1d->draw_controls();
        else                          view2d->draw_controls();

        rlImGuiEnd();
        EndDrawing();
    }

    view1d.reset();
    view2d.reset();

    rlImGuiShutdown();
    CloseWindow();
    return 0;
}
