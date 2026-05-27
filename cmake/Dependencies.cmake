# Third-party dependencies via FetchContent.
# Downloads go under build/_deps/ on first configure.

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

# ------------------------------------------------------------------
# doctest (testing) — header-only, single file
# ------------------------------------------------------------------
if(WAVELAB_BUILD_TESTS)
    FetchContent_Declare(doctest
        GIT_REPOSITORY https://github.com/doctest/doctest.git
        GIT_TAG v2.4.11
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(doctest)
endif()

# ------------------------------------------------------------------
# spdlog (logging)
# ------------------------------------------------------------------
set(SPDLOG_BUILD_SHARED OFF CACHE BOOL "" FORCE)
set(SPDLOG_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.14.1
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(spdlog)

# ------------------------------------------------------------------
# tomlplusplus (config files)
# ------------------------------------------------------------------
FetchContent_Declare(tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG v3.4.0
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(tomlplusplus)

# ------------------------------------------------------------------
# Eigen (small linear algebra; NOT the FDTD grid)
# ------------------------------------------------------------------
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(eigen)

# ------------------------------------------------------------------
# PocketFFT (header-only, single header)
# ------------------------------------------------------------------
FetchContent_Declare(pocketfft
    GIT_REPOSITORY https://github.com/mreineck/pocketfft.git
    GIT_TAG cpp
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(pocketfft)
add_library(pocketfft INTERFACE)
target_include_directories(pocketfft INTERFACE ${pocketfft_SOURCE_DIR})

# ------------------------------------------------------------------
# raylib + Dear ImGui + rlImGui (GUI stack — optional)
# ------------------------------------------------------------------
if(WAVELAB_BUILD_GUI)
    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(BUILD_GAMES OFF CACHE BOOL "" FORCE)
    set(CUSTOMIZE_BUILD ON CACHE BOOL "" FORCE)
    set(SUPPORT_MODULE_RAUDIO OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(raylib
        GIT_REPOSITORY https://github.com/raysan5/raylib.git
        GIT_TAG 5.5
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(raylib)

    FetchContent_Declare(imgui
        GIT_REPOSITORY https://github.com/ocornut/imgui.git
        GIT_TAG v1.91.5
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(imgui)

    add_library(imgui STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_demo.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp)
    target_include_directories(imgui PUBLIC ${imgui_SOURCE_DIR})

    FetchContent_Declare(rlimgui
        GIT_REPOSITORY https://github.com/raylib-extras/rlImGui.git
        GIT_TAG main
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(rlimgui)

    add_library(rlimgui STATIC ${rlimgui_SOURCE_DIR}/rlImGui.cpp)
    target_include_directories(rlimgui PUBLIC ${rlimgui_SOURCE_DIR})
    target_link_libraries(rlimgui PUBLIC raylib imgui)
endif()

# ------------------------------------------------------------------
# pybind11 (Python bindings)
# ------------------------------------------------------------------
if(WAVELAB_BUILD_PYTHON)
    set(PYBIND11_FINDPYTHON ON CACHE BOOL "" FORCE)
    FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.13.6
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(pybind11)
endif()
