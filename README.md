# wavelab

Wave-based molecular geometry engine. Treats a protein pocket as an
optical medium and probes it with a simulated light field (FDTD) to
extract geometric, energetic, and spectral features for docking.

See [`overview.md`](overview.md) for the full equation map.

## Status

**Phase 0** — project skeleton. Core types (`Grid<D>`, `Field<T,D>`),
CMake build, doctest unit tests, empty GUI window via raylib + ImGui.

The phased implementation plan lives at
`~/.claude/plans/composed-snuggling-corbato.md`.

## Building

### Prerequisites

- CMake ≥ 3.24
- C++20 compiler (GCC 13+, Clang 16+, MSVC 19.30+)
- For the GUI app on Linux/X11 (skip with `-DWAVELAB_BUILD_GUI=OFF`):

```sh
sudo apt install libgl1-mesa-dev libx11-dev libxrandr-dev \
    libxinerama-dev libxcursor-dev libxi-dev libxkbcommon-dev
```

### Configure + build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Headless / core-only build (no X11 deps required):

```sh
cmake -S . -B build -DWAVELAB_BUILD_GUI=OFF
cmake --build build -j
```

### Run

```sh
./build/bin/wavecli            # headless smoke test
./build/bin/wavelab            # GUI window (Phase 0: ImGui hello)
```

## Build options

| Option                       | Default | Effect                                       |
| ---------------------------- | ------- | -------------------------------------------- |
| `WAVELAB_BUILD_TESTS`        | ON      | Compile + register unit tests with CTest     |
| `WAVELAB_BUILD_GUI`          | ON      | Compile `wavelab` (raylib + ImGui)           |
| `WAVELAB_BUILD_CLI`          | ON      | Compile `wavecli` (headless runner)          |
| `WAVELAB_BUILD_BENCH`        | OFF     | Compile benchmark harness                    |
| `WAVELAB_DOUBLE_PRECISION`   | OFF     | Field scalar = `double` instead of `float`   |
| `WAVELAB_USE_OPENMP`         | ON      | Link OpenMP for parallel kernels             |
| `WAVELAB_SANITIZE`           | OFF     | Enable ASan + UBSan in Debug builds          |
| `WAVELAB_WERROR`             | ON      | Treat warnings as errors                     |

## Layout

```
src/      core/, fdtd/, medium/, source/, molecule/, score/, io/, viz/
apps/     wavelab/ (GUI), wavecli/ (headless)
tests/    unit/, integration/
bench/    perf harness
data/     configs/, pdb/, synthetic/
```
