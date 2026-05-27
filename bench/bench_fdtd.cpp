// bench_fdtd — wall-clock benchmark of FdtdCpuOmp<2> vs FdtdCuda<2>.
//
// No google-benchmark dependency. Times only the step() loop after a
// warm-up to factor out one-time allocation costs and JIT-compile +
// driver-init on the GPU side.
//
// Source injection is deliberately disabled in this benchmark: it costs
// a device⇄host roundtrip per step in the current CUDA path and would
// dominate at small grids, masking raw kernel throughput. The "GPU with
// sources" cost is a separate concern we'll measure in a follow-up.

#include "core/types.hpp"
#include "fdtd/fdtd_cpu_omp.hpp"
#include "fdtd/stability.hpp"
#include "medium/boundary.hpp"
#include "medium/medium.hpp"

#if WAVELAB_BENCH_HAVE_CUDA
#  include "fdtd/fdtd_cuda.hpp"
#endif

#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

template <typename Stepper>
double time_steps_seconds(Stepper& sim, int n_steps) {
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    sim.run(n_steps);
    // Force any pending async work (CUDA in particular). For CPU this
    // is a noop on the value path.
    (void)sim.current();   // forces device→host sync on GPU
    auto t1 = clock::now();
    return std::chrono::duration<double>(t1 - t0).count();
}

void seed_pulse(FdtdCpuOmp<2>& sim) {
    auto& curr = sim.current_mut();
    auto& prev = sim.previous_mut();
    Index const cx = sim.grid().shape[0] / 2;
    Index const cy = sim.grid().shape[1] / 2;
    curr(cx, cy) = 1.0_r;
    prev(cx, cy) = 1.0_r;
}

void bench_size(Index nx, Index ny) {
    Real const dx = 0.1_r;
    Real const c0 = 1.0_r;
    auto grid = Grid<2>::uniform(IVec<2>{nx, ny}, dx);
    Real const dt = 0.4_r * cfl_dt_max<2>(grid, c0);

    constexpr int warmup = 20;
    constexpr int measured = 200;

    std::printf("  grid %4lld x %4lld:\n",
                static_cast<long long>(nx), static_cast<long long>(ny));

    // --- CPU ---
    {
        FdtdCpuOmp<2> cpu(grid, c0, dt);
        cpu.set_boundary(std::make_shared<Dirichlet2D>());
        seed_pulse(cpu);
        cpu.run(warmup);
        double secs = time_steps_seconds(cpu, measured);
        double cells = static_cast<double>(nx) * static_cast<double>(ny);
        double glups = (cells * measured) / secs / 1e9;
        std::printf("    CPU+OMP : %7.3f s   %.3f GLups\n", secs, glups);
    }

#if WAVELAB_BENCH_HAVE_CUDA
    // --- GPU ---
    {
        FdtdCuda<2> gpu(grid, c0, dt);
        gpu.set_boundary(std::make_shared<Dirichlet2D>());
        // No host-mutate path on GPU — start from zero, kernel-warm
        gpu.run(warmup);
        double secs = time_steps_seconds(gpu, measured);
        double cells = static_cast<double>(nx) * static_cast<double>(ny);
        double glups = (cells * measured) / secs / 1e9;
        std::printf("    CUDA    : %7.3f s   %.3f GLups\n", secs, glups);
    }
#else
    std::printf("    CUDA    : (build with -DWAVELAB_BUILD_CUDA=ON)\n");
#endif
}

} // namespace

int main() {
    std::puts("wavelab FDTD benchmark — 2D stencil throughput");
    std::puts("(no sources / Dirichlet boundary; raw kernel throughput)");
    std::puts("");

    Index const sizes[] = {256, 512, 1024};
    for (Index n : sizes) bench_size(n, n);

    return 0;
}
