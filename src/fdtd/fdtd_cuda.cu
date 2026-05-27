// CUDA implementation of FdtdCuda<D>.
//
// Kernels target RTX 6000 Ada (SM 8.9) but the math is plain-vanilla
// 5-point (2D) / 7-point (3D) FDTD stencils with per-cell c, α. Block
// shape is 16x16 in 2D, 8x8x4 in 3D — modest occupancy, no shared-mem
// tiling yet (cache hit-rates on Ada are already high for this stencil).
//
// Device storage rotation mirrors the CPU path: three Real buffers
// swap pointers each step. Host mirrors for current()/previous() are
// pulled on demand via cudaMemcpy DeviceToHost.

#include "fdtd/fdtd_cuda.hpp"

#include "fdtd/stability.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

namespace wavelab {

namespace {

inline void cuda_check_(cudaError_t e, char const* what) {
    if (e != cudaSuccess) {
        throw std::runtime_error(std::string{"FdtdCuda: "} + what + ": "
                                 + cudaGetErrorString(e));
    }
}
#define WL_CUDA_CHECK(call) ::wavelab::cuda_check_((call), #call)

// ---------------------------------------------------------------------
// 2D stencil kernel
// ---------------------------------------------------------------------
__global__ void fdtd_step_2d_kernel(
        Real const* __restrict__ up,
        Real const* __restrict__ uc,
        Real* __restrict__       un,
        Real const* __restrict__ cf,
        Real const* __restrict__ af,
        int   nx, int ny,
        Real  dt, Real inv_dx) {
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < 1 || i >= nx - 1 || j < 1 || j >= ny - 1) return;

    int   idx   = i * ny + j;
    Real  c     = cf[idx];
    Real  alpha = af[idx];
    Real  lam   = c * dt * inv_dx;
    Real  lam2  = lam * lam;
    Real  a_co  = Real(2) - alpha * dt;
    Real  b_co  = Real(1) - alpha * dt;
    Real  lap   = uc[idx + ny] + uc[idx - ny]
                + uc[idx + 1]  + uc[idx - 1]
                - Real(4) * uc[idx];
    un[idx] = a_co * uc[idx] - b_co * up[idx] + lam2 * lap;
}

// Zero the four boundary edges (Dirichlet). Cheap enough to launch a
// single thread per boundary cell rather than per face.
__global__ void zero_boundary_2d_kernel(Real* un, int nx, int ny) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t < ny) {
        un[0 * ny + t]        = Real(0);
        un[(nx - 1) * ny + t] = Real(0);
    }
    if (t < nx) {
        un[t * ny + 0]        = Real(0);
        un[t * ny + (ny - 1)] = Real(0);
    }
}

// ---------------------------------------------------------------------
// 3D stencil kernel (7-point Laplacian)
// ---------------------------------------------------------------------
__global__ void fdtd_step_3d_kernel(
        Real const* __restrict__ up,
        Real const* __restrict__ uc,
        Real* __restrict__       un,
        Real const* __restrict__ cf,
        Real const* __restrict__ af,
        int   nx, int ny, int nz,
        Real  dt, Real inv_dx) {
    int i = blockIdx.z * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < 1 || i >= nx - 1) return;
    if (j < 1 || j >= ny - 1) return;
    if (k < 1 || k >= nz - 1) return;

    int   sy    = nz;             // y-stride
    int   sx    = ny * nz;        // x-stride
    int   idx   = i * sx + j * sy + k;
    Real  c     = cf[idx];
    Real  alpha = af[idx];
    Real  lam   = c * dt * inv_dx;
    Real  lam2  = lam * lam;
    Real  a_co  = Real(2) - alpha * dt;
    Real  b_co  = Real(1) - alpha * dt;
    Real  lap   = uc[idx + sx] + uc[idx - sx]
                + uc[idx + sy] + uc[idx - sy]
                + uc[idx + 1]  + uc[idx - 1]
                - Real(6) * uc[idx];
    un[idx] = a_co * uc[idx] - b_co * up[idx] + lam2 * lap;
}

// Zero the six faces (Dirichlet) for 3D. One small kernel per face pair.
__global__ void zero_boundary_3d_x_kernel(Real* un, int nx, int ny, int nz) {
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= ny || k >= nz) return;
    int sy = nz, sx = ny * nz;
    un[0 * sx + j * sy + k]         = Real(0);
    un[(nx - 1) * sx + j * sy + k]  = Real(0);
}
__global__ void zero_boundary_3d_y_kernel(Real* un, int nx, int ny, int nz) {
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nx || k >= nz) return;
    int sy = nz, sx = ny * nz;
    un[i * sx + 0 * sy + k]         = Real(0);
    un[i * sx + (ny - 1) * sy + k]  = Real(0);
}
__global__ void zero_boundary_3d_z_kernel(Real* un, int nx, int ny, int nz) {
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nx || j >= ny) return;
    int sy = nz, sx = ny * nz;
    un[i * sx + j * sy + 0]         = Real(0);
    un[i * sx + j * sy + (nz - 1)]  = Real(0);
}

} // namespace

// ---------------------------------------------------------------------
// PIMPL definition
// ---------------------------------------------------------------------

template <int D>
struct FdtdCuda<D>::Impl {
    Grid<D>   grid;
    Real      dt;
    Real      t  = Real{0};
    Index     step = 0;

    Real* d_up = nullptr;
    Real* d_uc = nullptr;
    Real* d_un = nullptr;
    Real* d_c  = nullptr;
    Real* d_a  = nullptr;
    std::size_t n_cells_bytes = 0;

    // Host mirrors for read-back accessors.
    mutable Field<Real, D> h_curr;
    mutable Field<Real, D> h_prev;
    mutable bool host_curr_stale = true;
    mutable bool host_prev_stale = true;

    std::vector<std::shared_ptr<Source<D>>> sources;
    std::shared_ptr<BoundaryCondition<D>>   boundary;
    Real cmax = Real{1};

    Impl(Grid<D> g, Medium<D> medium, Real dt_)
        : grid(g),
          dt(dt_),
          h_curr(g, Real{0}),
          h_prev(g, Real{0}) {
        cmax = medium.max_c();
        if (!cfl_satisfied<D>(grid, cmax, dt)) {
            throw std::invalid_argument(
                "FdtdCuda: dt violates CFL bound (Courant condition §3)");
        }
        n_cells_bytes = static_cast<std::size_t>(grid.num_cells()) * sizeof(Real);

        WL_CUDA_CHECK(cudaMalloc(&d_up, n_cells_bytes));
        WL_CUDA_CHECK(cudaMalloc(&d_uc, n_cells_bytes));
        WL_CUDA_CHECK(cudaMalloc(&d_un, n_cells_bytes));
        WL_CUDA_CHECK(cudaMalloc(&d_c,  n_cells_bytes));
        WL_CUDA_CHECK(cudaMalloc(&d_a,  n_cells_bytes));

        WL_CUDA_CHECK(cudaMemset(d_up, 0, n_cells_bytes));
        WL_CUDA_CHECK(cudaMemset(d_uc, 0, n_cells_bytes));
        WL_CUDA_CHECK(cudaMemset(d_un, 0, n_cells_bytes));

        WL_CUDA_CHECK(cudaMemcpy(d_c, medium.c.data(),
                                 n_cells_bytes, cudaMemcpyHostToDevice));
        WL_CUDA_CHECK(cudaMemcpy(d_a, medium.alpha.data(),
                                 n_cells_bytes, cudaMemcpyHostToDevice));
    }

    ~Impl() {
        if (d_up) cudaFree(d_up);
        if (d_uc) cudaFree(d_uc);
        if (d_un) cudaFree(d_un);
        if (d_c)  cudaFree(d_c);
        if (d_a)  cudaFree(d_a);
    }

    void download_(Field<Real, D>& dst, Real const* src_dev) const {
        cuda_check_(cudaMemcpy(dst.data(), src_dev, n_cells_bytes,
                               cudaMemcpyDeviceToHost),
                    "download field");
    }
    void upload_(Real* dst_dev, Field<Real, D> const& src) const {
        cuda_check_(cudaMemcpy(dst_dev, src.data(), n_cells_bytes,
                               cudaMemcpyHostToDevice),
                    "upload field");
    }

    void launch_stencil_() {
        if constexpr (D == 2) {
            int nx = static_cast<int>(grid.shape[0]);
            int ny = static_cast<int>(grid.shape[1]);
            Real inv_dx = Real(1) / grid.spacing[0];
            dim3 block(16, 16);
            dim3 gridDim(
                static_cast<unsigned>((ny + block.x - 1) / block.x),
                static_cast<unsigned>((nx + block.y - 1) / block.y));
            fdtd_step_2d_kernel<<<gridDim, block>>>(
                d_up, d_uc, d_un, d_c, d_a, nx, ny, dt, inv_dx);
            WL_CUDA_CHECK(cudaGetLastError());

            int maxn = std::max(nx, ny);
            int tpb = 256;
            int blocks = (maxn + tpb - 1) / tpb;
            zero_boundary_2d_kernel<<<blocks, tpb>>>(d_un, nx, ny);
            WL_CUDA_CHECK(cudaGetLastError());
        } else if constexpr (D == 3) {
            int nx = static_cast<int>(grid.shape[0]);
            int ny = static_cast<int>(grid.shape[1]);
            int nz = static_cast<int>(grid.shape[2]);
            Real inv_dx = Real(1) / grid.spacing[0];
            dim3 block(8, 8, 4);          // (z fastest, y, x slowest) per kernel
            dim3 gridDim(
                static_cast<unsigned>((nz + block.x - 1) / block.x),
                static_cast<unsigned>((ny + block.y - 1) / block.y),
                static_cast<unsigned>((nx + block.z - 1) / block.z));
            fdtd_step_3d_kernel<<<gridDim, block>>>(
                d_up, d_uc, d_un, d_c, d_a, nx, ny, nz, dt, inv_dx);
            WL_CUDA_CHECK(cudaGetLastError());

            // Dirichlet face zeroing — three face-pair kernels.
            dim3 face_block(16, 16);
            auto blocks_for = [](int a, int b) {
                return dim3(static_cast<unsigned>((a + 15) / 16),
                            static_cast<unsigned>((b + 15) / 16));
            };
            zero_boundary_3d_x_kernel<<<blocks_for(nz, ny), face_block>>>(d_un, nx, ny, nz);
            zero_boundary_3d_y_kernel<<<blocks_for(nz, nx), face_block>>>(d_un, nx, ny, nz);
            zero_boundary_3d_z_kernel<<<blocks_for(ny, nx), face_block>>>(d_un, nx, ny, nz);
            WL_CUDA_CHECK(cudaGetLastError());
        } else {
            static_assert(D == 2 || D == 3,
                "FdtdCuda<D>: D=1 not implemented");
        }
    }

    // Hard-source path: download u_next, let host-side sources mutate,
    // upload back. Acceptable cost for tests / point sources; will move
    // to device-side source kernels in a follow-up.
    void inject_sources_(Real t_new) {
        if (sources.empty()) return;
        Field<Real, D> tmp(grid);
        download_(tmp, d_un);
        for (auto const& s : sources) {
            if (s->active(t_new)) s->inject(tmp, t_new, dt);
        }
        upload_(d_un, tmp);
    }

    void rotate_() {
        // (prev, curr, next) <- (curr, next, prev)
        Real* t = d_up;
        d_up = d_uc;
        d_uc = d_un;
        d_un = t;
        host_curr_stale = true;
        host_prev_stale = true;
    }
};

// ---------------------------------------------------------------------
// Public-API method implementations
// ---------------------------------------------------------------------

template <int D>
FdtdCuda<D>::FdtdCuda(Grid<D> grid, Medium<D> medium, Real dt_)
    : impl_(std::make_unique<Impl>(std::move(grid), std::move(medium), dt_)) {}

template <int D>
FdtdCuda<D>::FdtdCuda(Grid<D> grid, Real c, Real dt_, Real gamma)
    : FdtdCuda(grid, Medium<D>::uniform(grid, c, gamma), dt_) {}

template <int D>
FdtdCuda<D>::~FdtdCuda() = default;

template <int D>
void FdtdCuda<D>::step() {
    impl_->launch_stencil_();
    impl_->inject_sources_(impl_->t + impl_->dt);
    impl_->rotate_();
    impl_->t += impl_->dt;
    ++impl_->step;
}

template <int D>
void FdtdCuda<D>::reset() {
    WL_CUDA_CHECK(cudaMemset(impl_->d_up, 0, impl_->n_cells_bytes));
    WL_CUDA_CHECK(cudaMemset(impl_->d_uc, 0, impl_->n_cells_bytes));
    WL_CUDA_CHECK(cudaMemset(impl_->d_un, 0, impl_->n_cells_bytes));
    impl_->t    = Real{0};
    impl_->step = 0;
    impl_->host_curr_stale = true;
    impl_->host_prev_stale = true;
}

template <int D>
Grid<D> const& FdtdCuda<D>::grid() const noexcept { return impl_->grid; }

template <int D>
Field<Real, D> const& FdtdCuda<D>::current() const noexcept {
    if (impl_->host_curr_stale) {
        impl_->download_(impl_->h_curr, impl_->d_uc);
        impl_->host_curr_stale = false;
    }
    return impl_->h_curr;
}

template <int D>
Field<Real, D> const& FdtdCuda<D>::previous() const noexcept {
    if (impl_->host_prev_stale) {
        impl_->download_(impl_->h_prev, impl_->d_up);
        impl_->host_prev_stale = false;
    }
    return impl_->h_prev;
}

template <int D>
Real FdtdCuda<D>::time() const noexcept { return impl_->t; }
template <int D>
Real FdtdCuda<D>::dt() const noexcept { return impl_->dt; }
template <int D>
Index FdtdCuda<D>::step_count() const noexcept { return impl_->step; }

template <int D>
void FdtdCuda<D>::add_source(std::shared_ptr<Source<D>> s) {
    impl_->sources.push_back(std::move(s));
}
template <int D>
void FdtdCuda<D>::set_boundary(std::shared_ptr<BoundaryCondition<D>> b) {
    impl_->boundary = std::move(b);
    // Boundary objects are currently honored only as a marker — the GPU
    // stencil zeros the frame (Dirichlet) regardless. Non-Dirichlet
    // boundaries on GPU are a follow-up.
}
template <int D>
void FdtdCuda<D>::clear_sources() { impl_->sources.clear(); }

// Explicit instantiations.
template class FdtdCuda<2>;
template class FdtdCuda<3>;

} // namespace wavelab
