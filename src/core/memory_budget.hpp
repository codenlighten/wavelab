#pragma once
//
// Memory-budget estimator for an FDTD scene.
//
// Three rotating Field<Real,D> buffers (u_prev, u_curr, u_next) plus
// two Medium fields (c, alpha) = 5 grids of `sizeof(Real) · N_cells`
// bytes. Add any caller-side density / charge / hydro grids on top.
//
// Use this BEFORE constructing the stepper to fail fast — a 512³ grid
// at single precision is 4·5·128 MB = ~2.5 GB just for the engine, and
// the user almost certainly didn't mean to allocate that on a laptop.
//

#include "core/grid.hpp"
#include "core/types.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>

namespace wavelab {

struct MemoryEstimate {
    std::size_t n_cells = 0;
    std::size_t bytes_per_field = 0;
    std::size_t n_fields = 5;          // u_prev, u_curr, u_next, c, alpha
    std::size_t total_bytes = 0;
};

template <int D>
inline MemoryEstimate estimate_engine_memory(Grid<D> const& grid,
                                             std::size_t n_extra_fields = 0) {
    MemoryEstimate e;
    e.n_cells         = static_cast<std::size_t>(grid.num_cells());
    e.bytes_per_field = e.n_cells * sizeof(Real);
    e.n_fields        = 5u + n_extra_fields;
    e.total_bytes     = e.bytes_per_field * e.n_fields;
    return e;
}

inline std::string format_bytes(std::size_t b) {
    char buf[64];
    if (b < 1024ull)              std::snprintf(buf, sizeof(buf), "%zu B", b);
    else if (b < 1024ull * 1024)  std::snprintf(buf, sizeof(buf), "%.1f KB", static_cast<double>(b) / 1024.0);
    else if (b < (1ull << 30))    std::snprintf(buf, sizeof(buf), "%.1f MB", static_cast<double>(b) / (1024.0 * 1024.0));
    else                          std::snprintf(buf, sizeof(buf), "%.2f GB", static_cast<double>(b) / (1024.0 * 1024.0 * 1024.0));
    return buf;
}

// Throws std::runtime_error if estimate exceeds `cap_bytes`.
template <int D>
inline void enforce_memory_budget(Grid<D> const& grid, std::size_t cap_bytes,
                                  std::size_t n_extra_fields = 0) {
    auto est = estimate_engine_memory<D>(grid, n_extra_fields);
    if (est.total_bytes > cap_bytes) {
        throw std::runtime_error(
            "engine memory estimate " + format_bytes(est.total_bytes)
            + " exceeds cap " + format_bytes(cap_bytes));
    }
}

} // namespace wavelab
