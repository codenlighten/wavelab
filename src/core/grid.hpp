#pragma once
//
// Grid<D>: a uniform Cartesian D-dimensional grid description.
//
// Does NOT own field data — it just describes layout: shape (cell counts),
// spacing (cell size in world units), origin (world coordinate of cell 0).
// Field<T,D> holds a Grid<D> and owns the storage.
//
// Layout convention: row-major (C-style). For D=3, the fastest-varying
// axis is the last (z). linear index = ((i*ny) + j)*nz + k.
//

#include "types.hpp"

#include <cassert>

namespace wavelab {

template <int D>
struct Grid {
    static_assert(D >= 1 && D <= 3, "Grid<D> currently supports D in {1,2,3}");

    IVec<D> shape{};    // number of cells per axis
    Vec<D>  spacing{};  // cell size per axis (world units, e.g. Angstroms)
    Vec<D>  origin{};   // world-coord of cell (0,0,...) lower corner

    // --- Construction -----------------------------------------------------

    constexpr Grid() = default;

    constexpr Grid(IVec<D> shape_, Vec<D> spacing_, Vec<D> origin_ = make_vec<D>(Real{0}))
        : shape(shape_), spacing(spacing_), origin(origin_) {}

    // Uniform-spacing convenience.
    static constexpr Grid uniform(IVec<D> shape_, Real h, Vec<D> origin_ = make_vec<D>(Real{0})) {
        return Grid(shape_, make_vec<D>(h), origin_);
    }

    // --- Queries ----------------------------------------------------------

    constexpr Index num_cells() const noexcept { return product<D>(shape); }

    constexpr bool in_bounds(IVec<D> const& idx) const noexcept {
        for (std::size_t i = 0; i < static_cast<std::size_t>(D); ++i) {
            if (idx[i] < 0 || idx[i] >= shape[i]) return false;
        }
        return true;
    }

    // Row-major linear index. No bounds check (use in_bounds in debug).
    constexpr Index linear_index(IVec<D> const& idx) const noexcept {
        if constexpr (D == 1) {
            return idx[0];
        } else if constexpr (D == 2) {
            return idx[0] * shape[1] + idx[1];
        } else { // D == 3
            return (idx[0] * shape[1] + idx[1]) * shape[2] + idx[2];
        }
    }

    // World coordinate of the *center* of cell `idx`.
    constexpr Vec<D> cell_center(IVec<D> const& idx) const noexcept {
        Vec<D> r{};
        for (std::size_t i = 0; i < static_cast<std::size_t>(D); ++i) {
            r[i] = origin[i] + (static_cast<Real>(idx[i]) + Real{0.5}) * spacing[i];
        }
        return r;
    }

    // World extent (length per axis).
    constexpr Vec<D> extent() const noexcept {
        Vec<D> r{};
        for (std::size_t i = 0; i < static_cast<std::size_t>(D); ++i) {
            r[i] = static_cast<Real>(shape[i]) * spacing[i];
        }
        return r;
    }

    // Smallest spacing (used for Courant: dt <= h/(c*sqrt(D))).
    constexpr Real min_spacing() const noexcept {
        Real h = spacing[0];
        for (std::size_t i = 1; i < static_cast<std::size_t>(D); ++i) {
            if (spacing[i] < h) h = spacing[i];
        }
        return h;
    }
};

} // namespace wavelab
