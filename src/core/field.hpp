#pragma once
//
// Field<T,D>: dense storage of scalar value T over a D-dimensional grid.
//
// Owns its storage. Row-major layout (matches Grid<D>::linear_index).
// operator() accepts D integer indices (variadic) or an IVec<D>.
// at(...) is bounds-checked; operator() is not (debug assert only).
//
// Storage backend: std::vector<T> for Phase 0/1. Switching to over-
// aligned storage for SIMD is a Phase 7 task and won't change this API.
//

#include "grid.hpp"
#include "types.hpp"

#include <algorithm>
#include <cassert>
#include <span>
#include <stdexcept>
#include <vector>

namespace wavelab {

template <typename T, int D>
class Field {
public:
    using value_type = T;
    static constexpr int dim = D;

    Field() = default;

    explicit Field(Grid<D> grid, T init = T{})
        : grid_(grid), data_(static_cast<std::size_t>(grid.num_cells()), init) {}

    // --- Grid + size -----------------------------------------------------

    Grid<D> const& grid() const noexcept { return grid_; }
    Index           size() const noexcept { return static_cast<Index>(data_.size()); }

    // --- Raw access (for hot loops / interop) ---------------------------

    T*       data()       noexcept { return data_.data(); }
    T const* data() const noexcept { return data_.data(); }

    std::span<T>       view()       noexcept { return {data_.data(), data_.size()}; }
    std::span<T const> view() const noexcept { return {data_.data(), data_.size()}; }

    // --- Indexed access -------------------------------------------------

    T& operator()(IVec<D> const& idx) noexcept {
        assert(grid_.in_bounds(idx));
        return data_[static_cast<std::size_t>(grid_.linear_index(idx))];
    }
    T const& operator()(IVec<D> const& idx) const noexcept {
        assert(grid_.in_bounds(idx));
        return data_[static_cast<std::size_t>(grid_.linear_index(idx))];
    }

    // Variadic convenience: f(i), f(i,j), f(i,j,k).
    template <typename... Args>
        requires (sizeof...(Args) == static_cast<std::size_t>(D))
              && (std::is_integral_v<std::decay_t<Args>> && ...)
    T& operator()(Args... a) noexcept {
        IVec<D> idx{ static_cast<Index>(a)... };
        return (*this)(idx);
    }

    template <typename... Args>
        requires (sizeof...(Args) == static_cast<std::size_t>(D))
              && (std::is_integral_v<std::decay_t<Args>> && ...)
    T const& operator()(Args... a) const noexcept {
        IVec<D> idx{ static_cast<Index>(a)... };
        return (*this)(idx);
    }

    // Bounds-checked access; throws std::out_of_range.
    T& at(IVec<D> const& idx) {
        if (!grid_.in_bounds(idx)) throw std::out_of_range("Field::at: index out of bounds");
        return data_[static_cast<std::size_t>(grid_.linear_index(idx))];
    }
    T const& at(IVec<D> const& idx) const {
        if (!grid_.in_bounds(idx)) throw std::out_of_range("Field::at: index out of bounds");
        return data_[static_cast<std::size_t>(grid_.linear_index(idx))];
    }

    // --- Bulk ops -------------------------------------------------------

    void fill(T v) noexcept { std::fill(data_.begin(), data_.end(), v); }
    void zero()    noexcept { fill(T{}); }

    void swap(Field& other) noexcept {
        std::swap(grid_, other.grid_);
        data_.swap(other.data_);
    }

private:
    Grid<D>        grid_{};
    std::vector<T> data_{};
};

} // namespace wavelab
