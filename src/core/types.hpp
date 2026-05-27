#pragma once
//
// Fundamental type aliases used across the wavelab engine.
//
// `Real` is the scalar floating-point type for all field values. It is
// configurable at build time via the CMake option WAVELAB_DOUBLE_PRECISION
// (default OFF -> float). Using a single alias lets us swap precision
// engine-wide without touching call sites.
//
// `Vec<D>` and `IVec<D>` are fixed-size D-dimensional vectors backed by
// std::array. They are deliberately small, trivially-copyable value
// types — passed by value, no heap, no virtuals.
//

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace wavelab {

#if WAVELAB_DOUBLE_PRECISION
using Real = double;
#else
using Real = float;
#endif

using Index = std::int64_t;

template <int D>
using Vec = std::array<Real, static_cast<std::size_t>(D)>;

template <int D>
using IVec = std::array<Index, static_cast<std::size_t>(D)>;

static_assert(std::is_trivially_copyable_v<Vec<3>>);
static_assert(std::is_trivially_copyable_v<IVec<3>>);

// --- Small helpers --------------------------------------------------------
//
// Loop indices use std::size_t so that indexing std::array (which takes
// size_t) does not trigger -Wsign-conversion.

template <int D>
constexpr Real dot(Vec<D> const& a, Vec<D> const& b) noexcept {
    Real s{};
    for (std::size_t i = 0; i < static_cast<std::size_t>(D); ++i) s += a[i] * b[i];
    return s;
}

template <int D>
constexpr Index product(IVec<D> const& v) noexcept {
    Index p = 1;
    for (std::size_t i = 0; i < static_cast<std::size_t>(D); ++i) p *= v[i];
    return p;
}

template <int D>
constexpr IVec<D> make_ivec(Index v) noexcept {
    IVec<D> r{};
    for (std::size_t i = 0; i < static_cast<std::size_t>(D); ++i) r[i] = v;
    return r;
}

template <int D>
constexpr Vec<D> make_vec(Real v) noexcept {
    Vec<D> r{};
    for (std::size_t i = 0; i < static_cast<std::size_t>(D); ++i) r[i] = v;
    return r;
}

// --- User-defined literal for Real scalars --------------------------------
//
// Use `0.1_r`, `1.5_r` etc. to write a literal whose type matches the
// engine's configured precision. Avoids -Wfloat-conversion noise on
// `Real{0.1}` when Real == float.

inline namespace literals {

constexpr Real operator""_r(long double v) noexcept {
    return static_cast<Real>(v);
}

constexpr Real operator""_r(unsigned long long v) noexcept {
    return static_cast<Real>(v);
}

} // namespace literals

} // namespace wavelab
