#pragma once
//
// Wave-field similarity / difference scoring (overview §18, §26).
//
//   D_wave         = ||W_a - W_b||₂
//   S_wave         = exp(-η · D_wave)            (overall similarity)
//   S_complementarity = (W_a · W_b) / (||W_a|| · ||W_b||)   (cosine, §26)
//
// Both `W_a` and `W_b` must share the same grid. We use the raw field
// values (u(x,t)); callers can equally pass energy-density fields if
// they want amplitude-invariant scoring.
//

#include "core/field.hpp"
#include "core/types.hpp"

#include <cmath>
#include <stdexcept>

namespace wavelab {

template <int D>
inline Real wave_difference(Field<Real, D> const& a,
                            Field<Real, D> const& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("wave_difference: field sizes differ");
    }
    Index const n = a.size();
    Real sum_sq = Real{0};
    Real const* pa = a.data();
    Real const* pb = b.data();
    #pragma omp parallel for reduction(+:sum_sq) schedule(static)
    for (Index k = 0; k < n; ++k) {
        Real const d = pa[k] - pb[k];
        sum_sq += d * d;
    }
    return std::sqrt(sum_sq);
}

template <int D>
inline Real wave_similarity(Field<Real, D> const& a,
                            Field<Real, D> const& b,
                            Real eta = 1.0_r) {
    return std::exp(-eta * wave_difference<D>(a, b));
}

// Cosine similarity (overview §26 shape complementarity). Returns 1 for
// identical fields, 0 for orthogonal, -1 for anti-correlated.
template <int D>
inline Real shape_complementarity(Field<Real, D> const& a,
                                  Field<Real, D> const& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("shape_complementarity: size mismatch");
    }
    Index const n = a.size();
    Real const* pa = a.data();
    Real const* pb = b.data();
    Real dot = Real{0};
    Real na  = Real{0};
    Real nb  = Real{0};
    #pragma omp parallel for reduction(+:dot,na,nb) schedule(static)
    for (Index k = 0; k < n; ++k) {
        dot += pa[k] * pb[k];
        na  += pa[k] * pa[k];
        nb  += pb[k] * pb[k];
    }
    Real const denom = std::sqrt(na) * std::sqrt(nb);
    if (denom == Real{0}) return Real{0};
    return dot / denom;
}

} // namespace wavelab
