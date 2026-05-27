#pragma once
//
// Binder calibration (overview §34).
//
// Build a prototype from K known-binder fingerprints:
//
//   F̄_binder = (1/K) Σ_k F_k
//
// Then for a new candidate fingerprint, compute cosine similarity
// against the prototype's spectral vector:
//
//   S_known = (F_new · F̄_binder) / (||F_new|| · ||F̄_binder||)
//
// Scalar features are averaged separately and exposed for downstream
// inspection / weighted combination. The cosine similarity here uses
// just the spectral component because (a) spectra are unitless and
// directly comparable, (b) scalars often live on incompatible scales.
//

#include "core/types.hpp"
#include "io/fingerprint.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace wavelab {

struct BinderPrototype {
    std::unordered_map<std::string, Real> scalar_mean;   // averaged scalars
    std::vector<Real>                     spectral_mean; // averaged spectrum
    std::vector<Real>                     spectral_freqs;
    Index                                 sample_count = 0;
};

// Aggregate K fingerprints into a prototype. All inputs must agree on
// spectral length; throws std::invalid_argument otherwise.
BinderPrototype build_prototype(std::vector<Fingerprint> const& binders);

// Cosine similarity between `candidate.spectral` and `proto.spectral_mean`.
// Returns 0 if either vector is zero-norm.
Real correlate_with_prototype(Fingerprint const& candidate,
                              BinderPrototype const& proto);

// Convenience: scalar agreement = exp(-||c - p||/scale) over named keys.
// Returns 1 when scalars match exactly, decays as they diverge. `scale`
// controls sensitivity; defaults to a generic 1.0.
Real scalar_agreement(Fingerprint const& candidate,
                      BinderPrototype const& proto,
                      Real scale = 1.0_r);

} // namespace wavelab
