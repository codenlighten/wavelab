#pragma once
//
// Spectral fingerprinting (overview §22, §29).
//
// Workflow:
//   1. Pick one or more probe cells in the simulation domain.
//   2. Record u(probe, t) into a time series during the run.
//   3. After the run, compute the one-sided power spectrum P(ω) via FFT.
//   4. Concatenate per-probe spectra into a fingerprint vector
//      F_spectral.
//
// Multi-frequency sweeps (§29): caller runs N independent sims at
// different source frequencies, collects N fingerprints, combines.
//
// Implementation uses PocketFFT (real-to-complex, length = series size).
//

#include "core/field.hpp"
#include "core/types.hpp"

#include <cstddef>
#include <vector>

namespace wavelab {

// Append u(probe) to `series` — call once per step during the sim run.
template <int D>
inline void probe_record(Field<Real, D> const& u, IVec<D> probe,
                         std::vector<Real>& series) {
    series.push_back(u(probe));
}

// One-sided power spectrum P[k] = |U[k]|² for k = 0..N/2.
// Input length N, output length N/2 + 1.
std::vector<Real> compute_power_spectrum(std::vector<Real> const& series);

// Frequency axis in Hz for a power spectrum computed from N real samples
// taken at uniform spacing `dt`. Output length N/2 + 1.
std::vector<Real> compute_frequency_axis(std::size_t n_samples, Real dt);

// Bandpass integration: total spectral energy in [f_lo, f_hi]. Useful
// to reduce a multi-bin spectrum to a few interpretable features.
Real spectral_band_energy(std::vector<Real> const& power,
                          std::vector<Real> const& freqs,
                          Real f_lo, Real f_hi) noexcept;

// Convenience: collapse a full power spectrum into a small fixed-size
// fingerprint by summing into `nbands` logarithmic frequency bins
// covering [f_min, f_max]. Returns length-`nbands` vector.
std::vector<Real> spectral_fingerprint_logbins(
    std::vector<Real> const& power,
    std::vector<Real> const& freqs,
    std::size_t nbands,
    Real f_min, Real f_max);

} // namespace wavelab
