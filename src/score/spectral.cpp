#include "score/spectral.hpp"

// PocketFFT's heavily-templated code legitimately uses pointers the
// compiler can't always prove non-null; suppress the spurious warnings
// just for this include.
#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wnull-dereference"
#  pragma GCC diagnostic ignored "-Wdouble-promotion"
#  pragma GCC diagnostic ignored "-Wfloat-conversion"
#  pragma GCC diagnostic ignored "-Wconversion"
#endif
#include <pocketfft_hdronly.h>
#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic pop
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>

namespace wavelab {

std::vector<Real> compute_power_spectrum(std::vector<Real> const& series) {
    std::size_t const N = series.size();
    if (N == 0) return {};
    std::size_t const M = N / 2 + 1;

    using cpx = std::complex<Real>;
    std::vector<cpx> spectrum(M);

    pocketfft::shape_t  shape{N};
    pocketfft::stride_t stride_in{static_cast<std::ptrdiff_t>(sizeof(Real))};
    pocketfft::stride_t stride_out{static_cast<std::ptrdiff_t>(sizeof(cpx))};
    pocketfft::shape_t  axes{0};

    pocketfft::r2c<Real>(shape, stride_in, stride_out, axes,
                         /*forward=*/true,
                         series.data(), spectrum.data(),
                         /*fct=*/Real{1});

    std::vector<Real> power(M);
    for (std::size_t k = 0; k < M; ++k) {
        Real const re = spectrum[k].real();
        Real const im = spectrum[k].imag();
        power[k] = re * re + im * im;
    }
    return power;
}

std::vector<Real> compute_frequency_axis(std::size_t n_samples, Real dt) {
    std::size_t const M = n_samples / 2 + 1;
    std::vector<Real> freqs(M);
    Real const df = (n_samples == 0 || dt == Real{0})
        ? Real{0}
        : Real{1} / (dt * static_cast<Real>(n_samples));
    for (std::size_t k = 0; k < M; ++k) {
        freqs[k] = static_cast<Real>(k) * df;
    }
    return freqs;
}

Real spectral_band_energy(std::vector<Real> const& power,
                          std::vector<Real> const& freqs,
                          Real f_lo, Real f_hi) noexcept {
    if (power.size() != freqs.size()) return Real{0};
    Real sum = Real{0};
    for (std::size_t k = 0; k < power.size(); ++k) {
        Real const f = freqs[k];
        if (f >= f_lo && f <= f_hi) sum += power[k];
    }
    return sum;
}

std::vector<Real> spectral_fingerprint_logbins(
        std::vector<Real> const& power,
        std::vector<Real> const& freqs,
        std::size_t nbands,
        Real f_min, Real f_max) {
    std::vector<Real> bins(nbands, Real{0});
    if (nbands == 0 || power.size() != freqs.size()) return bins;
    if (f_min <= Real{0}) f_min = 1e-6_r;
    if (f_max <= f_min) return bins;

    Real const log_lo = std::log(f_min);
    Real const log_hi = std::log(f_max);
    Real const dlog   = (log_hi - log_lo) / static_cast<Real>(nbands);
    if (dlog <= Real{0}) return bins;

    for (std::size_t k = 0; k < power.size(); ++k) {
        Real const f = freqs[k];
        if (f <= Real{0}) continue;
        Real const lf = std::log(f);
        if (lf < log_lo || lf >= log_hi) continue;
        auto b = static_cast<std::size_t>((lf - log_lo) / dlog);
        if (b >= nbands) b = nbands - 1;
        bins[b] += power[k];
    }
    return bins;
}

} // namespace wavelab
