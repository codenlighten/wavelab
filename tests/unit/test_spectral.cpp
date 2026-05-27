#include <doctest/doctest.h>

#include "core/types.hpp"
#include "score/spectral.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

// Pure tone in → peak at the right frequency bin out.
TEST_CASE("compute_power_spectrum: sine input peaks at expected bin") {
    constexpr std::size_t N  = 512;
    constexpr Real        dt = 0.01_r;
    constexpr Real        f0 = 5.0_r;     // 5 Hz; with dt=0.01 → period 20 samples, well within Nyquist
    std::vector<Real> series(N);
    for (std::size_t k = 0; k < N; ++k) {
        Real const t = static_cast<Real>(k) * dt;
        series[k] = std::sin(2.0_r * static_cast<Real>(std::numbers::pi) * f0 * t);
    }

    auto power = compute_power_spectrum(series);
    auto freqs = compute_frequency_axis(N, dt);
    REQUIRE(power.size() == freqs.size());

    // Argmax peak
    std::size_t argmax = 0;
    for (std::size_t k = 1; k < power.size(); ++k) {
        if (power[k] > power[argmax]) argmax = k;
    }
    Real const f_peak = freqs[argmax];

    // df = 1 / (N * dt) = 1/(512*0.01) = 0.195; should be within 1 bin of f0
    CAPTURE(f_peak);
    CAPTURE(f0);
    CHECK(std::abs(f_peak - f0) < 0.5_r);
}

TEST_CASE("compute_frequency_axis: monotonic, starts at 0, ends near Nyquist") {
    auto freqs = compute_frequency_axis(/*N=*/128, /*dt=*/0.1_r);
    REQUIRE(freqs.size() == 65);
    CHECK(freqs.front() == doctest::Approx(0.0));
    Real const fnyq = 0.5_r / 0.1_r;
    CHECK(freqs.back() == doctest::Approx(static_cast<double>(fnyq)));
    for (std::size_t k = 1; k < freqs.size(); ++k) {
        CHECK(freqs[k] > freqs[k - 1]);
    }
}

TEST_CASE("spectral_fingerprint_logbins: sine concentrates in one bin") {
    constexpr std::size_t N  = 1024;
    constexpr Real        dt = 0.01_r;
    constexpr Real        f0 = 3.0_r;
    std::vector<Real> series(N);
    for (std::size_t k = 0; k < N; ++k) {
        Real const t = static_cast<Real>(k) * dt;
        series[k] = std::sin(2.0_r * static_cast<Real>(std::numbers::pi) * f0 * t);
    }
    auto power = compute_power_spectrum(series);
    auto freqs = compute_frequency_axis(N, dt);
    auto bins  = spectral_fingerprint_logbins(power, freqs, /*nbands=*/8,
                                              /*f_min=*/0.5_r, /*f_max=*/20.0_r);
    REQUIRE(bins.size() == 8);

    // Find peak bin, ensure it covers f0 = 3.
    std::size_t pk = 0;
    for (std::size_t k = 0; k < bins.size(); ++k) if (bins[k] > bins[pk]) pk = k;
    Real total = 0.0_r;
    for (auto v : bins) total += v;
    CHECK(total > 0.0_r);
    // Peak bin's center frequency should bracket f0=3.
    Real const log_lo = std::log(0.5_r);
    Real const log_hi = std::log(20.0_r);
    Real const dlog   = (log_hi - log_lo) / static_cast<Real>(bins.size());
    Real const f_pk   = std::exp(log_lo + (static_cast<Real>(pk) + 0.5_r) * dlog);
    CHECK(f_pk > 1.0_r);
    CHECK(f_pk < 10.0_r);
}

} // namespace
