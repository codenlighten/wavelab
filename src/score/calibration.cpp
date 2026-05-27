#include "score/calibration.hpp"

#include <cmath>
#include <stdexcept>

namespace wavelab {

BinderPrototype build_prototype(std::vector<Fingerprint> const& binders) {
    BinderPrototype proto;
    if (binders.empty()) return proto;

    proto.spectral_mean.assign(binders.front().spectral.size(), Real{0});
    proto.spectral_freqs = binders.front().spectral_freqs;

    for (auto const& fp : binders) {
        if (fp.spectral.size() != proto.spectral_mean.size()) {
            throw std::invalid_argument(
                "build_prototype: spectral length mismatch across fingerprints");
        }
        for (std::size_t i = 0; i < fp.spectral.size(); ++i) {
            proto.spectral_mean[i] += fp.spectral[i];
        }
        for (auto const& [k, v] : fp.scalars) {
            proto.scalar_mean[k] += v;
        }
    }

    Real const inv = Real{1} / static_cast<Real>(binders.size());
    for (auto& v : proto.spectral_mean) v *= inv;
    for (auto& [k, v] : proto.scalar_mean) v *= inv;
    proto.sample_count = static_cast<Index>(binders.size());
    return proto;
}

Real correlate_with_prototype(Fingerprint const& candidate,
                              BinderPrototype const& proto) {
    if (candidate.spectral.size() != proto.spectral_mean.size()) return Real{0};
    Real dot = Real{0};
    Real nc  = Real{0};
    Real np  = Real{0};
    for (std::size_t i = 0; i < candidate.spectral.size(); ++i) {
        Real const c = candidate.spectral[i];
        Real const p = proto.spectral_mean[i];
        dot += c * p;
        nc  += c * c;
        np  += p * p;
    }
    Real const denom = std::sqrt(nc) * std::sqrt(np);
    if (denom == Real{0}) return Real{0};
    return dot / denom;
}

Real scalar_agreement(Fingerprint const& candidate,
                      BinderPrototype const& proto,
                      Real scale) {
    if (proto.scalar_mean.empty() || candidate.scalars.empty()) return Real{0};

    Real sum_sq = Real{0};
    std::size_t shared = 0;
    for (auto const& [k, vp] : proto.scalar_mean) {
        auto it = candidate.scalars.find(k);
        if (it == candidate.scalars.end()) continue;
        Real const d = it->second - vp;
        sum_sq += d * d;
        ++shared;
    }
    if (shared == 0) return Real{0};
    Real const dist = std::sqrt(sum_sq);
    return std::exp(-dist / scale);
}

} // namespace wavelab
