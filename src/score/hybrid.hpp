#pragma once
//
// Hybrid docking score combiner (overview §27, §35).
//
//   S_total = Σ w_i · S_i        weighted sum of individual scores
//   ΔG̅     = ΔG_classical + λ_w · ΔG_wave    binding-energy form
//
// The combiner is dumb on purpose — it just holds named weights and a
// dictionary of named scores, returns the weighted sum. The intelligence
// lives upstream (which scores you compute) and downstream (calibration
// of the weights against ground-truth binding data, §28).
//

#include "core/types.hpp"

#include <string>
#include <unordered_map>

namespace wavelab {

class HybridScore {
public:
    void set_weight(std::string const& term, Real w) { weights_[term] = w; }
    void set_score(std::string const& term, Real s)  { scores_[term]  = s; }

    Real weight(std::string const& term) const {
        auto it = weights_.find(term);
        return (it == weights_.end()) ? Real{0} : it->second;
    }
    Real score(std::string const& term) const {
        auto it = scores_.find(term);
        return (it == scores_.end()) ? Real{0} : it->second;
    }

    // Combine: Σ_i w_i · S_i over all scored terms. Terms with no
    // explicit weight contribute zero (defensive default).
    Real total() const {
        Real s = Real{0};
        for (auto const& [term, sval] : scores_) {
            auto it = weights_.find(term);
            if (it != weights_.end()) s += it->second * sval;
        }
        return s;
    }

    std::unordered_map<std::string, Real> const& weights() const { return weights_; }
    std::unordered_map<std::string, Real> const& scores()  const { return scores_;  }

private:
    std::unordered_map<std::string, Real> weights_;
    std::unordered_map<std::string, Real> scores_;
};

} // namespace wavelab
