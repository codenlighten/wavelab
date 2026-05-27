#include <doctest/doctest.h>

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"
#include "score/similarity.hpp"

#include <cmath>
#include <numbers>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

TEST_CASE("wave_similarity: identical fields score 1 exactly") {
    auto g = Grid<2>::uniform(IVec<2>{16, 16}, 0.1_r);
    Field<float, 2> a(g, 0.0f);
    for (Index i = 0; i < 16; ++i)
        for (Index j = 0; j < 16; ++j)
            a(i, j) = std::sin(static_cast<float>(i) * 0.3f
                             + static_cast<float>(j) * 0.5f);

    CHECK(wave_difference<2>(a, a) == doctest::Approx(0.0));
    CHECK(wave_similarity<2>(a, a, 1.0_r) == doctest::Approx(1.0));
    CHECK(shape_complementarity<2>(a, a) == doctest::Approx(1.0));
}

TEST_CASE("shape_complementarity: orthogonal patterns score ~0") {
    auto g = Grid<2>::uniform(IVec<2>{32, 32}, 0.1_r);
    Field<float, 2> a(g, 0.0f);
    Field<float, 2> b(g, 0.0f);

    // a: sin(kx) — pure x-direction wave
    // b: sin(ky) — pure y-direction wave (orthogonal in inner product)
    for (Index i = 0; i < 32; ++i) {
        for (Index j = 0; j < 32; ++j) {
            a(i, j) = std::sin(2.0f * static_cast<float>(std::numbers::pi)
                             * static_cast<float>(i) / 32.0f);
            b(i, j) = std::sin(2.0f * static_cast<float>(std::numbers::pi)
                             * static_cast<float>(j) / 32.0f);
        }
    }

    Real const s = shape_complementarity<2>(a, b);
    CHECK(std::abs(s) < 0.05_r);
}

TEST_CASE("wave_similarity decays with field difference") {
    auto g = Grid<2>::uniform(IVec<2>{8, 8}, 1.0_r);
    Field<float, 2> a(g, 1.0f);
    Field<float, 2> b(g, 1.0f);

    // Identical -> 1
    CHECK(wave_similarity<2>(a, b) == doctest::Approx(1.0));

    // Perturb b
    b(4, 4) = 2.0f;
    Real const s1 = wave_similarity<2>(a, b, 1.0_r);
    CHECK(s1 < 1.0_r);
    CHECK(s1 > 0.0_r);

    // Bigger perturbation -> lower similarity
    b(4, 4) = 5.0f;
    Real const s2 = wave_similarity<2>(a, b, 1.0_r);
    CHECK(s2 < s1);
}

} // namespace
