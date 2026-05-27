#include <doctest/doctest.h>

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"
#include "score/entropy.hpp"

#include <cmath>

using namespace wavelab;
using namespace wavelab::literals;

namespace {

TEST_CASE("entropy of uniform energy field ≈ log(M)") {
    auto g = Grid<2>::uniform(IVec<2>{16, 16}, 1.0_r);
    Field<Real, 2> e(g, 1.0_r);    // uniform energy
    Real const H     = energy_entropy(e);
    Real const log_M = std::log(static_cast<Real>(e.size()));
    CHECK(H == doctest::Approx(static_cast<double>(log_M)).epsilon(1e-5));
    CHECK(focus_score(e) == doctest::Approx(0.0).epsilon(1e-5));
}

TEST_CASE("entropy of a delta is 0; focus is 1") {
    auto g = Grid<2>::uniform(IVec<2>{8, 8}, 1.0_r);
    Field<Real, 2> e(g, 0.0_r);
    e(3, 4) = 10.0_r;
    CHECK(energy_entropy(e) == doctest::Approx(0.0).epsilon(1e-5));
    CHECK(focus_score(e) == doctest::Approx(1.0).epsilon(1e-5));
}

TEST_CASE("hotspot_concentration: full region == 1, empty == 0") {
    auto g = Grid<2>::uniform(IVec<2>{8, 8}, 1.0_r);
    Field<Real, 2> e(g, 1.0_r);

    Field<Real, 2> mask_full(g, 1.0_r);
    CHECK(hotspot_concentration(e, mask_full) == doctest::Approx(1.0));

    Field<Real, 2> mask_none(g, 0.0_r);
    CHECK(hotspot_concentration(e, mask_none) == doctest::Approx(0.0));

    // Half-area mask
    Field<Real, 2> mask_half(g, 0.0_r);
    for (Index i = 0; i < 4; ++i)
        for (Index j = 0; j < 8; ++j)
            mask_half(i, j) = 1.0_r;
    CHECK(hotspot_concentration(e, mask_half) == doctest::Approx(0.5));
}

TEST_CASE("hotspot_concentration_region: rectangular probe") {
    auto g = Grid<2>::uniform(IVec<2>{4, 4}, 1.0_r);
    Field<Real, 2> e(g, 1.0_r);

    // Top-left 2x2 quadrant = 4 cells out of 16 = 0.25
    Real const r = hotspot_concentration_region(e, IVec<2>{0, 0}, IVec<2>{1, 1});
    CHECK(r == doctest::Approx(0.25));
}

} // namespace
