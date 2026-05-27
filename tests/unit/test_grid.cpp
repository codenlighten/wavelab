#include <doctest/doctest.h>

#include "core/grid.hpp"

using namespace wavelab;
using namespace wavelab::literals;

TEST_CASE("Grid<2> shape and indexing") {
    auto g = Grid<2>::uniform(IVec<2>{4, 5}, 0.1_r);

    CHECK(g.num_cells() == 4 * 5);
    CHECK(g.min_spacing() == doctest::Approx(0.1));

    // Row-major: (i, j) -> i*ny + j
    CHECK(g.linear_index({0, 0}) == 0);
    CHECK(g.linear_index({0, 4}) == 4);
    CHECK(g.linear_index({1, 0}) == 5);
    CHECK(g.linear_index({3, 4}) == 19);

    CHECK(g.in_bounds({0, 0}));
    CHECK(g.in_bounds({3, 4}));
    CHECK_FALSE(g.in_bounds({4, 0}));
    CHECK_FALSE(g.in_bounds({-1, 0}));
}

TEST_CASE("Grid<3> indexing and cell_center") {
    auto g = Grid<3>(IVec<3>{2, 3, 4},
                     Vec<3>{0.5_r, 0.25_r, 0.1_r},
                     Vec<3>{0.0_r, 0.0_r, 0.0_r});

    CHECK(g.num_cells() == 2 * 3 * 4);
    CHECK(g.linear_index({1, 2, 3}) == (1 * 3 + 2) * 4 + 3);
    CHECK(g.min_spacing() == doctest::Approx(0.1));

    auto c = g.cell_center({0, 0, 0});
    CHECK(c[0] == doctest::Approx(0.25));   // 0.5 * 0.5
    CHECK(c[1] == doctest::Approx(0.125));
    CHECK(c[2] == doctest::Approx(0.05));

    auto ext = g.extent();
    CHECK(ext[0] == doctest::Approx(1.0));  // 2 * 0.5
    CHECK(ext[1] == doctest::Approx(0.75));
    CHECK(ext[2] == doctest::Approx(0.4));
}

TEST_CASE("Grid<1> degenerate case") {
    auto g = Grid<1>::uniform(IVec<1>{100}, 0.01_r);
    CHECK(g.num_cells() == 100);
    CHECK(g.linear_index({37}) == 37);
}
