#include <doctest/doctest.h>

#include "core/field.hpp"

#include <stdexcept>

using namespace wavelab;
using namespace wavelab::literals;

TEST_CASE("Field<float,1> construct, write, read") {
    auto g = Grid<1>::uniform(IVec<1>{8}, 1.0_r);
    Field<float, 1> f(g, 0.0f);

    CHECK(f.size() == 8);
    CHECK(f.data() != nullptr);

    for (Index i = 0; i < 8; ++i) f(i) = static_cast<float>(i * i);

    for (Index i = 0; i < 8; ++i) {
        CHECK(f(i) == doctest::Approx(static_cast<float>(i * i)));
    }
}

TEST_CASE("Field<float,2> variadic and IVec indexing agree") {
    auto g = Grid<2>::uniform(IVec<2>{4, 5}, 0.1_r);
    Field<float, 2> f(g, 0.0f);

    f(2, 3) = 7.5f;
    CHECK(f({2, 3}) == doctest::Approx(7.5f));
    CHECK(f.data()[2 * 5 + 3] == doctest::Approx(7.5f));
}

TEST_CASE("Field<double,3> fill and zero") {
    auto g = Grid<3>::uniform(IVec<3>{3, 3, 3}, 1.0_r);
    Field<double, 3> f(g, 0.0);

    f.fill(2.0);
    for (Index i = 0; i < f.size(); ++i) {
        CHECK(f.data()[i] == doctest::Approx(2.0));
    }

    f.zero();
    for (Index i = 0; i < f.size(); ++i) {
        CHECK(f.data()[i] == doctest::Approx(0.0));
    }
}

TEST_CASE("Field::at throws on out-of-bounds") {
    auto g = Grid<2>::uniform(IVec<2>{2, 2}, 1.0_r);
    Field<float, 2> f(g);

    CHECK_THROWS_AS(f.at({2, 0}), std::out_of_range);
    CHECK_THROWS_AS(f.at({0, -1}), std::out_of_range);
    CHECK_NOTHROW(f.at({1, 1}));
}

TEST_CASE("Field swap") {
    auto g = Grid<1>::uniform(IVec<1>{4}, 1.0_r);
    Field<float, 1> a(g, 1.0f);
    Field<float, 1> b(g, 2.0f);

    a.swap(b);

    CHECK(a(0) == doctest::Approx(2.0f));
    CHECK(b(0) == doctest::Approx(1.0f));
}
