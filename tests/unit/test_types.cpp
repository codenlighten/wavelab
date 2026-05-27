#include <doctest/doctest.h>

#include "core/types.hpp"

using namespace wavelab;

TEST_CASE("Vec and IVec basic construction") {
    Vec<3>  v = make_vec<3>(Real{1.5});
    IVec<3> i = make_ivec<3>(Index{4});

    CHECK(v[0] == doctest::Approx(1.5));
    CHECK(v[2] == doctest::Approx(1.5));
    CHECK(i[0] == 4);
    CHECK(i[2] == 4);
}

TEST_CASE("dot and product") {
    Vec<3>  a{1, 2, 3};
    Vec<3>  b{4, 5, 6};
    IVec<3> shape{2, 3, 4};

    CHECK(dot<3>(a, b) == doctest::Approx(1 * 4 + 2 * 5 + 3 * 6));
    CHECK(product<3>(shape) == 24);
}
