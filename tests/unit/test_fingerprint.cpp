#include <doctest/doctest.h>

#include "core/types.hpp"
#include "io/fingerprint.hpp"

#include <cmath>

using namespace wavelab;
using namespace wavelab::literals;

TEST_CASE("Fingerprint round-trip JSON") {
    Fingerprint fp;
    fp.scene_name        = "test scene \"with quotes\"";
    fp.scalars["S_wave"] = 0.875_r;
    fp.scalars["R_E"]    = 1.234_r;
    fp.scalars["H_E"]    = 3.14_r;
    fp.spectral          = {0.1_r, 0.2_r, 0.3_r, 0.4_r};
    fp.spectral_freqs    = {1.0_r, 2.0_r, 4.0_r, 8.0_r};
    fp.meta["engine"]    = "wavelab 0.0.1";
    fp.meta["date"]      = "2026-05-27";

    auto const text = fingerprint_to_json(fp);
    auto const rt   = fingerprint_from_json(text);

    CHECK(rt.scene_name == fp.scene_name);
    CHECK(rt.scalars.size() == fp.scalars.size());
    for (auto const& [k, v] : fp.scalars) {
        REQUIRE(rt.scalars.count(k) == 1u);
        CHECK(rt.scalars.at(k) == doctest::Approx(static_cast<double>(v)));
    }
    REQUIRE(rt.spectral.size() == fp.spectral.size());
    for (std::size_t i = 0; i < fp.spectral.size(); ++i) {
        CHECK(rt.spectral[i] == doctest::Approx(static_cast<double>(fp.spectral[i])));
    }
    CHECK(rt.meta.at("engine") == "wavelab 0.0.1");
    CHECK(rt.meta.at("date")   == "2026-05-27");
}

TEST_CASE("Fingerprint deterministic serialization (sorted keys)") {
    Fingerprint a;
    a.scalars["c"] = 3.0_r;
    a.scalars["a"] = 1.0_r;
    a.scalars["b"] = 2.0_r;
    a.meta["zz"]   = "z";
    a.meta["aa"]   = "a";

    Fingerprint b;
    b.scalars["b"] = 2.0_r;
    b.scalars["c"] = 3.0_r;
    b.scalars["a"] = 1.0_r;
    b.meta["aa"]   = "a";
    b.meta["zz"]   = "z";

    // Insertion order differs, but serialized output is the same.
    CHECK(fingerprint_to_json(a) == fingerprint_to_json(b));
}

TEST_CASE("Fingerprint handles NaN as null") {
    Fingerprint fp;
    fp.scalars["nan_val"] = static_cast<Real>(std::nan(""));
    auto const text = fingerprint_to_json(fp);
    CHECK(text.find("null") != std::string::npos);
    auto const rt = fingerprint_from_json(text);
    CHECK(std::isnan(static_cast<double>(rt.scalars.at("nan_val"))));
}
