// wavecli — headless batch runner. Phase 0: smoke test that the core
// library links and reports its configuration. Real pipeline wires in
// starting in Phase 4 (fingerprint emission).

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"

#include <cstdio>
#include <string>

int main(int argc, char** argv) {
    std::puts("wavelab CLI — Phase 0 skeleton");

#if WAVELAB_DOUBLE_PRECISION
    std::puts("  precision : double");
#else
    std::puts("  precision : float");
#endif
#if WAVELAB_HAVE_OPENMP
    std::puts("  openmp    : yes");
#else
    std::puts("  openmp    : no");
#endif

    using namespace wavelab::literals;
    auto grid = wavelab::Grid<3>::uniform(wavelab::IVec<3>{8, 8, 8}, 1.0_r);
    wavelab::Field<wavelab::Real, 3> u(grid, 0.0_r);
    u(4, 4, 4) = 1.0_r;

    std::printf("  grid      : %lld x %lld x %lld (%lld cells)\n",
                static_cast<long long>(grid.shape[0]),
                static_cast<long long>(grid.shape[1]),
                static_cast<long long>(grid.shape[2]),
                static_cast<long long>(grid.num_cells()));
    std::printf("  u(4,4,4)  : %g\n", static_cast<double>(u(4, 4, 4)));

    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h") {
            std::puts("\nusage: wavecli [scene.toml] [-o out.fp.json]");
            std::puts("(scene + fingerprint I/O are wired in Phase 4)");
            return 0;
        }
    }
    return 0;
}
