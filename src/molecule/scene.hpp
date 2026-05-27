#pragma once
//
// MolecularScene<D>: a collection of Atoms living inside a world-space
// axis-aligned box. The simulation grid is independent (built from
// box_min/box_max + a chosen spacing).
//
// All downstream pipeline pieces (synthetic generators, PDB/SDF parsers,
// density/charge/hydro splatting) produce or consume one of these. This
// is the single data-flow throat between molecule I/O and the engine.
//

#include "core/types.hpp"
#include "molecule/atom.hpp"

#include <vector>

namespace wavelab {

template <int D>
struct MolecularScene {
    std::vector<Atom<D>> atoms;
    Vec<D>               box_min{};   // domain lower corner (world coords)
    Vec<D>               box_max{};   // domain upper corner

    Vec<D> extent() const noexcept {
        Vec<D> e{};
        for (std::size_t i = 0; i < static_cast<std::size_t>(D); ++i) {
            e[i] = box_max[i] - box_min[i];
        }
        return e;
    }
};

} // namespace wavelab
