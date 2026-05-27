#pragma once
//
// BoundaryCondition<D>: applied to u^{t+1} after the FDTD stencil update
// and source injection. Three flavors (overview §31):
//
//   * Dirichlet      : u = 0 on ∂Ω (perfect reflector, pressure-release)
//   * Neumann        : ∂u/∂n = 0 (perfect reflector, hard wall)
//   * Mur (1st order): radiation BC, absorbs normally-incident waves
//
// Phase 1 implements the D=1 specializations directly. D=2, D=3 will be
// added in Phase 2 (with PML as the production absorber).
//

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/types.hpp"

#include <memory>
#include <stdexcept>

namespace wavelab {

template <int D>
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    // Modify u_next in place. u_curr is the pre-update field (u^t), needed
    // by Mur for its time-difference term.
    virtual void apply(Field<Real, D>&       u_next,
                       Field<Real, D> const& u_curr,
                       Real                  c,
                       Real                  dt) const = 0;
};

// ============================================================
// 1D specializations
// ============================================================

class Dirichlet1D final : public BoundaryCondition<1> {
public:
    void apply(Field<Real, 1>&       u_next,
               Field<Real, 1> const& /*u_curr*/,
               Real /*c*/, Real /*dt*/) const override {
        Index const n = u_next.grid().shape[0];
        u_next(Index{0})  = Real{0};
        u_next(n - 1)     = Real{0};
    }
};

class Neumann1D final : public BoundaryCondition<1> {
public:
    void apply(Field<Real, 1>&       u_next,
               Field<Real, 1> const& /*u_curr*/,
               Real /*c*/, Real /*dt*/) const override {
        Index const n = u_next.grid().shape[0];
        u_next(Index{0})  = u_next(Index{1});
        u_next(n - 1)     = u_next(n - 2);
    }
};

// 1st-order Mur absorbing boundary:
//
//   u_b^{t+1} = u_{b±1}^{t} + α (u_{b±1}^{t+1} - u_b^{t})
//
// where α = (c·dt - dx) / (c·dt + dx). Exact for normal incidence,
// 5–10% reflection at oblique incidence (an upgrade to PML lands in
// Phase 2 for the production 2D/3D pipeline).
//
class Mur1D final : public BoundaryCondition<1> {
public:
    void apply(Field<Real, 1>&       u_next,
               Field<Real, 1> const& u_curr,
               Real                  c,
               Real                  dt) const override {
        auto const& g = u_next.grid();
        Index const n  = g.shape[0];
        Real const  dx = g.spacing[0];
        Real const  a  = (c * dt - dx) / (c * dt + dx);

        // Left boundary: i = 0, interior neighbor i = 1
        u_next(Index{0}) = u_curr(Index{1}) + a * (u_next(Index{1}) - u_curr(Index{0}));
        // Right boundary: i = n-1, interior neighbor i = n-2
        u_next(n - 1) = u_curr(n - 2) + a * (u_next(n - 2) - u_curr(n - 1));
    }
};

// ============================================================
// 2D specializations
// ============================================================

class Dirichlet2D final : public BoundaryCondition<2> {
public:
    void apply(Field<Real, 2>&       u_next,
               Field<Real, 2> const& /*u_curr*/,
               Real /*c*/, Real /*dt*/) const override {
        auto const& g = u_next.grid();
        Index const nx = g.shape[0];
        Index const ny = g.shape[1];
        for (Index j = 0; j < ny; ++j) {
            u_next(Index{0},       j) = Real{0};
            u_next(nx - 1,         j) = Real{0};
        }
        for (Index i = 0; i < nx; ++i) {
            u_next(i, Index{0})       = Real{0};
            u_next(i, ny - 1)         = Real{0};
        }
    }
};

class Neumann2D final : public BoundaryCondition<2> {
public:
    void apply(Field<Real, 2>&       u_next,
               Field<Real, 2> const& /*u_curr*/,
               Real /*c*/, Real /*dt*/) const override {
        auto const& g = u_next.grid();
        Index const nx = g.shape[0];
        Index const ny = g.shape[1];
        for (Index j = 0; j < ny; ++j) {
            u_next(Index{0},  j) = u_next(Index{1},  j);
            u_next(nx - 1,    j) = u_next(nx - 2,    j);
        }
        for (Index i = 0; i < nx; ++i) {
            u_next(i, Index{0}) = u_next(i, Index{1});
            u_next(i, ny - 1)   = u_next(i, ny - 2);
        }
    }
};

// 1st-order Mur on all four faces. Corner cells are written twice (the
// last face write wins) — acceptable for Mur, irrelevant once PML covers
// the production case.
class Mur2D final : public BoundaryCondition<2> {
public:
    void apply(Field<Real, 2>&       u_next,
               Field<Real, 2> const& u_curr,
               Real                  c,
               Real                  dt) const override {
        auto const& g = u_next.grid();
        Index const nx = g.shape[0];
        Index const ny = g.shape[1];
        Real const  dx = g.spacing[0];
        Real const  dy = g.spacing[1];
        Real const  ax = (c * dt - dx) / (c * dt + dx);
        Real const  ay = (c * dt - dy) / (c * dt + dy);

        // Left (i=0) and right (i=nx-1) faces
        for (Index j = 1; j < ny - 1; ++j) {
            u_next(Index{0}, j) = u_curr(Index{1}, j)
                + ax * (u_next(Index{1}, j) - u_curr(Index{0}, j));
            u_next(nx - 1, j)   = u_curr(nx - 2, j)
                + ax * (u_next(nx - 2, j) - u_curr(nx - 1, j));
        }
        // Top (j=0) and bottom (j=ny-1) faces
        for (Index i = 1; i < nx - 1; ++i) {
            u_next(i, Index{0}) = u_curr(i, Index{1})
                + ay * (u_next(i, Index{1}) - u_curr(i, Index{0}));
            u_next(i, ny - 1)   = u_curr(i, ny - 2)
                + ay * (u_next(i, ny - 2) - u_curr(i, ny - 1));
        }
        // Corners: average of two face values (cheap heuristic).
        u_next(Index{0},  Index{0})  = Real{0.5} * (u_next(Index{1}, Index{0}) + u_next(Index{0}, Index{1}));
        u_next(nx - 1,    Index{0})  = Real{0.5} * (u_next(nx - 2, Index{0}) + u_next(nx - 1, Index{1}));
        u_next(Index{0},  ny - 1)    = Real{0.5} * (u_next(Index{1}, ny - 1) + u_next(Index{0}, ny - 2));
        u_next(nx - 1,    ny - 1)    = Real{0.5} * (u_next(nx - 2, ny - 1) + u_next(nx - 1, ny - 2));
    }
};

// ============================================================
// 3D specializations
// ============================================================

class Dirichlet3D final : public BoundaryCondition<3> {
public:
    void apply(Field<Real, 3>&       u_next,
               Field<Real, 3> const& /*u_curr*/,
               Real /*c*/, Real /*dt*/) const override {
        auto const& g = u_next.grid();
        Index const nx = g.shape[0];
        Index const ny = g.shape[1];
        Index const nz = g.shape[2];
        for (Index j = 0; j < ny; ++j)
            for (Index k = 0; k < nz; ++k) {
                u_next(Index{0},   j, k) = Real{0};
                u_next(nx - 1,     j, k) = Real{0};
            }
        for (Index i = 0; i < nx; ++i)
            for (Index k = 0; k < nz; ++k) {
                u_next(i, Index{0},   k) = Real{0};
                u_next(i, ny - 1,     k) = Real{0};
            }
        for (Index i = 0; i < nx; ++i)
            for (Index j = 0; j < ny; ++j) {
                u_next(i, j, Index{0}) = Real{0};
                u_next(i, j, nz - 1)   = Real{0};
            }
    }
};

class Neumann3D final : public BoundaryCondition<3> {
public:
    void apply(Field<Real, 3>&       u_next,
               Field<Real, 3> const& /*u_curr*/,
               Real /*c*/, Real /*dt*/) const override {
        auto const& g = u_next.grid();
        Index const nx = g.shape[0];
        Index const ny = g.shape[1];
        Index const nz = g.shape[2];
        for (Index j = 0; j < ny; ++j)
            for (Index k = 0; k < nz; ++k) {
                u_next(Index{0},   j, k) = u_next(Index{1}, j, k);
                u_next(nx - 1,     j, k) = u_next(nx - 2,   j, k);
            }
        for (Index i = 0; i < nx; ++i)
            for (Index k = 0; k < nz; ++k) {
                u_next(i, Index{0},   k) = u_next(i, Index{1}, k);
                u_next(i, ny - 1,     k) = u_next(i, ny - 2,   k);
            }
        for (Index i = 0; i < nx; ++i)
            for (Index j = 0; j < ny; ++j) {
                u_next(i, j, Index{0}) = u_next(i, j, Index{1});
                u_next(i, j, nz - 1)   = u_next(i, j, nz - 2);
            }
    }
};

// 1st-order Mur on all six faces. Edges/corners written multiple times
// (last write wins); for production use PML covers most absorption.
class Mur3D final : public BoundaryCondition<3> {
public:
    void apply(Field<Real, 3>&       u_next,
               Field<Real, 3> const& u_curr,
               Real                  c,
               Real                  dt) const override {
        auto const& g = u_next.grid();
        Index const nx = g.shape[0];
        Index const ny = g.shape[1];
        Index const nz = g.shape[2];
        Real const  dx = g.spacing[0];
        Real const  dy = g.spacing[1];
        Real const  dz = g.spacing[2];
        Real const  ax = (c * dt - dx) / (c * dt + dx);
        Real const  ay = (c * dt - dy) / (c * dt + dy);
        Real const  az = (c * dt - dz) / (c * dt + dz);

        for (Index j = 1; j < ny - 1; ++j)
            for (Index k = 1; k < nz - 1; ++k) {
                u_next(Index{0}, j, k) = u_curr(Index{1}, j, k)
                    + ax * (u_next(Index{1}, j, k) - u_curr(Index{0}, j, k));
                u_next(nx - 1, j, k) = u_curr(nx - 2, j, k)
                    + ax * (u_next(nx - 2, j, k) - u_curr(nx - 1, j, k));
            }
        for (Index i = 1; i < nx - 1; ++i)
            for (Index k = 1; k < nz - 1; ++k) {
                u_next(i, Index{0}, k) = u_curr(i, Index{1}, k)
                    + ay * (u_next(i, Index{1}, k) - u_curr(i, Index{0}, k));
                u_next(i, ny - 1, k) = u_curr(i, ny - 2, k)
                    + ay * (u_next(i, ny - 2, k) - u_curr(i, ny - 1, k));
            }
        for (Index i = 1; i < nx - 1; ++i)
            for (Index j = 1; j < ny - 1; ++j) {
                u_next(i, j, Index{0}) = u_curr(i, j, Index{1})
                    + az * (u_next(i, j, Index{1}) - u_curr(i, j, Index{0}));
                u_next(i, j, nz - 1) = u_curr(i, j, nz - 2)
                    + az * (u_next(i, j, nz - 2) - u_curr(i, j, nz - 1));
            }
    }
};

// ============================================================
// Dim-generic factories
// ============================================================
// Convenience helpers so callers don't branch on D when constructing
// a boundary. Returns a shared_ptr<BoundaryCondition<D>> ready to pass
// to Stepper<D>::set_boundary.

template <int D>
inline std::shared_ptr<BoundaryCondition<D>> make_dirichlet() {
    if constexpr (D == 1) return std::make_shared<Dirichlet1D>();
    else if constexpr (D == 2) return std::make_shared<Dirichlet2D>();
    else if constexpr (D == 3) return std::make_shared<Dirichlet3D>();
}

template <int D>
inline std::shared_ptr<BoundaryCondition<D>> make_neumann() {
    if constexpr (D == 1) return std::make_shared<Neumann1D>();
    else if constexpr (D == 2) return std::make_shared<Neumann2D>();
    else if constexpr (D == 3) return std::make_shared<Neumann3D>();
}

template <int D>
inline std::shared_ptr<BoundaryCondition<D>> make_mur() {
    if constexpr (D == 1) return std::make_shared<Mur1D>();
    else if constexpr (D == 2) return std::make_shared<Mur2D>();
    else if constexpr (D == 3) return std::make_shared<Mur3D>();
}

} // namespace wavelab
