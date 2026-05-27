"""End-to-end smoke test of the wavelab pybind11 module.

The 5-line Phase 6 gate: a Python script that runs the pipeline on a
real PDB and reads back scalar + spectral fingerprint data.
"""

import os
import sys


def main() -> int:
    pdb_path = os.environ.get(
        "WAVELAB_TEST_PDB",
        os.path.join(os.path.dirname(__file__), "..", "..", "data", "pdb", "1ubq.pdb"),
    )
    if not os.path.exists(pdb_path):
        print(f"SKIP: {pdb_path} not found", file=sys.stderr)
        return 0

    import wavelab

    # 1. Run a synthetic scene
    fp_dumb = wavelab.run_synthetic("dumbbell", nx=60, ny=60, steps=120)
    # 2. Run on the PDB
    fp_pdb  = wavelab.run_pdb(pdb_path, nx=60, ny=60, steps=120)
    # 3. Build a prototype from N binders
    proto   = wavelab.build_prototype([fp_dumb, fp_dumb])
    # 4. Correlate a candidate against the prototype
    r       = wavelab.correlate(fp_pdb, proto)

    print(f"dumbbell.total_energy = {fp_dumb['scalars']['total_energy']:.4f}")
    print(f"pdb.atoms             = {fp_pdb['meta']['atoms']}")
    print(f"pdb.total_energy      = {fp_pdb['scalars']['total_energy']:.4f}")
    print(f"pdb.entropy           = {fp_pdb['scalars']['entropy']:.4f}")
    print(f"pdb.spectral_len      = {len(fp_pdb['spectral'])}")
    print(f"correlate(pdb, dumbbell-proto) = {r:.4f}")

    # Sanity assertions for CI.
    assert fp_dumb["scalars"]["total_energy"] > 0.0
    assert fp_pdb["scalars"]["entropy"] > 0.0
    assert len(fp_pdb["spectral"]) == 16
    assert -1.0 <= r <= 1.0
    return 0


if __name__ == "__main__":
    sys.exit(main())
