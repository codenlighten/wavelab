#!/usr/bin/env python3
"""Validate the Vina pipeline by re-docking a known ligand into its
apo receptor and checking heavy-atom RMSD vs the crystal pose.

Phase 8.5 gate:
  * saquinavir (1HXB ROC) re-docks within RMSD 2.5 Å of crystal
  * benzamidine (3PTB BEN) re-docks within RMSD 1.5 Å of crystal
"""

from __future__ import annotations

import argparse
import math
import sys
import tempfile
from pathlib import Path
from typing import List, Tuple

from vina_dock import dock, parse_pdbqt_heavy_atoms, run


def parse_pdb_heavy_atoms(pdb: Path) -> List[Tuple[str, Tuple[float, float, float]]]:
    """Return [(atom_name, (x,y,z))] for non-hydrogen ATOM/HETATM records."""
    out = []
    with open(pdb) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            element = line[76:78].strip()
            if not element:
                # Fall back to atom-name first char
                atom = line[12:16].strip()
                element = atom[0] if atom else "X"
            if element.upper() == "H":
                continue
            name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            out.append((name, (x, y, z)))
    return out


def heavy_atom_rmsd(a: List[Tuple[str, Tuple[float, float, float]]],
                    b: List[Tuple[str, Tuple[float, float, float]]]) -> float:
    """Naive atom-pair RMSD: pairs atoms by ORDER (no permutation search).
    Fine when both molecules came from the same input + same atom
    ordering, which is the redock case here."""
    if len(a) != len(b):
        raise ValueError(f"atom count mismatch: {len(a)} vs {len(b)}")
    if not a:
        return 0.0
    sq = 0.0
    for (_, p), (_, q) in zip(a, b):
        sq += sum((p[i] - q[i]) ** 2 for i in range(3))
    return math.sqrt(sq / len(a))


def best_pose_rmsd(poses: List[Path], crystal: List[Tuple[str, Tuple[float, float, float]]]) -> Tuple[int, float]:
    """Score each pose; return (best_pose_idx, best_rmsd)."""
    best_idx = 0
    best = float("inf")
    for i, p in enumerate(poses):
        atoms = parse_pdbqt_heavy_atoms(p)
        try:
            r = heavy_atom_rmsd(atoms, crystal)
        except ValueError as e:
            print(f"  pose {i}: skip ({e})", file=sys.stderr)
            continue
        print(f"  pose {i}: rmsd = {r:.2f} Å  ({p.name})")
        if r < best:
            best = r
            best_idx = i
    return best_idx, best


def redock_one(name: str, receptor_pdb: Path, crystal_pdb: Path,
               center: Tuple[float, float, float], threshold_A: float,
               out_dir: Path) -> bool:
    print(f"\n=== Redock: {name} ===")
    print(f"  receptor: {receptor_pdb}")
    print(f"  crystal ligand: {crystal_pdb}")
    print(f"  search box center: {center}, threshold {threshold_A} Å")

    out_dir.mkdir(parents=True, exist_ok=True)

    # Dock — pass the PDB directly so ligand_to_pdbqt's alt-loc filter
    # triggers (otherwise obabel emits a dual-ROOT PDBQT that Vina
    # rejects).
    result = dock(receptor_pdb, crystal_pdb, center=center,
                  size_A=(20.0, 20.0, 20.0), num_modes=9,
                  exhaustiveness=8, seed=42,
                  out_dir=out_dir / f"{name}_vina")

    print(f"  vina scores: {[f'{s:+.2f}' for s in result.pose_scores]}")

    # Crystal coordinates (heavy atoms only). The crystal PDB still has
    # alt-loc A and B; filter to A so the count matches Vina's poses
    # (which were docked from the alt-A-only ligand).
    crystal_atoms_full = parse_pdb_heavy_atoms(crystal_pdb)
    # Re-parse with alt-loc filter to match what we docked.
    crystal_atoms: List[Tuple[str, Tuple[float, float, float]]] = []
    with open(crystal_pdb) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            alt = line[16] if len(line) > 16 else " "
            if alt != " " and alt != "A":
                continue
            elem = line[76:78].strip()
            if not elem:
                atom = line[12:16].strip()
                elem = atom[0] if atom else "X"
            if elem.upper() == "H":
                continue
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            crystal_atoms.append((elem, (x, y, z)))
    print(f"  crystal heavy atoms: {len(crystal_atoms)} (alt-A only;"
          f" full file had {len(crystal_atoms_full)})")

    best_idx, best_rmsd = best_pose_rmsd(result.pose_files, crystal_atoms)
    print(f"  BEST pose: #{best_idx} with rmsd {best_rmsd:.2f} Å")
    ok = best_rmsd <= threshold_A
    print(f"  {'PASS' if ok else 'FAIL'}: threshold {threshold_A} Å")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("/root/wavelab/wavelab"))
    ap.add_argument("--out-dir", type=Path, default=Path("/tmp/vina_redock"))
    args = ap.parse_args()

    root = args.root
    cases = [
        # (name, apo, crystal, center, threshold_A)
        ("saquinavir_1HXB",
         root / "data/pdb/apo/1HXB_apo.pdb",
         root / "data/ligands/extracted/1HXB_ROC.pdb",
         (18.257, -0.023, 11.416),
         2.5),
        ("benzamidine_3PTB",
         root / "data/pdb/apo/3PTB_apo.pdb",
         root / "data/ligands/extracted/3PTB_BEN.pdb",
         (-1.759, 14.461, 16.916),
         1.5),
    ]

    import traceback
    pass_count = 0
    for name, apo, crystal, center, thresh in cases:
        try:
            if redock_one(name, apo, crystal, center, thresh, args.out_dir):
                pass_count += 1
        except Exception as e:
            print(f"  ERROR on {name}: {e}", file=sys.stderr)
            traceback.print_exc()

    print()
    print(f"=== redock validation: {pass_count}/{len(cases)} passed ===")
    return 0 if pass_count == len(cases) else 1


if __name__ == "__main__":
    sys.exit(main())
