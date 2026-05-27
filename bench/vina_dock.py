#!/usr/bin/env python3
"""AutoDock Vina wrapper for wavelab.

Pipeline per call:
  1. Convert receptor PDB → PDBQT via obabel (rigid receptor).
  2. Convert ligand SDF / PDB → PDBQT via obabel (with Gasteiger
     partial charges + added hydrogens).
  3. Run vina with a search box centered at `center`, side length
     `size_A` (Å), keeping the top `num_modes` poses.
  4. Split the multi-model PDBQT output into individual SDFs.

This script is deliberately thin — it shells out to `vina` and
`obabel` rather than using their Python APIs. Easier to debug, no
deps beyond what `apt install autodock-vina openbabel` provides.

Caveats handled here:
  * receptor PDB stripped of HETATM (waters/ligands) before obabel
    conversion (obabel will silently happily PDBQT-ify a HETATM block
    as part of the receptor, which is wrong).
  * obabel `-xr` (rigid receptor) avoids the rotatable-bond detection
    obabel does by default, which can blow up for protein backbones.
  * deterministic seed so re-runs reproduce.

Pose-RMSD validation lives in `vina_redock_validate.py`.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple


@dataclass
class DockResult:
    pose_files: List[Path]   # one PDBQT per pose, ordered best→worst
    pose_scores: List[float]  # Vina affinity (kcal/mol) per pose
    docked_pdbqt: Path       # raw multi-model output


def parse_pdbqt_heavy_atoms(pdbqt: Path) -> List[Tuple[str, Tuple[float, float, float]]]:
    """Extract (element, (x,y,z)) for non-H atoms from a single-pose
    PDBQT (one ROOT block). PDBQT atom records use the same columns
    1-54 as PDB."""
    out: List[Tuple[str, Tuple[float, float, float]]] = []
    with open(pdbqt) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            # Vina PDBQT element/type in columns 77-78 (e.g., "C ", "OA", "HD")
            elem = line[76:78].strip() if len(line) >= 78 else ""
            # Skip explicit hydrogens (Vina includes polar Hs as HD/HS).
            if elem in ("H", "HD", "HS"):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            out.append((elem or "?", (x, y, z)))
    return out


def run(cmd: List[str], **kw) -> subprocess.CompletedProcess:
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(
            f"command failed: {' '.join(cmd)}\nstdout:\n{r.stdout}\nstderr:\n{r.stderr}"
        )
    return r


def receptor_to_pdbqt(pdb: Path, pdbqt_out: Path) -> None:
    # Strip HETATM (waters/ligands/ions) before conversion.
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp:
        with open(pdb) as src:
            for line in src:
                if line.startswith(("ATOM", "TER", "HEADER", "REMARK")):
                    tmp.write(line)
            tmp.write("END\n")
        clean_pdb = tmp.name
    try:
        # -xr = rigid receptor (no rotatable bonds detected)
        # --partialcharge gasteiger = explicit GA-charge assignment
        run(["obabel", clean_pdb, "-O", str(pdbqt_out), "-xr",
             "--partialcharge", "gasteiger"])
    finally:
        os.unlink(clean_pdb)


def ligand_to_pdbqt(lig_in: Path, pdbqt_out: Path,
                    alt_loc_keep: str = "A") -> None:
    # For PDB inputs, pre-filter alt-loc records (obabel treats A/B
    # conformers as separate molecules and emits two ROOT trees, which
    # Vina rejects as malformed ligand PDBQT).
    ext = lig_in.suffix.lower().lstrip(".")
    actual_in: Path = lig_in
    cleanup: Path | None = None
    if ext == "pdb":
        with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp:
            with open(lig_in) as src:
                for line in src:
                    if line.startswith(("ATOM", "HETATM")):
                        alt = line[16] if len(line) > 16 else " "
                        if alt != " " and alt != alt_loc_keep:
                            continue
                    tmp.write(line)
            actual_in = Path(tmp.name)
            cleanup = actual_in
    in_fmt = ext if ext else "sdf"
    try:
        run(["obabel", f"-i{in_fmt}", str(actual_in), "-O", str(pdbqt_out),
             "-h", "--partialcharge", "gasteiger"])
    finally:
        if cleanup is not None:
            try: cleanup.unlink()
            except FileNotFoundError: pass


def split_pdbqt_models(pdbqt: Path, out_dir: Path, stem: str) -> List[Path]:
    """Hand-split a multi-MODEL PDBQT into per-model PDBQT files. Each
    file has a single ROOT/ENDROOT block from one Vina pose."""
    out_dir.mkdir(parents=True, exist_ok=True)
    poses: List[Path] = []
    with open(pdbqt) as f:
        current_lines: List[str] = []
        in_model = False
        for line in f:
            if line.startswith("MODEL"):
                current_lines = []
                in_model = True
            elif line.startswith("ENDMDL"):
                if in_model:
                    idx = len(poses)
                    out = out_dir / f"{stem}_pose_{idx}.pdbqt"
                    with open(out, "w") as g:
                        g.writelines(current_lines)
                    poses.append(out)
                in_model = False
                current_lines = []
            elif in_model:
                current_lines.append(line)
    if not poses:
        # No MODEL/ENDMDL — single-pose file. Copy whole thing.
        out = out_dir / f"{stem}_pose_0.pdbqt"
        with open(pdbqt) as src, open(out, "w") as dst:
            dst.write(src.read())
        poses.append(out)
    return poses


def parse_vina_scores(stdout: str) -> List[float]:
    """Pull per-mode affinity (kcal/mol) from vina's text output."""
    scores: List[float] = []
    in_table = False
    for line in stdout.splitlines():
        if line.strip().startswith("mode |"):
            in_table = True
            continue
        if not in_table:
            continue
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                scores.append(float(parts[1]))
            except ValueError:
                pass
    return scores


def dock(receptor_pdb: Path,
         ligand_in: Path,
         center: Tuple[float, float, float],
         size_A: Tuple[float, float, float] = (20.0, 20.0, 20.0),
         num_modes: int = 9,
         exhaustiveness: int = 8,
         seed: int = 42,
         out_dir: Path = Path(".")) -> DockResult:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    stem = ligand_in.stem
    receptor_qt = out_dir / "receptor.pdbqt"
    ligand_qt = out_dir / f"{stem}.pdbqt"
    docked_qt = out_dir / f"{stem}_docked.pdbqt"

    receptor_to_pdbqt(receptor_pdb, receptor_qt)
    ligand_to_pdbqt(ligand_in, ligand_qt)

    cx, cy, cz = center
    sx, sy, sz = size_A
    vina_cmd = [
        "vina",
        "--receptor", str(receptor_qt),
        "--ligand", str(ligand_qt),
        "--out", str(docked_qt),
        "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
        "--size_x", f"{sx:.1f}", "--size_y", f"{sy:.1f}", "--size_z", f"{sz:.1f}",
        "--num_modes", str(num_modes),
        "--exhaustiveness", str(exhaustiveness),
        "--seed", str(seed),
    ]
    r = run(vina_cmd)
    scores = parse_vina_scores(r.stdout)

    pose_files = split_pdbqt_models(docked_qt, out_dir / f"{stem}_poses", stem)
    return DockResult(pose_files=pose_files,
                      pose_scores=scores,
                      docked_pdbqt=docked_qt)


def main() -> int:
    ap = argparse.ArgumentParser(description="Dock a ligand into a receptor with AutoDock Vina.")
    ap.add_argument("--receptor", type=Path, required=True, help="apo receptor PDB")
    ap.add_argument("--ligand", type=Path, required=True, help="ligand SDF or PDB")
    ap.add_argument("--center", type=str, required=True, help="X,Y,Z of search box center (Å)")
    ap.add_argument("--size", type=str, default="20,20,20", help="X,Y,Z of search box size (Å)")
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", type=Path, default=Path("vina_out"))
    args = ap.parse_args()

    cx, cy, cz = map(float, args.center.split(","))
    sx, sy, sz = map(float, args.size.split(","))
    result = dock(args.receptor, args.ligand,
                  center=(cx, cy, cz), size_A=(sx, sy, sz),
                  num_modes=args.num_modes,
                  exhaustiveness=args.exhaustiveness,
                  seed=args.seed,
                  out_dir=args.out_dir)
    print(f"Docked {args.ligand.name} into {args.receptor.name}")
    print(f"  poses: {len(result.pose_files)}")
    for i, (pf, score) in enumerate(zip(result.pose_files, result.pose_scores)):
        print(f"  pose {i}: score={score:+.2f} kcal/mol  ->  {pf.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
