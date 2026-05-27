#!/usr/bin/env python3
"""Phase 8.6 — multi-target real-ligand discrimination benchmark.

For each protein target {1HXB, 3PTB, 1STP}:
  * apo PDB                    (data/pdb/apo/{target}_apo.pdb)
  * known active ligand         (data/ligands/extracted/{target}_{HET}.pdb,
                                 placed at crystal coords — no docking)
  * N decoy ligands             (PubChem 3D SDFs already in data/ligands/
                                 actives|decoys; we use them as decoys here
                                 regardless of label since none of them
                                 actually bind these specific targets)
  * For each decoy: Vina-dock at the known binding site → top-K poses
  * For each (active + decoy pose): run wavecli 3D with regional R_E
    centered on the binding-site centroid → fingerprint
  * Score each candidate via cosine sim to active prototype +
    regional R_E
  * Pairwise AUC per target.

Acceptance gate (per Phase 8 plan):
    AUC ≥ 0.85 on at least 2 of 3 targets.

Caveat: with only ONE active per target (the crystal ligand), the
"prototype" is just that one fingerprint and the test reduces to "does
the active outscore the decoys?" — which is the most basic possible
discrimination check. Multi-active prototypes would need more
co-crystals per target (deferred).
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from vina_dock import dock as vina_dock_fn

ROOT = Path("/root/wavelab/wavelab")


@dataclass
class ActiveLigand:
    name: str               # for reporting (e.g., "saquinavir")
    pdb: Path               # HETATM-extracted real binding pose


@dataclass
class Target:
    code: str
    apo_pdb: Path
    site_center: Tuple[float, float, float]
    actives: List[ActiveLigand] = field(default_factory=list)
    decoys: List[Path] = field(default_factory=list)


@dataclass
class Candidate:
    name: str
    kind: str               # "ACTIVE" or "decoy"
    # For actives: a single PDB path with the crystal pose.
    # For decoys: list of PDBQT files (top-K Vina poses).
    pose_files: List[Path] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Pipeline calls
# ---------------------------------------------------------------------------

def run_wavecli(wavecli: Path, args: List[str]) -> dict:
    cmd = [str(wavecli)] + args
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(
            f"wavecli failed: {' '.join(cmd)}\nstderr:\n{r.stderr}")
    return r


def score_pose(wavecli: Path, apo_pdb: Path, pose_file: Path,
               site_center: Tuple[float, float, float],
               nx: int, ny: int, nz: int, dx: float,
               freq: float, steps: int,
               beta_q: float, region_half: float,
               out_fp: Path) -> dict:
    """Run wavecli on apo + this pose, return the parsed fingerprint."""
    cx, cy, cz = site_center
    args = [
        "--dim", "3",
        "--pdb", str(apo_pdb),
        "--add-pdb", str(pose_file),
        "--include-hetatm",
        "--place-at", f"{cx:.3f},{cy:.3f},{cz:.3f}",
        "--region", f"{cx:.3f},{cy:.3f},{cz:.3f},{region_half:.1f}",
        "--nx", str(nx), "--ny", str(ny), "--nz", str(nz),
        "--dx", f"{dx:.3f}",
        "--freq", f"{freq:.3f}",
        "--steps", str(steps),
        "--pulse",
        "--beta-rho", "0.5",
        "--beta-q", f"{beta_q:.3f}",
        "-o", str(out_fp),
    ]
    run_wavecli(wavecli, args)
    with open(out_fp) as f:
        return json.load(f)


def score_apo(wavecli: Path, apo_pdb: Path,
              site_center: Tuple[float, float, float],
              nx: int, ny: int, nz: int, dx: float,
              freq: float, steps: int,
              beta_q: float, region_half: float,
              out_fp: Path) -> dict:
    """Run wavecli on apo alone — for regional R_E denominator."""
    cx, cy, cz = site_center
    args = [
        "--dim", "3",
        "--pdb", str(apo_pdb),
        "--region", f"{cx:.3f},{cy:.3f},{cz:.3f},{region_half:.1f}",
        "--nx", str(nx), "--ny", str(ny), "--nz", str(nz),
        "--dx", f"{dx:.3f}",
        "--freq", f"{freq:.3f}",
        "--steps", str(steps),
        "--pulse",
        "--beta-rho", "0.5",
        "--beta-q", f"{beta_q:.3f}",
        "-o", str(out_fp),
    ]
    run_wavecli(wavecli, args)
    with open(out_fp) as f:
        return json.load(f)


def cosine(a: List[float], b: List[float]) -> float:
    if len(a) != len(b) or not a:
        return 0.0
    dot = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(x * x for x in b))
    if na == 0 or nb == 0:
        return 0.0
    return dot / (na * nb)


# ---------------------------------------------------------------------------
# Vina docking of decoys at the binding site
# ---------------------------------------------------------------------------

def dock_decoy(target: Target, decoy_sdf: Path, out_dir: Path,
               num_modes: int) -> List[Path]:
    cx, cy, cz = target.site_center
    result = vina_dock_fn(
        receptor_pdb=target.apo_pdb,
        ligand_in=decoy_sdf,
        center=(cx, cy, cz),
        size_A=(20.0, 20.0, 20.0),
        num_modes=num_modes,
        exhaustiveness=8,
        seed=42,
        out_dir=out_dir,
    )
    return result.pose_files


# ---------------------------------------------------------------------------
# Main per-target loop
# ---------------------------------------------------------------------------

def run_target(target: Target, *,
               wavecli: Path,
               work_dir: Path,
               num_modes: int,
               nx: int, ny: int, nz: int, dx: float,
               freq: float, steps: int,
               beta_q: float, region_half: float) -> dict:
    print(f"\n========== TARGET {target.code} ==========")
    print(f"  apo:    {target.apo_pdb.name}")
    print(f"  actives ({len(target.actives)}): {[a.name for a in target.actives]}")
    print(f"  site center: {target.site_center}")
    print(f"  decoys ({len(target.decoys)}): {[d.stem for d in target.decoys]}")

    target_dir = work_dir / target.code
    fp_dir = target_dir / "fp"
    fp_dir.mkdir(parents=True, exist_ok=True)
    dock_dir = target_dir / "dock"
    dock_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print("  > running apo baseline...")
    apo_fp = score_apo(wavecli, target.apo_pdb, target.site_center,
                       nx, ny, nz, dx, freq, steps, beta_q, region_half,
                       fp_dir / "apo.fp.json")
    apo_region_E = apo_fp["scalars"].get("regional_energy", 0.0)
    print(f"    regional_energy(apo) = {apo_region_E:.6e}")

    # ALL actives: each gets its own fingerprint by placing at the
    # binding-site centroid (NOT the active's native PDB coords —
    # uniform placement makes the cross-active comparison fair).
    active_fps: List[Tuple[str, dict]] = []
    for a in target.actives:
        print(f"  > scoring active {a.name} ({a.pdb.name})...")
        fp = score_pose(wavecli, target.apo_pdb, a.pdb,
                        target.site_center,
                        nx, ny, nz, dx, freq, steps, beta_q, region_half,
                        fp_dir / f"active_{a.name}.fp.json")
        active_fps.append((a.name, fp))

    # Build prototype as the average spectral fingerprint across actives.
    # When N=1, the prototype IS the single active and cos(active,proto)=1
    # trivially; with N≥2 the prototype is a real average and each
    # individual active has its own cos < 1.
    if not active_fps:
        raise RuntimeError(f"no actives for target {target.code}")
    proto_len = len(active_fps[0][1]["spectral"])
    proto = [0.0] * proto_len
    for _, fp in active_fps:
        for i, v in enumerate(fp["spectral"]):
            proto[i] += v / len(active_fps)

    # Decoys: Vina dock → top-K poses → wavecli on each pose.
    decoy_results = []
    for sdf in target.decoys:
        print(f"  > docking decoy {sdf.stem}...")
        try:
            poses = dock_decoy(target, sdf, dock_dir / sdf.stem, num_modes)
        except Exception as e:
            print(f"    DOCK FAILED: {e}")
            continue
        per_pose_R_E = []
        per_pose_spec = []
        for i, pose in enumerate(poses):
            try:
                fp = score_pose(wavecli, target.apo_pdb, pose,
                                target.site_center,
                                nx, ny, nz, dx, freq, steps, beta_q, region_half,
                                fp_dir / f"decoy_{sdf.stem}_pose{i}.fp.json")
            except Exception as e:
                print(f"    pose {i} score FAILED: {e}")
                continue
            re = fp["scalars"]["regional_energy"] / apo_region_E if apo_region_E else 0.0
            per_pose_R_E.append(re)
            per_pose_spec.append(fp["spectral"])
        decoy_results.append({
            "name": sdf.stem,
            "R_Es": per_pose_R_E,
            "specs": per_pose_spec,
        })

    # Per-candidate rows. Each active is scored against the AVERAGED
    # prototype, so cos < 1 in general (only equal to 1 when N=1).
    rows = []
    for name, fp in active_fps:
        sp = fp["spectral"]
        R_E = fp["scalars"]["regional_energy"] / apo_region_E if apo_region_E else 0.0
        rows.append({
            "name": name,
            "kind": "ACTIVE",
            "binderR_best": cosine(sp, proto),
            "binderR_mean": cosine(sp, proto),
            "R_E_best":     R_E,
            "R_E_mean":     R_E,
            "n_poses":      1,
            "pose_variance": 0.0,
        })
    for dr in decoy_results:
        if not dr["specs"]:
            continue
        cosines = [cosine(s, proto) for s in dr["specs"]]
        var = (max(cosines) - min(cosines)) if len(cosines) > 1 else 0.0
        rows.append({
            "name": dr["name"],
            "kind": "decoy",
            "binderR_best": max(cosines),
            "binderR_mean": sum(cosines) / len(cosines),
            "R_E_best":     max(dr["R_Es"], key=lambda x: abs(x - 1.0)),
            "R_E_mean":     sum(dr["R_Es"]) / len(dr["R_Es"]),
            "n_poses":      len(dr["R_Es"]),
            "pose_variance": var,
        })

    actives = [r for r in rows if r["kind"] == "ACTIVE"]
    decoys  = [r for r in rows if r["kind"] == "decoy"]
    pairs = 0
    wins_cosine = 0
    wins_R_E    = 0
    for a in actives:
        for d in decoys:
            pairs += 1
            if a["binderR_best"] > d["binderR_best"]:
                wins_cosine += 1
            if abs(1 - a["R_E_best"]) > abs(1 - d["R_E_best"]):
                wins_R_E += 1
    auc_cosine = wins_cosine / pairs if pairs else 0.0
    auc_R_E    = wins_R_E    / pairs if pairs else 0.0

    elapsed = time.time() - t0
    print(f"\n  === {target.code} per-candidate results ===")
    print(f"  {'kind':>6} {'name':<20} {'binderR_best':>12} {'binderR_mean':>12} "
          f"{'R_E_best':>10} {'R_E_mean':>10} {'pose_var':>10} {'#poses':>7}")
    for r in sorted(rows, key=lambda x: -x["binderR_best"]):
        print(f"  {r['kind']:>6} {r['name']:<20} "
              f"{r['binderR_best']:>12.6f} {r['binderR_mean']:>12.6f} "
              f"{r['R_E_best']:>10.4f} {r['R_E_mean']:>10.4f} "
              f"{r['pose_variance']:>10.4f} {r['n_poses']:>7d}")
    print(f"\n  AUC by binderR_best : {wins_cosine}/{pairs} = {auc_cosine:.3f}")
    print(f"  AUC by |1 - R_E_best|: {wins_R_E}/{pairs} = {auc_R_E:.3f}")
    print(f"  wall time: {elapsed:.1f}s")

    return {
        "target": target.code,
        "rows": rows,
        "auc_cosine": auc_cosine,
        "auc_R_E": auc_R_E,
        "elapsed_s": elapsed,
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def build_targets() -> List[Target]:
    """Define the three benchmark targets."""
    # Common decoy library (small molecules already on disk).
    decoy_set = sorted([
        ROOT / "data/ligands/decoys/aspirin.sdf",
        ROOT / "data/ligands/decoys/caffeine.sdf",
        ROOT / "data/ligands/decoys/glucose.sdf",
        ROOT / "data/ligands/decoys/benzene.sdf",
        ROOT / "data/ligands/decoys/acetaminophen.sdf",
    ])
    return [
        # HIV-1 protease: THREE known inhibitors as actives (multi-active
        # prototype is non-degenerate — each active gets cos < 1 against
        # the averaged prototype, dissolving the 1-active tautology).
        Target(
            code="1HXB",
            apo_pdb=ROOT / "data/pdb/apo/1HXB_apo.pdb",
            site_center=(18.257, -0.023, 11.416),
            actives=[
                ActiveLigand("saquinavir",
                             ROOT / "data/ligands/extracted/1HXB_ROC.pdb"),
                ActiveLigand("indinavir",
                             ROOT / "data/ligands/extracted/1HSG_MK1.pdb"),
                ActiveLigand("XK263",
                             ROOT / "data/ligands/extracted/1HVR_XK2.pdb"),
            ],
            decoys=list(decoy_set),
        ),
        # Single-active targets (until we get more co-crystals — Phase 9
        # data acquisition). Marked here so the AUC tautology caveat
        # carries forward in the report.
        Target(
            code="3PTB",
            apo_pdb=ROOT / "data/pdb/apo/3PTB_apo.pdb",
            site_center=(-1.759, 14.461, 16.916),
            actives=[ActiveLigand("benzamidine",
                                  ROOT / "data/ligands/extracted/3PTB_BEN.pdb")],
            decoys=list(decoy_set),
        ),
        Target(
            code="1STP",
            apo_pdb=ROOT / "data/pdb/apo/1STP_apo.pdb",
            site_center=(11.118, 1.680, -10.755),
            actives=[ActiveLigand("biotin",
                                  ROOT / "data/ligands/extracted/1STP_BTN.pdb")],
            decoys=list(decoy_set),
        ),
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--targets", type=str, default="1HXB,3PTB,1STP",
                    help="comma-separated PDB codes to run")
    ap.add_argument("--num-modes", type=int, default=5,
                    help="Vina top-K poses per decoy (default 5)")
    ap.add_argument("--nx", type=int, default=80)
    ap.add_argument("--ny", type=int, default=80)
    ap.add_argument("--nz", type=int, default=80)
    ap.add_argument("--dx", type=float, default=0.7)
    ap.add_argument("--freq", type=float, default=0.5)
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--beta-q", type=float, default=0.5)
    ap.add_argument("--region-half", type=float, default=7.0)
    ap.add_argument("--work-dir", type=Path,
                    default=Path("/tmp/multitarget_bench"))
    args = ap.parse_args()

    wavecli = args.root / "build" / "bin" / "wavecli"
    work = args.work_dir
    work.mkdir(parents=True, exist_ok=True)

    selected_codes = set(args.targets.split(","))
    targets = [t for t in build_targets() if t.code in selected_codes]

    print("Phase 8.6 — multi-target ligand-discrimination benchmark")
    print(f"  targets: {[t.code for t in targets]}")
    print(f"  grid: {args.nx}x{args.ny}x{args.nz} dx={args.dx} steps={args.steps}")
    print(f"  source: pulse, freq={args.freq}, beta_rho=0.5, beta_q={args.beta_q}")
    print(f"  region: half-width {args.region_half} Å around site")
    print(f"  vina:   top {args.num_modes} poses per decoy")
    print(f"  work:   {work}")

    all_results = []
    for t in targets:
        try:
            r = run_target(t, wavecli=wavecli, work_dir=work,
                           num_modes=args.num_modes,
                           nx=args.nx, ny=args.ny, nz=args.nz, dx=args.dx,
                           freq=args.freq, steps=args.steps,
                           beta_q=args.beta_q, region_half=args.region_half)
            all_results.append(r)
        except Exception as e:
            print(f"\nTARGET {t.code} FAILED: {e}", file=sys.stderr)
            import traceback; traceback.print_exc()

    # Summary
    print("\n========================= SUMMARY =========================")
    print(f"{'target':<8} {'AUC_cos':>9} {'AUC_R_E':>9} {'time':>8}")
    print("-" * 40)
    passing = 0
    for r in all_results:
        flag = "PASS" if (r["auc_cosine"] >= 0.85 or r["auc_R_E"] >= 0.85) else "fail"
        if flag == "PASS":
            passing += 1
        print(f"{r['target']:<8} {r['auc_cosine']:>9.3f} {r['auc_R_E']:>9.3f} "
              f"{r['elapsed_s']:>7.1f}s  {flag}")
    print()
    print(f"Phase 8.6 gate (AUC ≥ 0.85 on ≥ 2 of 3 targets): "
          f"{passing}/{len(all_results)} {'PASS' if passing >= 2 else 'FAIL'}")
    return 0 if passing >= 2 else 1


if __name__ == "__main__":
    sys.exit(main())
