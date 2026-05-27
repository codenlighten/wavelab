#!/usr/bin/env python3
"""Pocket + ligand benchmark — overview §15 R_E concept.

For each candidate ligand:
  1. Run wavecli on the pocket alone (cached after first call).
  2. Run wavecli on the pocket + ligand (translated to pocket centroid
     unless --place-at is given).
  3. Compute:
       * R_E    = total_energy(p+l) / total_energy(p)
       * pert   = 1 - cos(spectrum(p+l), spectrum(p))   ∈ [0, 2]

Caveat: this benchmark places all ligands at a single chosen point
(default: pocket centroid). It's not a docking score — it's "how
differently does each ligand perturb the same probe configuration?"
Discrimination signals here support the §15 architecture, but real
binding-affinity ranking would need docked poses.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path


def cosine(a: list[float], b: list[float]) -> float:
    if len(a) != len(b) or not a:
        return 0.0
    dot = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(x * x for x in b))
    if na == 0 or nb == 0:
        return 0.0
    return dot / (na * nb)


def run_wavecli(wavecli: Path, args: list[str]) -> None:
    r = subprocess.run([str(wavecli)] + args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(
            f"wavecli failed: {' '.join(args)}\nstdout:\n{r.stdout}\nstderr:\n{r.stderr}"
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("/root/wavelab/wavelab"))
    ap.add_argument("--pdb", type=Path, default=Path("data/pdb/1ubq.pdb"),
                    help="pocket PDB (relative to --root)")
    ap.add_argument("--nx", type=int, default=120)
    ap.add_argument("--ny", type=int, default=120)
    ap.add_argument("--dx", type=float, default=0.4)
    ap.add_argument("--steps", type=int, default=400)
    ap.add_argument("--freq", type=float, default=0.5)
    ap.add_argument("--pulse", action="store_true", default=True)
    ap.add_argument("--no-pulse", action="store_false", dest="pulse")
    ap.add_argument("--beta-rho", type=float, default=0.5)
    ap.add_argument("--beta-q", type=float, default=0.5)
    ap.add_argument("--place-at", type=str, default=None,
                    help="X,Y,Z world coord to center the ligand at "
                         "(default: pocket centroid)")
    args = ap.parse_args()

    root = args.root
    wavecli = root / "build" / "bin" / "wavecli"
    pocket = root / args.pdb
    if not pocket.exists():
        print(f"pocket PDB not found: {pocket}", file=sys.stderr)
        return 2

    data = root / "data" / "ligands"
    work = root / "data" / "ligands" / "_fp_pkt"
    work.mkdir(parents=True, exist_ok=True)

    actives = sorted((data / "actives").glob("*.sdf"))
    decoys = sorted((data / "decoys").glob("*.sdf"))

    # Probe defaults to (nx/2, ny/2) so we sample the field right where
    # the ligand was placed (default placement = scene centroid, which
    # maps to ~ mid-grid in cell coords). That gives much sharper signal
    # than the standard downstream probe at 3*nx/4.
    probe_i = args.nx // 2
    probe_j = args.ny // 2

    common = ["--nx", str(args.nx), "--ny", str(args.ny),
              "--dx", str(args.dx), "--steps", str(args.steps),
              "--freq", str(args.freq),
              "--beta-rho", str(args.beta_rho),
              "--probe", f"{probe_i},{probe_j}"]
    if args.beta_q != 0:
        common += ["--beta-q", str(args.beta_q)]
    if args.pulse:
        common.append("--pulse")
    if args.place_at:
        common += ["--place-at", args.place_at]

    print(f"pocket: {pocket}")
    print(f"actives: {[p.stem for p in actives]}")
    print(f"decoys:  {[p.stem for p in decoys]}")
    print(f"common args: {' '.join(common)}")
    print()

    # 1. Pocket alone (cached).
    pocket_fp_path = work / "pocket_alone.fp.json"
    print("running pocket alone...")
    run_wavecli(wavecli, ["--pdb", str(pocket), "-o", str(pocket_fp_path)] + common)
    with open(pocket_fp_path) as f:
        pocket_fp = json.load(f)
    pocket_E = pocket_fp["scalars"]["total_energy"]
    pocket_spec = pocket_fp["spectral"]

    # 2. Pocket + each ligand.
    print("running pocket + each ligand...")
    rows = []
    for sdf in actives + decoys:
        is_active = sdf in actives
        out = work / (sdf.stem + ".fp.json")
        run_wavecli(wavecli, ["--pdb", str(pocket), "--add-sdf", str(sdf),
                              "-o", str(out)] + common)
        with open(out) as f:
            fp = json.load(f)
        E_pl = fp["scalars"]["total_energy"]
        R_E = E_pl / pocket_E if pocket_E else 0.0
        pert = 1.0 - cosine(fp["spectral"], pocket_spec)
        rows.append({
            "name": sdf.stem,
            "kind": "ACTIVE" if is_active else "decoy",
            "R_E": R_E,
            "perturbation": pert,
            "energy": E_pl,
            "atoms": int(fp["meta"].get("atoms", 0)),
        })

    # Sort by perturbation (largest first — most distinctive ligands).
    rows.sort(key=lambda r: r["perturbation"], reverse=True)

    print()
    print(f"{'rank':>4} {'kind':>6} {'ligand':<15} {'R_E':>10} "
          f"{'perturb':>14} {'atoms':>6}")
    print("-" * 70)
    for i, r in enumerate(rows, 1):
        marker = " **" if r["kind"] == "ACTIVE" else "   "
        print(f"{i:>4} {r['kind']:>6} {r['name']:<15} "
              f"{r['R_E']:>10.6f} {r['perturbation']:>12.4e}{marker} "
              f"{r['atoms']:>6d}")

    # Discriminate: build prototype from active p+l fingerprints, score
    # each candidate by cosine sim to that prototype. (Different from
    # the pocket-alone-vs-p+l perturbation — this is family-clustering
    # in the pocket-context.)
    def cluster_score(fp_list: list[dict]) -> list[float]:
        n = len(fp_list[0]["spectral"])
        proto = [0.0] * n
        for fp in fp_list:
            for i, v in enumerate(fp["spectral"]):
                proto[i] += v / len(fp_list)
        return proto

    active_fps = []
    for sdf in actives:
        with open(work / (sdf.stem + ".fp.json")) as f:
            active_fps.append(json.load(f))
    proto = cluster_score(active_fps)

    pairs = 0
    correct = 0
    for row_a in rows:
        if row_a["kind"] != "ACTIVE":
            continue
        with open(work / (row_a["name"] + ".fp.json")) as f:
            fp_a = json.load(f)
        sim_a = cosine(fp_a["spectral"], proto)
        for row_b in rows:
            if row_b["kind"] != "decoy":
                continue
            with open(work / (row_b["name"] + ".fp.json")) as f:
                fp_b = json.load(f)
            sim_b = cosine(fp_b["spectral"], proto)
            pairs += 1
            if sim_a > sim_b:
                correct += 1
    auc = correct / pairs if pairs else 0.0

    print()
    print(f"pocket-context binder AUC (clustering p+l prototypes): "
          f"{correct}/{pairs} = {auc:.3f}")

    # Summary stats.
    act_pert = [r["perturbation"] for r in rows if r["kind"] == "ACTIVE"]
    dec_pert = [r["perturbation"] for r in rows if r["kind"] == "decoy"]
    print(f"perturbation  actives mean: {sum(act_pert)/len(act_pert):.6f}   "
          f"decoys mean: {sum(dec_pert)/len(dec_pert):.6f}")
    act_RE = [r["R_E"] for r in rows if r["kind"] == "ACTIVE"]
    dec_RE = [r["R_E"] for r in rows if r["kind"] == "decoy"]
    print(f"R_E           actives mean: {sum(act_RE)/len(act_RE):.4f}   "
          f"decoys mean: {sum(dec_RE)/len(dec_RE):.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
