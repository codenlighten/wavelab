#!/usr/bin/env python3
"""Real ligand benchmark — does wavelab discriminate a chemical family
from a diverse decoy set via spectral-fingerprint correlation?

Family: 5 steroids (cholesterol, testosterone, estradiol, cortisol,
progesterone) — share the rigid cyclopentanoperhydrophenanthrene
4-ring scaffold.

Decoys: 5 chemically diverse small molecules (aspirin, caffeine,
glucose, benzene, acetaminophen).

Hypothesis: prototype built from the 5 steroids' fingerprints will
correlate more strongly with each steroid than with each decoy.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path


def cosine(a: list[float], b: list[float]) -> float:
    if len(a) != len(b):
        return 0.0
    dot = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(x * x for x in b))
    if na == 0 or nb == 0:
        return 0.0
    return dot / (na * nb)


def build_prototype(fps: list[dict]) -> list[float]:
    if not fps:
        return []
    n = len(fps[0]["spectral"])
    proto = [0.0] * n
    for fp in fps:
        for i, v in enumerate(fp["spectral"]):
            proto[i] += v / len(fps)
    return proto


def run_wavecli(wavecli: Path, sdf: Path, out: Path,
                common_args: list[str]) -> None:
    cmd = [str(wavecli), "--sdf", str(sdf), "-o", str(out)] + common_args
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(
            f"wavecli failed on {sdf}:\nstdout:\n{r.stdout}\nstderr:\n{r.stderr}"
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("/root/wavelab/wavelab"),
                    help="wavelab repo root on the runner")
    ap.add_argument("--nx", type=int, default=80)
    ap.add_argument("--ny", type=int, default=80)
    ap.add_argument("--dx", type=float, default=0.4)
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--freq", type=float, default=3.0)
    args = ap.parse_args()

    root = args.root
    wavecli = root / "build" / "bin" / "wavecli"
    data = root / "data" / "ligands"
    work = data / "_fp"
    work.mkdir(parents=True, exist_ok=True)

    actives = sorted((data / "actives").glob("*.sdf"))
    decoys = sorted((data / "decoys").glob("*.sdf"))
    if not actives:
        print(f"no actives in {data / 'actives'}", file=sys.stderr)
        return 2
    if not decoys:
        print(f"no decoys in {data / 'decoys'}", file=sys.stderr)
        return 2

    common = ["--nx", str(args.nx), "--ny", str(args.ny),
              "--dx", str(args.dx), "--steps", str(args.steps),
              "--freq", str(args.freq)]

    print("actives:", [p.stem for p in actives])
    print("decoys: ", [p.stem for p in decoys])
    print("common args:", " ".join(common))
    print()

    print("running wavecli on all ligands...")
    all_fps: dict[str, dict] = {}
    for sdf in actives + decoys:
        out = work / (sdf.stem + ".fp.json")
        run_wavecli(wavecli, sdf, out, common)
        with open(out) as f:
            all_fps[sdf.stem] = json.load(f)
    print("done\n")

    active_names = {p.stem for p in actives}
    active_fps = [all_fps[n] for n in (p.stem for p in actives)]
    proto = build_prototype(active_fps)

    scored = []
    for name, fp in all_fps.items():
        is_active = name in active_names
        r = cosine(fp["spectral"], proto)
        scored.append({
            "score": r,
            "name": name,
            "kind": "ACTIVE" if is_active else "decoy",
            "energy": fp["scalars"].get("total_energy", float("nan")),
            "entropy": fp["scalars"].get("entropy", float("nan")),
            "atoms": int(fp["meta"].get("atoms", 0)),
        })
    scored.sort(key=lambda r: r["score"], reverse=True)

    print(f"{'rank':>4} {'kind':>6} {'ligand':<15} {'binderR':>9} "
          f"{'atoms':>6} {'total_E':>10} {'entropy':>9}")
    print("-" * 70)
    for i, row in enumerate(scored, 1):
        marker = " **" if row["kind"] == "ACTIVE" else "   "
        print(f"{i:>4} {row['kind']:>6} {row['name']:<15} "
              f"{row['score']:>7.4f}{marker} "
              f"{row['atoms']:>6d} {row['energy']:>10.3f} {row['entropy']:>9.4f}")

    # AUC-style discrimination: fraction of (active, decoy) pairs where
    # the active outscores the decoy.
    pairs = 0
    correct = 0
    for a in scored:
        if a["kind"] != "ACTIVE":
            continue
        for b in scored:
            if b["kind"] != "decoy":
                continue
            pairs += 1
            if a["score"] > b["score"]:
                correct += 1
    auc = correct / pairs if pairs else 0.0
    print()
    print(f"discrimination AUC: {correct}/{pairs} active>decoy = {auc:.3f}")
    print(f"  1.0 = perfect separation, 0.5 = random, <0.5 = inverted")

    # Mean score per group for an additional easy-to-read summary.
    act_mean = sum(r["score"] for r in scored if r["kind"] == "ACTIVE") / len(active_fps)
    dec_mean = sum(r["score"] for r in scored if r["kind"] == "decoy") / len(decoys)
    print(f"mean binderR  actives: {act_mean:.4f}   decoys: {dec_mean:.4f}   "
          f"separation: {act_mean - dec_mean:+.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
