# Phase 8 multi-target benchmark — first 3-target run

Run on the GPU droplet (RTX 6000 Ada, but using `FdtdCpuOmp<3>` —
GPU 3D backend exists but `wavecli` still uses CPU path). 64³ grid,
dx 0.85 Å, 150 steps, pulse source freq 0.5, β_rho=0.5, β_q=0.5,
7 Å region half-width.

| Target | Active (HETATM) | AUC_cos | AUC_R_E | active R_E | decoy R_E range | Verdict |
|---|---|---|---|---|---|---|
| **1HXB** | saquinavir (49 heavy atoms, 31 torsions) | **1.000** | **1.000** | 0.41 | 0.90–0.94 | **PASS** |
| **3PTB** | benzamidine (9 atoms, rigid) | **1.000** | 0.200 | 0.98 | 0.90–0.99 | PASS (cosine) |
| **1STP** | biotin (16 atoms) | **1.000** | 0.600 | 0.86 | 0.79–0.94 | PASS (cosine) |

**Phase 8.6 gate (AUC ≥ 0.85 on ≥ 2 of 3 targets): 3/3 PASS.**

## Honest reading

The AUC_cosine "1.000" requires interpretation. With ONE active per
target, the prototype IS the active's fingerprint, so
`cos(active, prototype) = 1.0` by construction. AUC is 1.0 iff no
decoy hits exactly 1.0 in cosine.

| Target | active cos | decoy cos range | separation |
|---|---|---|---|
| 1HXB | 1.000000 | 0.983793–0.993722 | real (~0.01) |
| 3PTB | 1.000000 | 0.999944–1.000000 | computational noise (~5e-5) |
| 1STP | 1.000000 | 1.000000 | none (cosines tied) |

So **only 1HXB shows real spectral-fingerprint discrimination**.
3PTB and 1STP "pass" the AUC gate on float-roundoff ordering — the
fingerprints are essentially indistinguishable from the active's
fingerprint at this grid/step/source configuration. A multi-active
prototype (more co-crystal PDBs per target) would dissolve the
tautology and give a real comparison.

The regional R_E result is more chemically honest:
- 1HXB: saquinavir is a much bigger ligand (49 heavy atoms) than
  the decoys (9–18 atoms) → bigger pocket perturbation → R_E drops
  to 0.41 vs decoys at 0.90+. Strong signal.
- 3PTB: benzamidine is SMALLER than most decoys → its perturbation
  is small, comparable to (or less than) the decoys. R_E AUC is 0.2.
- 1STP: biotin is mid-sized; mixed signal (R_E AUC 0.6).

## What this validates

The Phase 8 architecture works end-to-end on real co-crystal data:

1. **3D wavecli** (8.1) consumes real protein + ligand PDBs and
   produces fingerprints without 2D-slicing artifacts.
2. **PDB parser hardening** (8.4) extracts the right ligand from
   each co-crystal (ROC/MK1/XK2/BEN/BTN), filters alt-locs, drops
   waters/ions cleanly.
3. **AutoDock Vina** (8.5) docks decoys at known binding sites,
   producing top-K poses we score directly through wavelab.
4. **Regional R_E** (8.3a) localizes the signal — total R_E was
   pinned at ~1.0 for everything in the previous benchmark; regional
   R_E swings 0.41–0.99 across these candidates.
5. **GPU 3D** (8.2) compiled cleanly and parity-tested vs CPU 3D
   (FdtdCuda<3> agrees with FdtdCpuOmp<3> within 0.5% per-cell
   error). Not yet wired into `wavecli` — that's a one-line CMake
   switch + scoring substitution for the next iteration.

## What's NOT validated

* **Spectral fingerprint discrimination for small ligands.** When
  the active and decoys are similar sizes (benzamidine in 3PTB,
  biotin in 1STP), the spectral cosine cannot distinguish them. The
  16-log-bin spectrum is too coarse, OR the wavelength isn't matched
  to the molecular scale, OR (most likely) we need multi-probe
  spectra (deferred to 8.3b — Fingerprint v2 schema bump).
* **Multi-active prototype.** With one crystal per target the
  "prototype" is degenerate. Pulling 3–5 co-crystals per target
  (1HXB / 1HSG / 1HVR all for HIV protease, then averaging) would
  give an honest cosine separation across the active set.

## Headline trajectory across Phase 8

| Configuration | Best result |
|---|---|
| Phase 7: shape+polar fingerprint, ligand-alone, tight steroid family | AUC 1.000 (small synthetic set) |
| Phase 8 start: 2D pocket+ligand, total R_E | not discriminative (R_E pinned near 1.0) |
| Phase 8.3a: 2D pocket+ligand, regional R_E on 1ubq | 16,000× signal amplification per single-ligand test |
| Phase 8.6: 3D pocket+ligand, real co-crystals × 3 targets | **3/3 gate pass; 1 target with genuine discrimination** |

## Reproduce

On the droplet:

```sh
cd /root/wavelab/wavelab
python3 bench/multitarget_benchmark.py --num-modes 3 \
  --nx 64 --ny 64 --nz 64 --dx 0.85 --steps 150
```

Total wall time ~5 min (90 s per target × 3) on the contended droplet.
With CUDA 3D wired in, this should drop to under a minute total.
