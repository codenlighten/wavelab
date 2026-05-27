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

## Phase 9.1 follow-up: multi-active prototype for HIV-1 protease

The single-active concern materialized when we re-ran 1HXB with all
three known HIV protease inhibitors as actives — saquinavir (ROC),
indinavir (MK1, from 1HSG), and XK263 (XK2, from 1HVR) — all placed
at the 1HXB binding-site centroid, with the averaged spectral
fingerprint as a non-degenerate prototype.

Result with N=3 actives (5 decoys × 3 actives = 15 pairs):

| Rank | Kind   | Candidate     | binderR_best | R_E_best |
| ---: | :----- | :------------ | -----------: | -------: |
| 1    | decoy  | acetaminophen |    0.999511  |   0.8979 |
| 2    | ACTIVE | XK263         |    0.999179  |   0.8656 |
| 3    | ACTIVE | indinavir     |    0.998824  |   0.8917 |
| 4    | decoy  | aspirin       |    0.998233  |   0.9093 |
| 5    | decoy  | caffeine      |    0.997879  |   0.9000 |
| 6    | decoy  | glucose       |    0.997194  |   0.9179 |
| 7    | ACTIVE | saquinavir    |    0.996696  |   0.4136 |
| 8    | decoy  | benzene       |    0.995090  |   0.9394 |

**AUC_cosine:  9/15 = 0.600**  (was 1.000 with N=1, but tautological)
**AUC_R_E:    15/15 = 1.000**  (regional energy perturbation perfectly
                                separates the three drugs from the five
                                decoys)

The 1.000 → 0.600 collapse in cosine AUC under multi-active is the
honest reading: **the spectral cosine fingerprint, at this grid /
source / probe configuration, is not strongly discriminative for HIV
protease inhibitors against this 5-molecule decoy set.** The R_E AUC
of 1.000 IS real signal — saquinavir / indinavir / XK263 all displace
more binding-region energy than acetaminophen / aspirin / benzene /
caffeine / glucose. The mechanism is unsurprising: drug-class
inhibitors have 40–50 heavy atoms; decoys have 9–24.

So the engine is doing science:
* **R_E is the working discriminator** for the bigger-active case.
  It's a clean physical scalar — "how much does this ligand disrupt
  the energy field in the binding pocket" — and it ranks correctly.
* **Spectral cosine is weak here** because the protein dominates the
  fingerprint and the ligands' spectral contributions all sit on the
  same "background" envelope.

Most important Phase 9 follow-ups:
1. **Multi-active prototypes for 3PTB and 1STP** (need more
   co-crystals than the single 3PTB+BEN / 1STP+BTN we have today).
2. **Higher-resolution probes for small ligands** — benzamidine
   and biotin are only 9 and 16 atoms; the current 0.85 Å grid +
   0.5 Hz source has wavelength much larger than the ligand. Smaller
   dx and higher freq (better wavelength-to-feature matching) would
   give them a chance.
3. **Hybrid score combining cosine + R_E** — `HybridScore` already
   exists in src/score/hybrid.hpp but isn't yet wired into the
   benchmark. With weights chosen on calibration data, combined
   cosine+R_E should outperform either alone.

## Phase 9.2 follow-up: small-ligand grid tune for 3PTB (trypsin)

Hypothesis: 3PTB's R_E AUC of 0.2 at the low-resolution config
(`dx=0.85, freq=0.5`, wavelength only 2.4 cells) was resolution-
limited; matching wavelength to the ~2.5 Å feature scale where
benzamidine vs aspirin/caffeine differences live should reveal
real signal.

Tuned config: `dx=0.3, freq=0.4, steps=1500, nx=128`. Wavelength
2.5 Å = 8.3 cells/wavelength (well-resolved). Steps bumped from 400
to 1500 so the wave actually propagates from the corner source to
the binding site (source-to-pocket ≈ 43 Å, needs ≥ 43 time units;
old 400 steps × 0.069 dt = 27 time units was insufficient and gave
all-machine-epsilon regional energies).

Result (single-active prototype caveat still applies):

| Rank | Kind   | Candidate     | binderR_best | R_E_best | wall-time |
| ---: | :----- | :------------ | -----------: | -------: | --------: |
| 1    | ACTIVE | benzamidine   |    1.000000  |   1.0728 | 856s      |
| 2    | decoy  | benzene       |    0.999998  |   1.0456 |           |
| 3    | decoy  | glucose       |    0.999998  |   1.0777 |           |
| 4    | decoy  | acetaminophen |    0.999996  |   1.0738 |           |
| 5    | decoy  | aspirin       |    0.999964  |   1.0669 |           |
| 6    | decoy  | caffeine      |    0.999956  |   1.0477 |           |

**AUC_R_E: 3/5 = 0.600** (was 0.200 at low res)
**AUC_cosine: 5/5 = 1.000** (still tautological — only 1 active)

Two qualitative observations:
1. **Improvement is real but modest.** R_E AUC tripled (0.2 → 0.6),
   confirming that 3PTB's low-res failure was *partly* resolution-
   limited. Benzamidine now ranks above caffeine and benzene; it
   still loses to glucose and acetaminophen, which are similarly-
   sized polar molecules.
2. **R_E sign flipped vs the low-res HIV protease case.** At
   `dx=0.3`, all ligands give `R_E > 1` — they *increase* the
   binding-region energy via refractive concentration of the wave.
   At `dx=0.85` (low res), the same physics gave `R_E < 1` —
   ligands acted as absorbers/dispersers. Different physical regime,
   driven by whether the wavelength is large or small relative to
   the atom-spacing features. Not a bug; a real engine response.

The fundamental message: **R_E magnitude tracks ligand size + polarity
more than it tracks "is this molecule a binder for THIS protein."**
A trypsin-specific inhibitor (benzamidine) and a non-specific polar
molecule (glucose, acetaminophen) of comparable size give comparable
R_E. The engine is doing physics; it isn't doing biology.

What WOULD give biology-relevant discrimination:
* **Multi-probe spectra (Fingerprint v2, deferred 8.3b).** A single
  probe captures only one slice of the scattered field; multiple
  probes around the pocket would expose how the ligand reshapes the
  full scattering pattern.
* **Pocket-shape-matched scoring**, not point-source scoring. The
  current setup uses one source far from the pocket — a uniform
  plane-wave illumination from many angles, or a source PLACED IN
  the pocket, would give different (and probably more biologically-
  meaningful) responses.
* **Calibrated multi-target prototypes per ligand class.** Build the
  "trypsin-like binder" prototype from 4-5 known trypsin inhibitors
  with varied scaffolds (the multi-active approach extended to
  small ligands).

## Phase 9.4 follow-up: in-pocket source

Hypothesis: corner-source-far-from-pocket geometry is the reason
spectral discrimination is weak; placing the source AT the
binding-site centroid puts the ligand in the source's near-field
where each atom matters proportionally more.

Result (all three targets, in-pocket source, low-res 64³ grid):

| Target | actives (N) | AUC_cos       | AUC_R_E       | active R_E | decoy R_E range |
| :----- | ----------: | ------------: | ------------: | ---------: | --------------: |
| 1HXB   | 3 (multi)   | 0.333 (↓0.27) | 0.533 (↓0.47) | varied     | varied          |
| 3PTB   | 1           | 1.000 (taut.) | 0.200         | 0.0200     | 0.0100–0.0282   |
| 1STP   | 1           | 1.000 (taut.) | 0.200         | 0.0253     | 0.0102–0.0317   |

Signal magnitude EXPLODED: R_E is now 0.01–0.03 across the board
(vs 0.4–1.0 with corner-source). The ligand is absorbing/scattering
~97–99% of the source's near-field energy. But the residual 1–3%
discrimination doesn't favor actives.

3PTB R_E ranking by |1 - R_E| (biggest perturbation first):
   caffeine 0.990 → aspirin 0.989 → glucose 0.988 → acetaminophen 0.983
   → benzamidine 0.980 → benzene 0.972

**Benzamidine is 5th of 6 by perturbation magnitude.** It beats
only benzene (the smallest decoy). The four decoys with more atoms
than benzamidine all out-perturb it.

**The honest mechanism:** in-pocket regional R_E measures
*atoms-in-the-source-near-field*. Caffeine has 24 atoms in the
binding region; benzamidine has 9. The bigger molecule blocks more
of the source's outward propagation, regardless of whether it's a
real binder for THIS pocket.

So in-pocket source amplifies the **physics** signal (much bigger
energy displacement per ligand) but doesn't surface **biology**
signal (which ligand specifically fits this pocket). At this engine
chemistry — density + polar via β_q — the engine genuinely doesn't
know about hydrogen-bond geometry, π-stacking, or shape
complementarity beyond raw atom-density occupancy.

What this rules out:
* In-pocket source isn't a quick fix for small-ligand discrimination.
  The geometry change amplifies signal but tracks the same "atoms
  in the way" physics.

What's still on the table for actually discriminating biology:
* **Multi-probe spectral with probes around the pocket boundary**
  (8.3b deferred). Different from regional R_E, which is a scalar.
  A multi-probe spectrum captures HOW the ligand reshapes the
  outgoing wave at multiple angles — closer to a "scattering
  pattern" measurement than a "blocking" measurement.
* **Hydrogen-bond-geometry-aware chemistry features.** The current
  per-element pseudo-charges miss carbonyl polarity, amide H-bond
  donors, aromatic π-systems. Real Gasteiger or QM partials would
  let polarity actually distinguish "amidine that protonates and
  H-bonds to Asp" from "benzene that doesn't."
* **Probe at SPECIFIC chemistry locations**, not just centroid.
  Trypsin's S1 pocket has Asp189 at the bottom — a probe there
  would respond to ligands' positive charge specifically.

## Reproduce

On the droplet:

```sh
cd /root/wavelab/wavelab
python3 bench/multitarget_benchmark.py --num-modes 3 \
  --nx 64 --ny 64 --nz 64 --dx 0.85 --steps 150
```

Total wall time ~5 min (90 s per target × 3) on the contended droplet.
With CUDA 3D wired in, this should drop to under a minute total.
