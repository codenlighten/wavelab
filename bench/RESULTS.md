# Ligand benchmark — first results

Run on the GPU droplet (RTX 6000 Ada, but CPU `FdtdCpuOmp<2>` is what
the current `wavecli` pipeline uses; CUDA path is exposed only via
the test/bench harnesses for now).

## Setup

- **Actives (5 steroids):** cholesterol, testosterone, estradiol,
  cortisol, progesterone. Share the cyclopentanoperhydrophenanthrene
  4-ring core but differ substantially in pendant groups —
  cholesterol's 8-carbon side chain makes it nearly 30% larger than
  estradiol by atom count.
- **Decoys (5 chemically diverse small molecules):** aspirin,
  caffeine, glucose, benzene, acetaminophen.
- **Source data:** PubChem 3D conformer SDFs, all atoms (heavy + H).
- **Pipeline:** wavecli `--sdf` → density splat → refractive-index
  medium → 2D FDTD → 16-bin log-spaced spectral fingerprint at probe
  past the scatterer.

## Discrimination

`binderR` = cosine similarity between the candidate's spectral
fingerprint and the mean spectral fingerprint of the 5 actives.
`AUC` = fraction of (active, decoy) pairs where the active outscores
the decoy (1.0 = perfect, 0.5 = random).

| Grid     | dx  | freq | source mode | AUC  | mean(active) | mean(decoy) |
| -------- | --- | ---- | ----------- | ---- | ------------ | ----------- |
| 80×80    | 0.4 | 3.0  | harmonic    | 0.32 | 0.9872       | 0.9905      |
| 150×150  | 0.2 | 0.5  | harmonic    | 0.84 | (1.0000 nominal)\* | (1.0000 nominal)\* |
| 150×150  | 0.2 | 0.5  | pulse       | 0.68 | 0.9997       | 0.9998      |

\* The harmonic-source spectrum is dominated by a single peak at the
source frequency, so cosine similarity asymptotes near 1.0 even though
the discrimination is real (it lives in the higher-decimal places).

The first row (`dx=0.4, freq=3`) was a setup bug — `λ = c/f = 0.33`
is *smaller* than `dx = 0.4`, so the wave was severely undersampled
and the engine produced essentially identical noisy fingerprints
regardless of scatterer. Lesson: enforce `dx ≤ λ/10` in any future
benchmark driver.

The two well-resolved rows give qualitatively useful discrimination
but not the >0.9 that a "tight" binder family would yield.

## Per-ligand ranking (pulse, AUC 0.68)

```
rank   kind ligand               binderR  atoms    total_E   entropy
---------------------------------------------------------------------------
   1 ACTIVE cortisol         0.9999807 **     50    305.286    7.6141
   2 ACTIVE progesterone     0.9999562 **     48    292.910    7.6975
   3  decoy caffeine         0.9999501        24    422.186    8.2035
   4 ACTIVE testosterone     0.9999489 **     46    326.398    7.6703
   5  decoy aspirin          0.9999274        21    379.487    8.0792
   6 ACTIVE estradiol        0.9999267 **     44    299.751    7.5614
   7  decoy glucose          0.9998506        24    331.880    7.9620
   8  decoy acetaminophen    0.9997630        20    323.989    7.8313
   9  decoy benzene          0.9995740        12    371.736    8.0595
  10 ACTIVE cholesterol      0.9987194 **     57    421.425    8.4284
```

Note that cholesterol — the largest active by atom count and the only
one with a long flexible side chain — ranks dead last among the
actives. The smaller, ring-only steroids cluster well. The smaller
decoys (benzene, acetaminophen) score lower than the medium-size ones
(caffeine, glucose, aspirin) because their fingerprints have less
in common with the prototype's average. This is shape-driven scoring
in action.

## Lighting up the polar channel — second pass

Wired up the `β_q · |Q(x)|` term from §9. `splat_charge` now runs
when `--beta-q > 0`, using electronegativity-derived per-element
pseudo-charges (H +0.1, C 0, N -0.4, O -0.5, S -0.1, P +0.5,
halogens negative). Real per-molecule Gasteiger / QM partials are
still TODO — these are crude, but they capture the polar-vs-non-polar
ordering correctly.

Sweep over `β_q` × source mode (same 5 actives / 5 decoys, same grid,
`β_rho = 0.5` throughout):

| β_q | source   | AUC  | Δ vs shape-only |
| --- | -------- | ---- | --------------- |
| 0.0 | harmonic | 0.60 | —               |
| 0.5 | harmonic | 0.60 | +0.00           |
| 1.0 | harmonic | 0.56 | −0.04           |
| 2.0 | harmonic | 0.52 | −0.08           |
| 0.0 | pulse    | 0.68 | —               |
| 0.5 | pulse    | **0.76** | **+0.08**   |
| 1.0 | pulse    | 0.76 | +0.08           |
| 2.0 | pulse    | 0.76 | +0.08           |

* **Pulse + polar gives the lift, harmonic doesn't.** Single-frequency
  probing only sees the polar contrast as a shifted wave speed; the
  spectrum stays a delta at the source frequency. Broadband probing
  picks up the polar-induced multipath structure across many
  frequencies.
* **Lift saturates at β_q ≈ 0.5.** Beyond that, the polar channel
  doesn't add new information on this 10-ligand set.

The single ligand that moved most: **caffeine fell from rank #3 to
#6** with the polar channel on. Caffeine has 4 nitrogens in the
xanthine ring — very polar — so it now looks distinct from the
mostly-hydrocarbon-scaffold steroids. Cholesterol stays at the
bottom of the actives at rank #10; that's a real chemistry failure
(its 8-carbon alkyl side chain makes it shape-different from the
other steroids), not an engine failure.

## What's still measured vs what isn't

After this work the engine sees:
- **shape + size** (always, via `ρ(x)`)
- **polarity** (when `--beta-q > 0`, via `|Q(x)|`)
- **hydrophobicity** (when `--beta-h > 0`, via `H(x)`)

What we still don't have for serious use:
- **Per-molecule Gasteiger / QM partial charges.** Per-element
  pseudo-charges miss carbonyl-C polarization, the difference between
  protonated and free amine, conjugation effects, etc. A small
  Gasteiger pass or piping RDKit / OpenBabel charges into the SDF
  parser would fix this.
- **The pocket as the scoring substrate.** All benchmarks above are
  ligand-only. Real binding affinity asks "does this ligand fit
  THIS pocket" — that's the §15 R_E concept (pocket+ligand vs pocket
  alone). Needs a small driver change to splat both into one medium.

## Reproduce

```sh
# On the droplet (RTX 6000 Ada or any box that built wavelab):
cd /root/wavelab/wavelab
python3 bench/ligand_benchmark.py \
    --nx 150 --ny 150 --dx 0.2 \
    --freq 0.5 --steps 400 \
    --pulse --beta-rho 0.5
```

Or harmonic mode (the row that scored 0.84):

```sh
python3 bench/ligand_benchmark.py \
    --nx 150 --ny 150 --dx 0.2 \
    --freq 0.5 --steps 400
```

Driver script honors `--root`, `--nx`, `--ny`, `--dx`, `--steps`,
`--freq`, `--pulse`, `--beta-rho`.
