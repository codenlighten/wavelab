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

## What's actually being measured

The current engine consumes only atomic positions and van-der-Waals
radii (Bondi 1964). Charge and hydrophobicity are loaded into the
`Atom` struct from the element parameter table but currently are NOT
wired into the medium assembly (only `ρ(x)` contributes to `n(x)` —
the `Q(x)` and `H(x)` terms in overview §9 are `β_q = β_h = 0` by
default).

So this benchmark is a **shape + size** fingerprint, not a chemistry
fingerprint. Two ways to make it more useful:

1. **Tighten the family** — pick actives that are genuinely
   shape-similar at the engine's resolution (e.g. the 4-ring
   scaffold only, no pendant chains; or all kinase inhibitors of a
   similar size class).
2. **Light up the charge / hydro channels** — wire `β_q > 0` so
   polar pockets shape the medium differently from hydrophobic
   pockets. The infrastructure is already there
   (`build_medium_from_fields` accepts the optional `charge` and
   `hydro` fields); just needs partial-charge data per atom and the
   wavecli flags to enable it.

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
