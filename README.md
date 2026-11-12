[![tests](https://github.com/dhwan0112/cm3191-expt4-nmr/actions/workflows/tests.yml/badge.svg)](https://github.com/dhwan0112/cm3191-expt4-nmr/actions/workflows/tests.yml)

# cm3191-expt4-nmr

Python scripts that read processed Bruker ¹H NMR spectra and reproduce the
numbers and figures used for CM3191 Organic Experiment 4, the Grignard
synthesis of (±)-linalool (NUS, AY26/27).

The scripts read the processed spectrum (`pdata/1/1r`) directly. They do not
need TopSpin or MestReNova, and they do not redo the Fourier transform or
phasing.

## What it does

- Reads a Bruker experiment folder or a `.zip` of one.
- References the axis to residual CHCl₃ at δ 7.26 and flattens the baseline.
- Picks line positions and reports first-order coupling constants.
- Integrates one resolved signal per species and applies the mass correction
  from the course handout *Correcting Compound Mass and Yield for Residual
  Solvent Using ¹H NMR*, extended from two to three components.
- Draws the purified spectrum, and crude against purified at true relative
  scale with labelled zoom insets.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Run

Raw data are not in this repository. Put the two spectrometer exports in
`data/` (ignored by git), then:

```bash
python scripts/run_expt4.py \
    --crude    data/CM3191-Week_11-Thurs-400-Oct29-2026-C05.zip \
    --purified data/CM3191-Week_11-Thurs-400-Oct29-2026-P05.zip \
    --out results
```

This prints a summary and writes `results/results.json`,
`results/nmr_purified.png` and `results/nmr_crude_vs_purified.png`.
Masses and the ketone volume default to the values on the signed datasheet
and can be changed with `--crude-mass`, `--purified-mass` and `--ketone-ml`.

Expected summary for the Expt 4 data:

| Quantity | Value |
|---|---|
| Ketone | 2.98 mmol (0.44 mL × 0.855 g mL⁻¹ ÷ 126.20 g mol⁻¹) |
| Theoretical linalool | 0.460 g |
| Purified yield | 33.7% (0.1549 g) |
| Crude mole ratio, linalool : side product : EtOAc | 1 : 0.57 : 0.08 |
| Linalool in crude | ≈ 66 wt%, 0.311 g |
| Column recovery (linalool basis) | ≈ 50% |

## Method

Mole ratios use one well-resolved, non-overlapping resonance per species,
divided by its number of protons:

| Species | Signal | Window (ppm) | nH |
|---|---|---|---|
| Linalool | CH=CH₂ internal H, δ 5.91 | 5.80–6.00 | 1 |
| Side product | CH(OH), δ 3.81 | 3.70–3.92 | 1 |
| Ethyl acetate | OCH₂, δ 4.12 | 4.05–4.20 | 2 |

The baseline is a quadratic fitted through signal-free regions, and integrals
are plain sums over the windows.

## Limits

- The spectra were acquired with `zg30`, D1 = 1 s and 8 scans. That is not a
  quantitative experiment, so treat the mole ratios as ±10%.
- Only the three listed species are counted. Water, salts or grease are not,
  so the linalool weight fraction is an upper bound.
- The side product is assigned tentatively as 6-methyl-5-hepten-2-ol.
- The integration windows and baseline regions are set for these spectra.
  Check them before using the code on another compound.
- The results have not been cross-checked against MestReNova.

## Tests

```bash
pytest
```

The tests use synthetic spectra and the worked example from the course
handout, so they run without the raw data.

## Layout

```
nmrtools/bruker.py     read pdata/1/1r, referencing, baseline
nmrtools/analysis.py   lines, coupling constants, integrals, mass correction
nmrtools/plots.py      the two figures
scripts/run_expt4.py   command-line entry point for Expt 4
tests/                 synthetic-data tests
```

## License

MIT. See [LICENSE](LICENSE).

The course submission itself (the lab report and datasheets) is not included in this repository.
