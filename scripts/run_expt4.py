#!/usr/bin/env python3
"""Reproduce the 1H NMR numbers and figures for CM3191 Organic Expt 4,
the Grignard synthesis of (±)-linalool.

Example
-------
    python scripts/run_expt4.py \
        --crude    data/CM3191-Week_11-Thurs-400-Oct29-2026-C05.zip \
        --purified data/CM3191-Week_11-Thurs-400-Oct29-2026-P05.zip \
        --out results
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nmrtools import Component, Spectrum, coupling_constants, integrate, load_spectrum, mass_correction, pick_lines
from nmrtools.plots import plot_crude_vs_purified, plot_purified

MW = {"linalool": 154.25, "ketone": 126.20, "side product": 128.21, "EtOAc": 88.11}

# One well-resolved, non-overlapping resonance per species (handout, section 1).
# (label, window low, window high, number of H)
QUANT_SIGNALS = {
    "linalool": ("CH=CH2 internal H, δ 5.91", 5.80, 6.00, 1),
    "side product": ("CH(OH), δ 3.81", 3.70, 3.92, 1),
    "EtOAc": ("OCH2, δ 4.12", 4.05, 4.20, 2),
}


def prepare(path: Path) -> Spectrum:
    """Load, reference to CHCl3 at 7.26 ppm, and flatten the baseline."""
    return load_spectrum(path).referenced().baseline_offset().baseline_poly()


def purified_assignment(spec: Spectrum) -> list[dict]:
    """Shifts and coupling constants of the purified product."""
    rows = []
    for label, lo, hi in (("Ha", 5.86, 5.96), ("Hb", 5.18, 5.25), ("Hc", 5.03, 5.095)):
        centre, js = coupling_constants(pick_lines(spec, lo, hi, n=4), spec.sf)
        rows.append({"label": label, "delta": round(centre, 3), "mult": "dd", "J_Hz": [round(j, 1) for j in js]})
    hd = pick_lines(spec, 5.095, 5.16, min_rel_height=0.3)
    rows.append({"label": "Hd", "delta": round(sum(hd) / len(hd), 3), "mult": "m", "J_Hz": []})
    he = pick_lines(spec, 1.90, 2.20, min_rel_height=0.2)
    rows.append({"label": "He", "delta": [round(min(he), 2), round(max(he), 2)], "mult": "m", "J_Hz": []})
    for label, lo, hi in (("Hf (CH3)", 1.665, 1.70), ("Hg (CH3)", 1.59, 1.615)):
        rows.append({"label": label, "delta": round(pick_lines(spec, lo, hi, n=1)[0], 3), "mult": "s", "J_Hz": []})
    hh = pick_lines(spec, 1.50, 1.595, min_rel_height=0.2)
    rows.append({"label": "Hh", "delta": [round(min(hh), 2), round(max(hh), 2)], "mult": "m", "J_Hz": []})
    rows.append({"label": "Hi (CH3)", "delta": round(pick_lines(spec, 1.26, 1.30, n=1)[0], 3), "mult": "s", "J_Hz": []})
    return rows


def quantify(spec: Spectrum, total_mass: float) -> dict:
    """Integrate the quantification signals and apply the mass correction to `total_mass` (g)."""
    comps = [
        Component(name, integrate(spec, lo, hi), n_h, MW[name]) for name, (_, lo, hi, n_h) in QUANT_SIGNALS.items()
    ]
    comps = [c if c.integral > 0 else Component(c.name, 0.0, c.n_h, c.mw) for c in comps]
    return mass_correction(comps, total_mass, target="linalool")


def main() -> None:
    """Parse the command line, run the analysis, and write results.json and the two figures."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--crude", type=Path, required=True, help="Bruker folder or .zip of the crude sample")
    p.add_argument("--purified", type=Path, required=True, help="Bruker folder or .zip of the purified sample")
    p.add_argument("--out", type=Path, default=Path("results"))
    p.add_argument("--ketone-ml", type=float, default=0.44)
    p.add_argument("--ketone-density", type=float, default=0.855, help="g/mL")
    p.add_argument("--crude-mass", type=float, default=0.4725, help="g")
    p.add_argument("--purified-mass", type=float, default=0.1549, help="g")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    crude, purified = prepare(args.crude), prepare(args.purified)

    ketone_mmol = args.ketone_ml * args.ketone_density / MW["ketone"] * 1000
    theoretical_g = ketone_mmol / 1000 * MW["linalool"]

    q_crude = quantify(crude, args.crude_mass)
    q_pure = quantify(purified, args.purified_mass)
    lin = q_crude["components"]["linalool"]
    side = q_crude["components"]["side product"]

    results = {
        "basis": {
            "ketone_mmol": round(ketone_mmol, 3),
            "theoretical_linalool_g": round(theoretical_g, 4),
            "note": "0.44 mL x 0.855 g/mL / 126.20 g/mol; the manual table rounds this to 3.0 mmol",
        },
        "acquisition": {"crude": crude.meta, "purified": purified.meta},
        "isolated": {
            "crude_mass_g": args.crude_mass,
            "crude_yield_pct": round(args.crude_mass / theoretical_g * 100, 1),
            "purified_mass_g": args.purified_mass,
            "purified_yield_pct": round(args.purified_mass / theoretical_g * 100, 1),
        },
        "crude_composition": q_crude,
        "purified_composition": q_pure,
        "derived": {
            "linalool_in_crude_g": round(lin["mass_g"], 3),
            "linalool_in_crude_wt_pct": round(lin["wt_fraction"] * 100, 1),
            "corrected_yield_pct": round(lin["mass_g"] / theoretical_g * 100, 1),
            "side_product_mmol": round(side["mmol"], 2),
            "column_recovery_pct": round(args.purified_mass / lin["mass_g"] * 100, 1),
            "mass_balance_pct_of_ketone": round((lin["mmol"] + side["mmol"]) / ketone_mmol * 100, 1),
        },
        "purified_assignment": purified_assignment(purified),
        "quant_signals": {
            k: {"signal": v[0], "window_ppm": [v[1], v[2]], "n_H": v[3]} for k, v in QUANT_SIGNALS.items()
        },
        "caveats": [
            "zg30 with a short relaxation delay is not quantitative; treat mole ratios as +/- 10 %.",
            "Only the listed components are counted, so the linalool weight fraction is an upper bound.",
            "The side product is tentatively 6-methyl-5-hepten-2-ol (MW 128.21).",
        ],
    }
    (args.out / "results.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))

    labels = {
        r["delta"]: r["label"].split()[0] for r in results["purified_assignment"] if isinstance(r["delta"], float)
    }
    plot_purified(purified, args.out / "nmr_purified.png", labels)
    plot_crude_vs_purified(crude, purified, args.out / "nmr_crude_vs_purified.png")

    # -- console summary ----------------------------------------------------
    d, iso = results["derived"], results["isolated"]
    ratio = " : ".join(f"{q_crude['components'][k]['mole_ratio']:.2f}" for k in QUANT_SIGNALS)
    print(f"Ketone                {ketone_mmol:.2f} mmol   theoretical linalool {theoretical_g:.3f} g")
    print(f"Crude                 {args.crude_mass:.4f} g ({iso['crude_yield_pct']} %)")
    print(f"Purified              {args.purified_mass:.4f} g ({iso['purified_yield_pct']} %)")
    print(f"Crude mole ratio      linalool : side product : EtOAc = {ratio}")
    print(
        f"Linalool in crude     {d['linalool_in_crude_wt_pct']} wt%  = {d['linalool_in_crude_g']} g  (corrected yield {d['corrected_yield_pct']} %)"
    )
    print(
        f"Side product          {d['side_product_mmol']} mmol   mass balance {d['mass_balance_pct_of_ketone']} % of ketone"
    )
    print(f"Column recovery       {d['column_recovery_pct']} % (linalool basis)")
    print("Purified assignment:")
    for r in results["purified_assignment"]:
        js = ", ".join(f"{j:.1f}" for j in r["J_Hz"])
        delta = f"{r['delta'][0]:.2f}–{r['delta'][1]:.2f}" if isinstance(r["delta"], list) else f"{r['delta']:.2f}"
        print(f"  {r['label']:9s} δ {delta:11s} {r['mult']:2s} {js}")
    print(f"Wrote {args.out}/results.json and two figures.")


if __name__ == "__main__":
    main()
