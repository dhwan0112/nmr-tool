"""Peak positions, coupling constants, integrals and the mass correction
described in the CM3191 handout "Correcting Compound Mass and Yield for
Residual Solvent Using 1H NMR"."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks

from .bruker import Spectrum


def _refine(spec: Spectrum, i: int) -> float:
    """Sub-point peak position from a parabola through three points."""
    a, b, c = spec.y[i - 1], spec.y[i], spec.y[i + 1]
    denom = a - 2 * b + c
    d = 0.0 if denom == 0 else 0.5 * (a - c) / denom
    return float(spec.ppm[i] + d * (spec.ppm[i + 1] - spec.ppm[i]))


def pick_lines(
    spec: Spectrum,
    lo: float,
    hi: float,
    n: int | None = None,
    min_rel_height: float = 0.15,
) -> list[float]:
    """Line positions (ppm, descending) between lo and hi.

    `n` keeps only the n tallest lines. `min_rel_height` is relative to the
    tallest point in the window.
    """
    idx = np.where(spec.window(lo, hi))[0]
    peaks, _ = find_peaks(spec.y[idx], height=min_rel_height * spec.y[idx].max())
    peaks = idx[peaks]
    if n is not None:
        peaks = sorted(peaks, key=lambda i: -spec.y[i])[:n]
    return sorted((_refine(spec, i) for i in peaks), reverse=True)


def coupling_constants(lines: list[float], sf: float) -> tuple[float, tuple[float, ...]]:
    """Centre (ppm) and coupling constants (Hz, largest first) of a first-order
    doublet (2 lines) or doublet of doublets (4 lines)."""
    lines = sorted(lines, reverse=True)
    centre = float(np.mean(lines))
    if len(lines) == 2:
        return centre, ((lines[0] - lines[1]) * sf,)
    if len(lines) == 4:
        j_small = ((lines[0] - lines[1]) + (lines[2] - lines[3])) / 2 * sf
        j_large = ((lines[0] - lines[2]) + (lines[1] - lines[3])) / 2 * sf
        return centre, (j_large, j_small)
    raise ValueError(f"expected 2 or 4 lines, got {len(lines)}")


def integrate(spec: Spectrum, lo: float, hi: float) -> float:
    """Sum of intensities between lo and hi (baseline assumed already at zero).

    The unit is arbitrary; only ratios of integrals from one spectrum mean
    anything.
    """
    return float(spec.y[spec.window(lo, hi)].sum())


@dataclass(frozen=True)
class Component:
    """One species in the mixture, for the mass correction."""

    name: str
    integral: float  # raw integral of the chosen resonance
    n_h: int  # number of protons behind that resonance
    mw: float  # g/mol

    @property
    def n(self) -> float:
        """Normalised molar value N = I / nH (handout, step 1)."""
        return self.integral / self.n_h


def mass_correction(components: list[Component], total_mass: float, target: str | None = None) -> dict:
    """Handout steps 1-4, for any number of components.

    Returns mole ratios relative to the target (default: first component),
    weight fractions, and the corrected mass and amount of every component.
    The result only accounts for the listed components: anything without a
    counted 1H signal (water, salts, grease) makes the target's weight
    fraction an upper bound.
    """
    target = target or components[0].name
    n_target = next(c.n for c in components if c.name == target)
    rel_mass = {c.name: c.n * c.mw for c in components}  # step 2
    total = sum(rel_mass.values())
    out: dict = {"target": target, "total_mass_g": total_mass, "components": {}}
    for c in components:
        wt = rel_mass[c.name] / total  # step 3
        mass = wt * total_mass  # step 4
        out["components"][c.name] = {
            "mole_ratio": c.n / n_target,
            "wt_fraction": wt,
            "mass_g": mass,
            "mmol": mass / c.mw * 1000,
        }
    return out
