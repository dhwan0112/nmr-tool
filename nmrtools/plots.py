"""Figures: the purified spectrum in two panels, and crude vs purified with
clearly labelled zoom insets."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager as fm

from .bruker import Spectrum

INK, ACCENT, WARM, MUTED, PAPER = "#1E2A32", "#2F6690", "#A4462F", "#5B6770", "#FAFAFA"


def _use_font(preferred: str = "IBM Plex Sans") -> None:
    """Use the preferred font when it is installed; otherwise keep the default."""
    available = {f.name for f in fm.fontManager.ttflist}
    plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.8})
    if preferred in available:
        plt.rcParams["font.family"] = preferred


def _clean_axis(ax, ticks) -> None:
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.set_xticks(ticks)
    ax.tick_params(axis="x", colors=INK, labelsize=9, length=3)


def plot_purified(spec: Spectrum, out: str | Path, labels: dict[float, str] | None = None) -> Path:
    """Purified spectrum: alkene region (6.1-4.9 ppm) and aliphatic region
    (2.3-1.1 ppm), each scaled to its own tallest peak."""
    _use_font()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), dpi=300, gridspec_kw={"width_ratios": [1, 1.2]})
    fig.patch.set_alpha(0)
    for ax, (hi, lo), ticks in zip(
        axes, ((6.1, 4.9), (2.3, 1.1)), ([6.0, 5.8, 5.6, 5.4, 5.2, 5.0], [2.2, 2.0, 1.8, 1.6, 1.4, 1.2])
    ):
        ax.patch.set_alpha(0)
        m = spec.window(lo, hi)
        y = spec.y[m] / spec.y[m].max()
        ax.plot(spec.ppm[m], y, color=ACCENT, lw=0.7)
        ax.set_xlim(hi, lo)
        ax.set_ylim(-0.03, 1.25)
        _clean_axis(ax, ticks)
        ax.set_xlabel("δ (ppm)", color=INK, fontsize=9.5, labelpad=1)
        for pos, text in (labels or {}).items():
            if lo <= pos <= hi:
                near = spec.window(pos - 0.03, pos + 0.03)
                top = spec.y[near].max() / spec.y[m].max()
                ax.text(pos, min(top + 0.04, 1.15), text, color=INK, fontsize=8, ha="center", va="bottom")
    fig.tight_layout()
    out = Path(out)
    fig.savefig(out, dpi=300, transparent=True)
    plt.close(fig)
    return out


def plot_crude_vs_purified(
    crude: Spectrum,
    purified: Spectrum,
    out: str | Path,
    norm_window: tuple[float, float] = (1.26, 1.30),
    zoom: tuple[float, float] = (3.60, 4.30),
    zoom_factor: float = 10,
) -> Path:
    """Stacked crude and purified spectra at true relative scale.

    Each trace is scaled to its own peak in `norm_window` (the linalool
    C(OH)CH3 singlet). The small signals in `zoom` are shown in separate,
    labelled inset boxes at `zoom_factor`, never magnified inside the trace.
    """
    _use_font()
    traces = {
        "Crude": (crude.normalised(*norm_window), WARM, 1.25),
        "Purified": (purified.normalised(*norm_window), ACCENT, 0.0),
    }
    lo, hi = 0.45, 6.2

    fig = plt.figure(figsize=(5.3, 3.55), dpi=300)
    fig.patch.set_alpha(0)
    ax = fig.add_axes([0.02, 0.13, 0.96, 0.85])
    ax.patch.set_alpha(0)
    for name, (spec, colour, offset) in traces.items():
        m = spec.window(lo, hi)
        ax.plot(spec.ppm[m], np.clip(spec.y[m], -0.02, 1.12) + offset, color=colour, lw=0.7)
        ax.text(6.15, offset + 0.16, name, color=colour, fontsize=10.5, fontweight="bold", ha="left", va="bottom")
    ax.set_xlim(hi, lo)
    ax.set_ylim(-0.08, 2.5)
    _clean_axis(ax, [6, 5, 4, 3, 2, 1])
    ax.set_xlabel("δ (ppm)", color=INK, fontsize=9.5, labelpad=1)

    # Signals that are visible at true scale are labelled on the crude trace.
    c, _, c_off = traces["Crude"]
    line = {"arrowstyle": "-", "color": WARM, "lw": 0.5, "shrinkA": 1, "shrinkB": 1}
    h_etoac = c.y[c.window(2.03, 2.06)].max()
    h_doublet = c.y[c.window(1.17, 1.21)].max()
    ax.annotate(
        "EtOAc 2.05 s",
        xy=(2.045, c_off + h_etoac + 0.02),
        xytext=(2.22, c_off + 0.80),
        color=WARM,
        fontsize=8,
        ha="center",
        arrowprops=line,
    )
    ax.annotate(
        "CH(OH)CH₃\n1.19 d",
        xy=(1.17, c_off + h_doublet * 0.9),
        xytext=(1.04, c_off + 0.50),
        color=WARM,
        fontsize=8,
        ha="left",
        va="bottom",
        linespacing=1.0,
        arrowprops=line,
    )

    zlo, zhi = zoom
    for (name, (spec, colour, _)), (bx, by) in zip(traces.items(), ((0.32, 0.62), (0.32, 0.20))):
        ia = fig.add_axes([bx, by, 0.28, 0.21])
        ia.set_facecolor(PAPER)
        m = spec.window(zlo - 0.02, zhi + 0.02)
        ia.plot(spec.ppm[m], spec.y[m] * zoom_factor, color=colour, lw=0.7)
        ia.set_xlim(zhi, zlo)
        ia.set_ylim(-0.04, 0.50)
        ia.set_yticks([])
        ia.set_xticks([4.2, 4.0, 3.8, 3.6])
        ia.tick_params(axis="x", labelsize=7, length=2, colors=MUTED, pad=1)
        for spine in ia.spines.values():
            spine.set_color(MUTED)
            spine.set_linewidth(0.6)
        ia.text(
            0.0,
            1.03,
            f"zoom of {zlo:.1f}–{zhi:.1f} ppm, ×{zoom_factor:g} vertical",
            transform=ia.transAxes,
            fontsize=7.3,
            color=MUTED,
            va="bottom",
            ha="left",
        )
        if name == "Crude":
            ia.text(4.13, 0.30, "EtOAc 4.12 q", fontsize=6.8, color=WARM, ha="center", va="bottom")
            ia.text(3.78, 0.30, "CH(OH) 3.81 m", fontsize=6.8, color=WARM, ha="center", va="bottom")
        else:
            ia.text(3.85, 0.11, "3.85 s, unassigned", fontsize=6.8, color=MUTED, ha="center", va="bottom")

    out = Path(out)
    fig.savefig(out, dpi=300, transparent=True)
    plt.close(fig)
    return out
