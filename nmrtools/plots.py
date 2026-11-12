"""Figures written as image files: one panel per ppm range, and spectra
stacked one above the other with optional magnified insets."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
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


COLOURS = (ACCENT, WARM, MUTED, INK)
FORMATS = (".png", ".svg", ".pdf")


def _check_format(out: str | Path) -> Path:
    out = Path(out)
    if out.suffix.lower() not in FORMATS:
        raise ValueError(f"unsupported figure format {out.suffix!r}; use {', '.join(FORMATS)}")
    return out


def _peak_height(spec: Spectrum, pos: float, half_width: float = 0.03) -> float:
    """Tallest point within half_width ppm of pos, or 0 outside the axis."""
    near = spec.window(pos - half_width, pos + half_width)
    return float(spec.y[near].max()) if near.any() else 0.0


def plot_panels(
    spectra: dict[str, Spectrum],
    ranges: list[tuple[float, float]],
    out: str | Path,
    labels: dict[float, str] | None = None,
    windows: list[tuple[float, float, str]] | tuple = (),
) -> Path:
    """Draw one panel per ppm range, each scaled to its own tallest peak.

    The format follows the suffix of `out` (.png, .svg or .pdf). With several
    spectra the traces are overlaid and keep their relative heights. `labels`
    maps a shift to a text placed above the first spectrum's peak there.
    `windows` are (lo, hi, label) integration windows, drawn as shaded bands.
    """
    out = _check_format(out)
    ranges = [(min(r), max(r)) for r in ranges]
    for lo, hi in ranges:
        for name, spec in spectra.items():
            if not spec.window(lo, hi).any():
                raise ValueError(f"{name} has no points between {lo} and {hi} ppm")
    _use_font()
    fig, axes = plt.subplots(
        1,
        len(ranges),
        figsize=(7.2, 3.2),
        dpi=300,
        squeeze=False,
        gridspec_kw={"width_ratios": [hi - lo for lo, hi in ranges]},
    )
    first = next(iter(spectra.values()))
    for ax, (lo, hi) in zip(axes[0], ranges):
        top = max(float(spec.y[spec.window(lo, hi)].max()) for spec in spectra.values())
        top = top if top > 0 else 1.0
        for (name, spec), colour in zip(spectra.items(), COLOURS * len(spectra)):
            m = spec.window(lo, hi)
            ax.plot(spec.ppm[m], spec.y[m] / top, color=colour, lw=0.7, label=name)
        ax.set_xlim(hi, lo)
        ax.set_ylim(-0.03, 1.25)
        ax.set_yticks([])
        for side in ("left", "right", "top"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="x", colors=INK, labelsize=9, length=3)
        ax.set_xlabel("δ (ppm)", color=INK, fontsize=9.5, labelpad=1)
        for pos, text in (labels or {}).items():
            if lo <= pos <= hi:
                height = min(_peak_height(first, pos) / top + 0.04, 1.15)
                ax.text(pos, height, text, color=INK, fontsize=8, ha="center", va="bottom")
        for wlo, whi, label in windows:
            if lo <= wlo and whi <= hi:
                ax.axvspan(wlo, whi, color=ACCENT, alpha=0.12, lw=0)
                ax.text((wlo + whi) / 2, 1.17, label, color=MUTED, fontsize=7.5, ha="center", va="bottom")
    if len(spectra) > 1:
        axes[0][0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out, dpi=300, facecolor="white")
    plt.close(fig)
    return out


def _quietest_start(traces: list[Spectrum], lo: float, hi: float, width: float) -> float:
    """Left edge (axes fraction) of the band of `width` with the least signal in any trace."""
    best, best_height = 0.14, np.inf
    for start in np.linspace(0.14, 0.98 - width, 60):
        left, right = hi - start * (hi - lo), hi - (start + width) * (hi - lo)
        height = max(float(t.y[t.window(right, left)].max()) for t in traces)
        if height < best_height - 1e-9:
            best, best_height = float(start), height
    return best


def plot_stack(
    spectra: dict[str, Spectrum],
    lo: float,
    hi: float,
    out: str | Path,
    norm: tuple[float, float] | None = None,
    zoom: tuple[float, float] | None = None,
    zoom_factor: float = 10,
    labels: dict[str, dict[float, str]] | None = None,
) -> Path:
    """Stack spectra one above the other at true relative scale.

    Each trace is scaled to its own tallest point in `norm` (default: the
    plotted range), so traces are compared relative to that signal. `zoom`
    adds, for every trace, a separate labelled inset of that range magnified
    `zoom_factor` times; the inset is placed over the quietest part of the
    traces. `labels` maps a spectrum name to {shift: text}; a label inside
    the zoom range is written in the inset, any other one above its peak.
    """
    out = _check_format(out)
    lo, hi = sorted((lo, hi))
    labels = labels or {}
    step = 1.25
    traces: dict[str, Spectrum] = {}
    for name, spec in spectra.items():
        ref = spec.window(*(norm or (lo, hi)))
        if not spec.window(lo, hi).any() or not ref.any() or spec.y[ref].max() <= 0:
            raise ValueError(f"{name} has no positive signal to scale to in the requested range")
        traces[name] = spec.normalised(*(norm or (lo, hi)))
    n = len(traces)
    _use_font()

    height_in = 1.05 + step * n
    fig = plt.figure(figsize=(5.3, height_in), dpi=300)
    bottom = 0.46 / height_in
    ax = fig.add_axes([0.02, bottom, 0.96, 0.98 - bottom])
    y_lo, y_hi = -0.08, step * n
    offsets = {name: (n - 1 - i) * step for i, name in enumerate(traces)}
    for (name, spec), colour in zip(traces.items(), COLOURS * n):
        m = spec.window(lo, hi)
        off = offsets[name]
        ax.plot(spec.ppm[m], np.clip(spec.y[m], -0.02, 1.12) + off, color=colour, lw=0.7)
        ax.text(hi - 0.008 * (hi - lo), off + 0.16, name, color=colour, fontsize=10.5, fontweight="bold", va="bottom")
        for pos, text in labels.get(name, {}).items():
            if lo <= pos <= hi and not (zoom and min(zoom) <= pos <= max(zoom)):
                top = min(_peak_height(spec, pos) + 0.04, 1.14)
                ax.text(pos, off + top, text, color=colour, fontsize=8, ha="center", va="bottom")
    ax.set_xlim(hi, lo)
    ax.set_ylim(y_lo, y_hi)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="x", colors=INK, labelsize=9, length=3)
    ax.set_xlabel("δ (ppm)", color=INK, fontsize=9.5, labelpad=1)

    if zoom:
        zlo, zhi = sorted(zoom)
        width = 0.29
        start = _quietest_start(list(traces.values()), lo, hi, width)
        for (name, spec), colour in zip(traces.items(), COLOURS * n):
            y0 = (offsets[name] + 0.16 - y_lo) / (y_hi - y_lo)
            ia = ax.inset_axes([start, y0, width, 0.64 / (y_hi - y_lo)])
            ia.set_facecolor(PAPER)
            m = spec.window(zlo, zhi)
            if m.any():
                ia.plot(spec.ppm[m], spec.y[m] * zoom_factor, color=colour, lw=0.7)
            ia.set_xlim(zhi, zlo)
            ia.set_ylim(-0.04, 0.50)
            ia.set_yticks([])
            ia.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
            ia.tick_params(axis="x", labelsize=7, length=2, colors=MUTED, pad=1)
            for spine in ia.spines.values():
                spine.set_color(MUTED)
                spine.set_linewidth(0.6)
            ia.set_title(
                f"zoom of {zlo:g}–{zhi:g} ppm, ×{zoom_factor:g} vertical", fontsize=7, color=MUTED, loc="left", pad=2
            )
            for pos, text in labels.get(name, {}).items():
                if zlo <= pos <= zhi:
                    top = min(_peak_height(spec, pos) * zoom_factor + 0.03, 0.42)
                    ia.text(pos, top, text, color=colour, fontsize=6.8, ha="center", va="bottom")

    fig.savefig(out, dpi=300, facecolor="white")
    plt.close(fig)
    return out
