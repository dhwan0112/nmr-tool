"""Read a processed Bruker 1D spectrum (pdata/<n>/1r) without TopSpin or MNova.

Only the already-processed real spectrum is read. Fourier transformation,
apodisation and phasing are whatever the spectrometer software applied.
"""

from __future__ import annotations

import re
import tempfile
import warnings
import zipfile
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

_PARAM = re.compile(r"^##\$(\w+)=\s*(.*)$")

# Signal-free regions of the linalool spectra in CDCl3 (ppm).
DEFAULT_BASELINE_REGIONS = ((8.3, 9.9), (6.2, 6.9), (4.4, 4.8), (2.8, 3.4), (0.15, 0.6))


def read_jcamp(path: str | Path) -> dict[str, str]:
    """Return the ##$KEY= value pairs of a Bruker parameter file as strings."""
    params: dict[str, str] = {}
    lines = Path(path).read_text(errors="ignore").splitlines()
    for i, line in enumerate(lines):
        m = _PARAM.match(line.strip())
        if not m:
            continue
        key, value = m.group(1), m.group(2).strip()
        is_array = value.startswith("(0..")
        if is_array or (value.startswith("<") and not value.endswith(">")):
            # Array values start on the next line and may wrap over several;
            # a <string> may also wrap. Collect up to the next record.
            rest = []
            for more in lines[i + 1 :]:
                if more.startswith(("##", "$$")):
                    break
                rest.append(more.strip())
            value = " ".join(rest) if is_array else " ".join([value, *rest])
        params[key] = value
    return params


@dataclass(frozen=True)
class Spectrum:
    """A 1D spectrum on a ppm axis (descending, as plotted)."""

    ppm: np.ndarray
    y: np.ndarray
    sf: float  # spectrometer frequency of the processed data, MHz
    meta: dict = field(default_factory=dict)

    @property
    def hz_per_point(self) -> float:
        """Digital resolution of the axis in Hz per point."""
        return float(abs(self.ppm[1] - self.ppm[0]) * self.sf)

    def window(self, lo: float, hi: float) -> np.ndarray:
        """Boolean mask for lo <= ppm <= hi."""
        lo, hi = sorted((lo, hi))
        return (self.ppm >= lo) & (self.ppm <= hi)

    def refined_position(self, i: int) -> float:
        """Sub-point position (ppm) of the maximum at index i, from a parabola
        through that point and its two neighbours."""
        if not 0 < i < len(self.y) - 1:
            return float(self.ppm[i])
        a, b, c = self.y[i - 1], self.y[i], self.y[i + 1]
        denom = a - 2 * b + c
        d = 0.0 if denom == 0 else 0.5 * (a - c) / denom
        return float(self.ppm[i] + d * (self.ppm[i + 1] - self.ppm[i]))

    # -- corrections; each returns a new Spectrum ---------------------------
    def referenced(self, target: float = 7.26, search: tuple[float, float] = (7.20, 7.35)) -> Spectrum:
        """Shift the axis so the tallest line in `search` sits at `target`.

        The line position is the tallest point refined with
        `refined_position`. Default: residual CHCl3 in CDCl3 at 7.26 ppm.
        """
        idx = np.where(self.window(*search))[0]
        found = self.refined_position(int(idx[np.argmax(self.y[idx])]))
        shift = target - found
        meta = dict(self.meta, reference_found=found, reference_shift=shift)
        return replace(self, ppm=self.ppm + shift, meta=meta)

    def baseline_offset(self, region: tuple[float, float] = (8.5, 10.0)) -> Spectrum:
        """Subtract the median of a signal-free region (constant offset)."""
        return replace(self, y=self.y - np.median(self.y[self.window(*region)]))

    def baseline_poly(
        self,
        regions: tuple[tuple[float, float], ...] = DEFAULT_BASELINE_REGIONS,
        order: int = 2,
    ) -> Spectrum:
        """Subtract a polynomial fitted through signal-free regions.

        The default regions are empty in the linalool spectra in CDCl3. Check
        them before using this on another compound.
        """
        m = np.zeros_like(self.ppm, dtype=bool)
        for lo, hi in regions:
            m |= self.window(lo, hi)
        coeff = np.polyfit(self.ppm[m], self.y[m], order)
        return replace(self, y=self.y - np.polyval(coeff, self.ppm))

    def normalised(self, lo: float, hi: float) -> Spectrum:
        """Scale so the tallest point between lo and hi equals 1."""
        return replace(self, y=self.y / self.y[self.window(lo, hi)].max())


def load_spectrum(path: str | Path, procno: int = 1, expno: int | str | None = None) -> Spectrum:
    """Load a processed spectrum from a Bruker experiment folder or a .zip of one.

    `expno` picks the experiment when the folder or archive holds several
    (the folder above pdata, e.g. 10 in sample/10/pdata/1/1r). Without it the
    first one in sorted order is loaded and a warning lists the others.
    """
    path = Path(path)
    if path.is_file() and path.suffix.lower() == ".zip":
        with tempfile.TemporaryDirectory() as tmp:
            with zipfile.ZipFile(path) as zf:
                zf.extractall(tmp)
            return _load_dir(Path(tmp), procno, expno, source=path.name)
    return _load_dir(path, procno, expno, source=path.name)


def _find_pdata(root: Path, procno: int, expno: int | str | None = None) -> Path:
    """The pdata/<procno> folder to read. The experiment number is the name of
    the folder two levels above it."""
    hits = sorted(p.parent for p in root.rglob("1r") if p.parent.name == str(procno))
    if not hits:
        raise FileNotFoundError(f"no pdata/{procno}/1r below {root}")
    found = ", ".join(h.parent.parent.name for h in hits)
    if expno is not None:
        chosen = [h for h in hits if h.parent.parent.name == str(expno)]
        if not chosen:
            raise FileNotFoundError(f"no experiment {expno} with pdata/{procno}/1r; found {found}")
        return chosen[0]
    if len(hits) > 1:
        warnings.warn(
            f"{len(hits)} experiments found ({found}); loading {hits[0].parent.parent.name}. "
            "Pass expno to choose another.",
            stacklevel=3,
        )
    return hits[0]


def _load_dir(root: Path, procno: int, expno: int | str | None, source: str) -> Spectrum:
    pdata = _find_pdata(root, procno, expno)
    procs = read_jcamp(pdata / "procs")
    acqus_path = pdata.parent.parent / "acqus"
    acqus = read_jcamp(acqus_path) if acqus_path.exists() else {}

    little = procs.get("BYTORDP", "0") == "0"
    if procs.get("DTYPP", "0") == "2":  # double-precision processed data
        dtype = "<f8" if little else ">f8"
    else:
        dtype = "<i4" if little else ">i4"
    y = np.fromfile(pdata / "1r", dtype=dtype).astype(float)
    y *= 2.0 ** float(procs.get("NC_proc", 0))

    si = int(procs["SI"])
    sf = float(procs["SF"])
    offset = float(procs["OFFSET"])  # ppm of the first point
    sw_ppm = float(procs["SW_p"]) / sf
    ppm = offset - np.arange(si) * sw_ppm / si

    d = acqus.get("D", "").split()
    meta = {
        "source": source,
        "expno": pdata.parent.parent.name,
        "pulprog": acqus.get("PULPROG", "").strip("<>"),
        "ns": int(acqus["NS"]) if "NS" in acqus else None,
        "d1_s": float(d[1]) if len(d) > 1 else None,
        "td": int(acqus["TD"]) if "TD" in acqus else None,
        "sw_hz": float(acqus["SW_h"]) if "SW_h" in acqus else None,
        "sfo1_mhz": float(acqus["SFO1"]) if "SFO1" in acqus else None,
        "solvent": acqus.get("SOLVENT", "").strip("<>"),
    }
    return Spectrum(ppm=ppm, y=y, sf=sf, meta=meta)
