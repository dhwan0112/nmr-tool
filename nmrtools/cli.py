"""Interactive command-line shell for looking at processed Bruker 1H spectra.

Start it with `python -m nmrtools [spectrum ...]`. Type `help` for the commands.
Figures are written as image files (PNG, SVG or PDF); nothing is drawn in the
terminal.
"""

from __future__ import annotations

import argparse
import cmd
import json
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import IO

from .analysis import Component, coupling_constants, integrate, mass_correction, pick_lines
from .bruker import Spectrum, load_spectrum

_KEYWORD = re.compile(r"^([a-z_]+)=(.*)$")

SUMMARY = """\
Spectra
  load <path> [name]            load a Bruker folder or .zip (no corrections applied)
  list                          loaded spectra
  use <name>                    make a spectrum the current one
  info                          acquisition parameters and axis of the current spectrum
Corrections
  ref [target [lo hi]]          put the tallest point in lo-hi at target (default 7.26 in 7.20-7.35)
  baseline offset [lo hi]       subtract the median of lo-hi (default 8.5-10.0)
  baseline poly [order] [lo:hi ...]  subtract a polynomial through signal-free regions
  prepare                       ref, baseline offset, baseline poly with the defaults
Reading the spectrum
  peaks <lo> <hi> [n=N] [min=0.15]   line positions in ppm and Hz
  j <lo> <hi> [n=4]             centre and coupling constants of a doublet (n=2) or dd (n=4)
Integrals and mass correction
  integrate <lo> <hi> [name=X] [nh=1] [mw=M]   add an integral to the table
  integrals                     table with I, I/nH and ratios to the first row
  remove <name>                 remove one integral
  clear                         remove all integrals of the current spectrum
  preset <file.json>            replace the table with the integrals listed in a preset file
  mass <total g> [target=X]     mass correction over the integrals that have mw
Output
  label <ppm> <text>            name a peak of the current spectrum for the figures
  labels                        the peak labels of the current spectrum
  unlabel <ppm>|all             remove one label, or all of them
  plot [lo hi ...] [out=file.png|.svg|.pdf] [with=name,...] [open]
                                one panel per lo-hi pair, each scaled to its tallest peak
  stack [name ...] [lo hi] [norm=lo:hi] [zoom=lo:hi] [x=10] [out=file] [open]
                                spectra one above the other, with magnified insets of zoom
  export <file.json>            write the corrections, labels and integral tables as JSON
  help [command]                this summary, or the usage of one command
  quit                          leave (also Ctrl-D)
"""


class CliError(Exception):
    """A user-facing error: wrong arguments or nothing loaded."""


@dataclass
class Integral:
    """One row of the integral table."""

    name: str
    lo: float
    hi: float
    n_h: int = 1
    mw: float | None = None


def split_args(line: str) -> tuple[list[str], dict[str, str]]:
    """Split a command line into positional words and key=value options."""
    words: list[str] = []
    options: dict[str, str] = {}
    for token in shlex.split(line):
        m = _KEYWORD.match(token)
        if m:
            options[m.group(1)] = m.group(2)
        else:
            words.append(token)
    return words, options


def _pair(text: str, usage: str) -> tuple[float, float]:
    """Parse a lo:hi option value."""
    parts = text.split(":")
    if len(parts) != 2:
        raise CliError(f"expected lo:hi, got {text!r}. usage: {usage}")
    lo, hi = sorted(_floats(parts, usage))
    return lo, hi


def _floats(words: list[str], usage: str) -> list[float]:
    try:
        return [float(w) for w in words]
    except ValueError:
        raise CliError(f"expected numbers. usage: {usage}") from None


class NmrShell(cmd.Cmd):
    """The interactive shell. Every `do_<name>` method is one command."""

    intro = "nmrtools shell. Type help for commands, quit to leave."
    prompt = "nmr> "

    def __init__(self, stdin: IO[str] | None = None, stdout: IO[str] | None = None) -> None:
        super().__init__(stdin=stdin, stdout=stdout)
        self.spectra: dict[str, Spectrum] = {}
        self.integrals: dict[str, list[Integral]] = {}
        self.labels: dict[str, dict[float, str]] = {}
        self.current: str | None = None
        self.errors = 0

    # -- plumbing ---------------------------------------------------------
    def say(self, text: str = "") -> None:
        """Print one line to the shell's output."""
        print(text, file=self.stdout)

    def onecmd(self, line: str) -> bool:
        """Run one command; report errors instead of leaving the shell."""
        try:
            return bool(super().onecmd(line))
        except CliError as exc:
            self.say(f"error: {exc}")
        except (ValueError, KeyError, IndexError, OSError, TypeError) as exc:
            self.say(f"error: {type(exc).__name__}: {exc}")
        self.errors += 1
        return False

    def emptyline(self) -> bool:
        """Do nothing on an empty line (the default would repeat the last command)."""
        return False

    def default(self, line: str) -> None:
        """Report an unknown command."""
        raise CliError(f"unknown command {line.split()[0]!r}. Type help.")

    def _spec(self) -> Spectrum:
        if self.current is None:
            raise CliError("no spectrum loaded. usage: load <path> [name]")
        return self.spectra[self.current]

    def _replace(self, spec: Spectrum) -> None:
        self.spectra[self.current] = spec

    def _rows(self) -> list[Integral]:
        self._spec()
        return self.integrals[self.current]

    # -- spectra ----------------------------------------------------------
    def do_load(self, line: str) -> None:
        """load <path> [name]  -- load a Bruker experiment folder or a .zip of one."""
        words, _ = split_args(line)
        if not 1 <= len(words) <= 2:
            raise CliError("usage: load <path> [name]")
        path = Path(words[0]).expanduser()
        if not path.exists():
            raise CliError(f"no such file or folder: {path}")
        name = words[1] if len(words) == 2 else f"s{len(self.spectra) + 1}"
        self.spectra[name] = load_spectrum(path)
        self.integrals[name] = []
        self.labels[name] = {}
        self.current = name
        spec = self.spectra[name]
        self.say(
            f"{name}: {len(spec.ppm)} points, {spec.ppm.max():.2f} to {spec.ppm.min():.2f} ppm, SF {spec.sf:.2f} MHz"
        )

    def do_list(self, line: str) -> None:
        """list  -- loaded spectra; * marks the current one."""
        if not self.spectra:
            self.say("no spectra loaded")
        for name, spec in self.spectra.items():
            mark = "*" if name == self.current else " "
            self.say(f"{mark} {name:10s} {spec.meta.get('source', '')}")

    def do_use(self, line: str) -> None:
        """use <name>  -- make a loaded spectrum the current one."""
        name = line.strip()
        if name not in self.spectra:
            raise CliError(f"no spectrum named {name!r}. loaded: {', '.join(self.spectra) or 'none'}")
        self.current = name

    def do_info(self, line: str) -> None:
        """info  -- acquisition parameters and axis of the current spectrum."""
        spec = self._spec()
        self.say(f"name          {self.current}")
        self.say(f"axis          {spec.ppm.max():.3f} to {spec.ppm.min():.3f} ppm, {len(spec.ppm)} points")
        self.say(f"resolution    {spec.hz_per_point:.3f} Hz/point at SF {spec.sf:.4f} MHz")
        for key, value in spec.meta.items():
            self.say(f"{key:13s} {value}")

    # -- corrections ------------------------------------------------------
    def do_ref(self, line: str) -> None:
        """ref [target [lo hi]]  -- shift the axis so the tallest point in lo-hi sits at target."""
        words, _ = split_args(line)
        nums = _floats(words, "ref [target [lo hi]]")
        if len(nums) not in (0, 1, 3):
            raise CliError("usage: ref [target [lo hi]]")
        target = nums[0] if nums else 7.26
        search = (nums[1], nums[2]) if len(nums) == 3 else (7.20, 7.35)
        spec = self._spec()
        if not spec.window(*search).any():
            raise CliError(f"no points between {search[0]} and {search[1]} ppm")
        spec = spec.referenced(target, search)
        self._replace(spec)
        self.say(
            f"reference found at {spec.meta['reference_found']:.4f}, shifted by {spec.meta['reference_shift']:+.4f} ppm"
        )

    def do_baseline(self, line: str) -> None:
        """baseline offset [lo hi] | baseline poly [order] [lo:hi ...]  -- flatten the baseline."""
        words, _ = split_args(line)
        usage = "baseline offset [lo hi] | baseline poly [order] [lo:hi ...]"
        if not words or words[0] not in ("offset", "poly"):
            raise CliError(f"usage: {usage}")
        spec = self._spec()
        if words[0] == "offset":
            nums = _floats(words[1:], usage)
            if len(nums) not in (0, 2):
                raise CliError(f"usage: {usage}")
            region = (nums[0], nums[1]) if nums else (8.5, 10.0)
            if not spec.window(*region).any():
                raise CliError(f"no points between {region[0]} and {region[1]} ppm")
            self._replace(spec.baseline_offset(region))
            self.say(f"subtracted the median of {region[0]}-{region[1]} ppm")
            return
        rest = words[1:]
        order = 2
        if rest and ":" not in rest[0]:
            order = int(_floats(rest[:1], usage)[0])
            rest = rest[1:]
        regions = tuple(tuple(_floats(r.split(":"), usage)) for r in rest)
        if any(len(r) != 2 for r in regions):
            raise CliError(f"regions are written lo:hi. usage: {usage}")
        if regions:
            self._replace(spec.baseline_poly(regions, order))
            self.say(f"subtracted an order-{order} polynomial through {len(regions)} regions")
        else:
            self._replace(spec.baseline_poly(order=order))
            self.say(f"subtracted an order-{order} polynomial through the default regions")

    def do_prepare(self, line: str) -> None:
        """prepare  -- ref, baseline offset and baseline poly with their defaults."""
        self.do_ref("")
        self.do_baseline("offset")
        self.do_baseline("poly")

    # -- reading ----------------------------------------------------------
    def _lines(self, line: str, usage: str, default_n: int | None) -> tuple[list[float], Spectrum]:
        words, options = split_args(line)
        nums = _floats(words, usage)
        if len(nums) != 2:
            raise CliError(f"usage: {usage}")
        spec = self._spec()
        if not spec.window(nums[0], nums[1]).any():
            raise CliError(f"no points between {nums[0]} and {nums[1]} ppm")
        n = int(options["n"]) if "n" in options else default_n
        min_rel = float(options.get("min", 0.15))
        return pick_lines(spec, nums[0], nums[1], n=n, min_rel_height=min_rel), spec

    def do_peaks(self, line: str) -> None:
        """peaks <lo> <hi> [n=N] [min=0.15]  -- line positions; n keeps the N tallest, min is relative height."""
        lines, spec = self._lines(line, "peaks <lo> <hi> [n=N] [min=0.15]", None)
        if not lines:
            self.say("no lines found")
            return
        self.say("   #   ppm       Hz        gap to next (Hz)")
        for i, pos in enumerate(lines):
            gap = f"{(pos - lines[i + 1]) * spec.sf:8.2f}" if i + 1 < len(lines) else ""
            self.say(f"  {i + 1:2d}  {pos:7.4f}  {pos * spec.sf:9.2f}  {gap}")

    def do_j(self, line: str) -> None:
        """j <lo> <hi> [n=4]  -- centre and J of a first-order doublet (n=2) or doublet of doublets (n=4)."""
        lines, spec = self._lines(line, "j <lo> <hi> [n=4]", 4)
        if len(lines) not in (2, 4):
            raise CliError(f"found {len(lines)} lines; need 2 (n=2) or 4. Check the window with peaks.")
        centre, js = coupling_constants(lines, spec.sf)
        mult = "d" if len(js) == 1 else "dd"
        self.say(f"δ {centre:.3f} {mult}, J = {', '.join(f'{j:.1f}' for j in js)} Hz")

    # -- integrals --------------------------------------------------------
    def do_integrate(self, line: str) -> None:
        """integrate <lo> <hi> [name=X] [nh=1] [mw=M]  -- add an integral; mw (g/mol) enables `mass`."""
        usage = "integrate <lo> <hi> [name=X] [nh=1] [mw=M]"
        words, options = split_args(line)
        nums = _floats(words, usage)
        if len(nums) != 2:
            raise CliError(f"usage: {usage}")
        rows = self._rows()
        lo, hi = sorted(nums)
        if not self._spec().window(lo, hi).any():
            raise CliError(f"no points between {lo} and {hi} ppm")
        name = options.get("name", f"i{len(rows) + 1}")
        row = Integral(name, lo, hi, int(options.get("nh", 1)), float(options["mw"]) if "mw" in options else None)
        if row.n_h < 1:
            raise CliError("nh must be at least 1")
        rows[:] = [r for r in rows if r.name != name] + [row]
        self.do_integrals("")

    def do_integrals(self, line: str) -> None:
        """integrals  -- the integral table of the current spectrum."""
        rows = self._rows()
        if not rows:
            self.say("no integrals. usage: integrate <lo> <hi> [name=X] [nh=1] [mw=M]")
            return
        spec = self._spec()
        values = [integrate(spec, r.lo, r.hi) for r in rows]
        first = values[0] / rows[0].n_h
        self.say(f"  {'name':14s} {'window (ppm)':15s} {'nH':>2s} {'integral':>12s} {'ratio':>7s}  mw")
        for r, v in zip(rows, values):
            ratio = f"{v / r.n_h / first:7.3f}" if first else "    n/a"
            mw = f"{r.mw:g}" if r.mw is not None else "-"
            self.say(f"  {r.name:14s} {r.lo:6.3f}-{r.hi:<8.3f} {r.n_h:2d} {v:12.4g} {ratio}  {mw}")
        self.say("ratio = (integral / nH) relative to the first row")

    def do_remove(self, line: str) -> None:
        """remove <name>  -- remove one integral from the table."""
        rows = self._rows()
        name = line.strip()
        if name not in [r.name for r in rows]:
            raise CliError(f"no integral named {name!r}")
        rows[:] = [r for r in rows if r.name != name]

    def do_clear(self, line: str) -> None:
        """clear  -- remove all integrals of the current spectrum."""
        self._rows().clear()

    def do_preset(self, line: str) -> None:
        """preset <file.json>  -- replace the integral table with the integrals listed in a preset file.

        The file holds {"integrals": [{"name": ..., "window_ppm": [lo, hi],
        "n_H": 1, "mw": 154.25}, ...]}; n_H and mw are optional. See
        examples/cm3191/expt4.json.
        """
        words, _ = split_args(line)
        if len(words) != 1:
            raise CliError("usage: preset <file.json>")
        path = Path(words[0]).expanduser()
        if not path.is_file():
            raise CliError(f"no such preset file: {path}")
        rows = self._rows()
        try:
            entries = json.loads(path.read_text(encoding="utf-8"))["integrals"]
            new = [
                Integral(
                    str(e["name"]),
                    *sorted(float(v) for v in e["window_ppm"]),
                    n_h=int(e.get("n_H", 1)),
                    mw=None if e.get("mw") is None else float(e["mw"]),
                )
                for e in entries
            ]
        except (ValueError, KeyError, TypeError) as exc:
            raise CliError(f"{path} is not a preset file: {type(exc).__name__}: {exc}") from None
        rows[:] = new
        self.do_integrals("")

    def _mass(self, total_mass: float, target: str | None) -> dict:
        spec = self._spec()
        rows = [r for r in self._rows() if r.mw is not None]
        if not rows:
            raise CliError("no integrals with mw. add them with integrate ... mw=M")
        if target is not None and target not in [r.name for r in rows]:
            raise CliError(f"target {target!r} is not an integral with mw")
        comps = [Component(r.name, max(integrate(spec, r.lo, r.hi), 0.0), r.n_h, r.mw) for r in rows]
        return mass_correction(comps, total_mass, target=target)

    def do_mass(self, line: str) -> None:
        """mass <total g> [target=X]  -- split a weighed mass over the integrals that have mw."""
        words, options = split_args(line)
        nums = _floats(words, "mass <total g> [target=X]")
        if len(nums) != 1:
            raise CliError("usage: mass <total g> [target=X]")
        res = self._mass(nums[0], options.get("target"))
        skipped = [r.name for r in self._rows() if r.mw is None]
        self.say(f"  {'name':14s} {'mole ratio':>10s} {'wt %':>7s} {'mass (g)':>9s} {'mmol':>7s}")
        for name, c in res["components"].items():
            self.say(
                f"  {name:14s} {c['mole_ratio']:10.3f} {c['wt_fraction'] * 100:7.1f} {c['mass_g']:9.4f} {c['mmol']:7.3f}"
            )
        self.say(f"mole ratios are relative to {res['target']}; negative integrals count as zero")
        if skipped:
            self.say(f"not counted (no mw): {', '.join(skipped)}")
        self.say("only the listed species are counted, so each weight fraction is an upper bound")

    # -- output -----------------------------------------------------------
    def do_label(self, line: str) -> None:
        """label <ppm> <text>  -- name a peak of the current spectrum; plot and stack write the text above it."""
        words, _ = split_args(line)
        if len(words) < 2:
            raise CliError("usage: label <ppm> <text>")
        (pos,) = _floats(words[:1], "label <ppm> <text>")
        spec = self._spec()
        if not float(spec.ppm.min()) <= pos <= float(spec.ppm.max()):
            raise CliError(f"{pos:g} ppm is outside the spectrum")
        self.labels[self.current][pos] = " ".join(words[1:])
        self.say(f"{self.current}: {pos:.3f} ppm = {self.labels[self.current][pos]}")

    def do_labels(self, line: str) -> None:
        """labels  -- the peak labels of the current spectrum."""
        self._spec()
        rows = self.labels[self.current]
        if not rows:
            self.say("no labels. usage: label <ppm> <text>")
        for pos in sorted(rows, reverse=True):
            self.say(f"  {pos:7.3f}  {rows[pos]}")

    def do_unlabel(self, line: str) -> None:
        """unlabel <ppm>|all  -- remove the label at that shift, or every label of the current spectrum."""
        words, _ = split_args(line)
        self._spec()
        rows = self.labels[self.current]
        if words == ["all"]:
            rows.clear()
            return
        if len(words) != 1:
            raise CliError("usage: unlabel <ppm>|all")
        (pos,) = _floats(words, "unlabel <ppm>|all")
        hits = [p for p in rows if abs(p - pos) < 5e-4]
        if not hits:
            raise CliError(f"no label at {pos:g} ppm")
        for p in hits:
            del rows[p]

    def _open(self, out: Path) -> None:
        opener = shutil.which("xdg-open") or shutil.which("open")
        if opener is None:
            raise CliError("no xdg-open or open on this system; open the file yourself")
        subprocess.Popen([opener, str(out)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def do_plot(self, line: str) -> None:
        """plot [lo hi ...] [out=file.png|.svg|.pdf] [with=name,...] [open]  -- write a figure of the current spectrum.

        Every lo-hi pair becomes one panel, scaled to the tallest peak in it,
        so `plot 4.9 6.1 1.1 2.3` shows two regions side by side. Peak labels
        are written above their peaks and integral windows are shaded.
        `with=` overlays other loaded spectra on the same absolute scale.
        `open` hands the file to the desktop's default viewer.
        """
        from .plots import plot_panels  # matplotlib is slow to import; only load it when drawing

        usage = "plot [lo hi ...] [out=file] [with=name,...] [open]"
        words, options = split_args(line)
        open_after = "open" in words
        nums = _floats([w for w in words if w != "open"], usage)
        if len(nums) % 2:
            raise CliError(f"usage: {usage}")
        spec = self._spec()
        ranges = [(min(a, b), max(a, b)) for a, b in zip(nums[::2], nums[1::2])]
        if not ranges:
            ranges = [(max(float(spec.ppm.min()), -0.5), min(float(spec.ppm.max()), 10.5))]
        names = [self.current] + [n for n in options.get("with", "").split(",") if n]
        missing = [n for n in names if n not in self.spectra]
        if missing:
            raise CliError(f"no spectrum named {', '.join(missing)}")
        tag = "_".join(f"{hi:g}-{lo:g}" for lo, hi in ranges)
        out = Path(options.get("out", f"nmr_{self.current}_{tag}.png")).expanduser()
        windows = [(r.lo, r.hi, r.name) for r in self._rows()]
        plot_panels({n: self.spectra[n] for n in names}, ranges, out, self.labels[self.current], windows)
        self.say(f"wrote {out}")
        if open_after:
            self._open(out)

    def do_stack(self, line: str) -> None:
        """stack [name ...] [lo hi] [norm=lo:hi] [zoom=lo:hi] [x=10] [out=file] [open]  -- compare spectra in one figure.

        The named spectra (default: all loaded, in load order) are drawn one
        above the other. Each is scaled to its tallest point in `norm`
        (default: the plotted range), so heights are relative to that signal.
        `zoom=` adds a labelled inset of that range for every spectrum,
        magnified `x` times. Peak labels of each spectrum are drawn too.
        """
        from .plots import plot_stack

        usage = "stack [name ...] [lo hi] [norm=lo:hi] [zoom=lo:hi] [x=10] [out=file] [open]"
        words, options = split_args(line)
        open_after = "open" in words
        words = [w for w in words if w != "open"]
        names = [w for w in words if w in self.spectra] or list(self.spectra)
        nums = _floats([w for w in words if w not in self.spectra], usage)
        if len(nums) not in (0, 2):
            raise CliError(f"usage: {usage}")
        if not names:
            raise CliError("no spectrum loaded. usage: load <path> [name]")
        if nums:
            lo, hi = sorted(nums)
        else:
            lo = max(max(float(self.spectra[n].ppm.min()) for n in names), -0.5)
            hi = min(min(float(self.spectra[n].ppm.max()) for n in names), 10.5)
        norm = _pair(options["norm"], usage) if "norm" in options else None
        zoom = _pair(options["zoom"], usage) if "zoom" in options else None
        (factor,) = _floats([options.get("x", "10")], usage)
        if factor <= 0:
            raise CliError("x must be positive")
        out = Path(options.get("out", f"nmr_stack_{'_'.join(names)}_{hi:g}-{lo:g}.png")).expanduser()
        plot_stack(
            {n: self.spectra[n] for n in names},
            lo,
            hi,
            out,
            norm=norm,
            zoom=zoom,
            zoom_factor=factor,
            labels={n: self.labels[n] for n in names},
        )
        self.say(f"wrote {out}")
        if open_after:
            self._open(out)

    def do_export(self, line: str) -> None:
        """export <file.json>  -- write the loaded spectra, their corrections and integral tables as JSON."""
        words, _ = split_args(line)
        if len(words) != 1:
            raise CliError("usage: export <file.json>")
        data = {}
        for name, spec in self.spectra.items():
            data[name] = {
                "meta": spec.meta,
                "sf_mhz": spec.sf,
                "labels": [{"ppm": pos, "text": text} for pos, text in sorted(self.labels[name].items())],
                "integrals": [
                    {
                        "name": r.name,
                        "window_ppm": [r.lo, r.hi],
                        "n_H": r.n_h,
                        "mw": r.mw,
                        "integral": integrate(spec, r.lo, r.hi),
                    }
                    for r in self.integrals[name]
                ],
            }
        out = Path(words[0]).expanduser()
        out.write_text(json.dumps(data, indent=2, ensure_ascii=False))
        self.say(f"wrote {out}")

    def do_help(self, arg: str) -> None:
        """help [command]  -- the command summary, or the usage of one command."""
        if arg:
            super().do_help(arg)
        else:
            self.say(SUMMARY.rstrip())

    def do_quit(self, line: str) -> bool:
        """quit  -- leave the shell."""
        return True

    do_exit = do_quit

    def do_EOF(self, line: str) -> bool:
        """Leave on Ctrl-D."""
        if self.stdin.isatty():
            self.say()
        return True


def main(argv: list[str] | None = None) -> int:
    """Entry point: load any spectra given on the command line, then run commands or start the shell."""
    parser = argparse.ArgumentParser(prog="nmrtools", description=__doc__.splitlines()[0])
    parser.add_argument("spectra", nargs="*", type=Path, help="Bruker folders or .zip files to load as s1, s2, ...")
    parser.add_argument(
        "-c",
        "--command",
        action="append",
        default=[],
        metavar="CMD",
        help="run CMD and exit instead of starting the shell; repeatable, and ';' separates commands",
    )
    args = parser.parse_args(argv)

    shell = NmrShell()
    for path in args.spectra:
        shell.onecmd(f"load {shlex.quote(str(path))}")

    if args.command:
        for chunk in args.command:
            for command in chunk.split(";"):
                if command.strip() and shell.onecmd(command.strip()):
                    return 1 if shell.errors else 0
        return 1 if shell.errors else 0

    interactive = sys.stdin.isatty()
    if interactive:
        try:
            import readline  # noqa: F401 - importing it gives the prompt history and line editing
        except ImportError:
            pass
    else:
        shell.prompt = ""
        shell.use_rawinput = False
    shell.cmdloop(shell.intro if interactive else "")
    return 1 if shell.errors and not interactive else 0


if __name__ == "__main__":
    raise SystemExit(main())
