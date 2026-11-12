import numpy as np
import pytest

from nmrtools import Component, coupling_constants, integrate, load_spectrum, mass_correction, pick_lines
from nmrtools.bruker import Spectrum

SF = 400.23


def lorentz(x, x0, fwhm, area):
    hw = fwhm / 2
    return area / np.pi * hw / ((x - x0) ** 2 + hw**2)


def synthetic(peaks, sf=SF, n=65536, lo=-1.0, hi=11.0):
    """peaks: list of (ppm, width in Hz, area)."""
    ppm = np.linspace(hi, lo, n)
    y = sum(lorentz(ppm, p, w / sf, a) for p, w, a in peaks)
    return Spectrum(ppm=ppm, y=y, sf=sf)


def test_handout_worked_example():
    """Benzil with residual ethyl acetate, from the course handout."""
    res = mass_correction(
        [Component("benzil", 4.00, 4, 210.23), Component("EtOAc", 0.65, 2, 88.11)],
        total_mass=0.520,
    )
    benzil = res["components"]["benzil"]
    assert res["components"]["EtOAc"]["mole_ratio"] == pytest.approx(0.325)
    assert benzil["wt_fraction"] * 100 == pytest.approx(88.0, abs=0.05)
    assert benzil["mass_g"] == pytest.approx(0.458, abs=0.0005)


def test_three_component_expt4_numbers():
    res = mass_correction(
        [
            Component("linalool", 1.00, 1, 154.25),
            Component("side product", 0.57, 1, 128.21),
            Component("EtOAc", 0.16, 2, 88.11),
        ],
        total_mass=0.4725,
    )
    lin = res["components"]["linalool"]
    assert lin["wt_fraction"] * 100 == pytest.approx(65.8, abs=0.1)
    assert lin["mass_g"] == pytest.approx(0.311, abs=0.001)


def test_coupling_constants_dd():
    c, (j1, j2) = 5.912, (17.3, 10.7)
    lines = [c + s1 * j1 / 2 / SF + s2 * j2 / 2 / SF for s1 in (1, -1) for s2 in (1, -1)]
    centre, js = coupling_constants(lines, SF)
    assert centre == pytest.approx(c)
    assert js == pytest.approx((17.3, 10.7))


def test_reference_integrate_and_pick():
    spec = synthetic([(7.284, 1.0, 1.0), (1.30, 1.0, 3.0), (5.90, 1.0, 1.0)])
    spec = spec.referenced()  # CHCl3 found at 7.284 -> 7.26
    assert spec.meta["reference_shift"] == pytest.approx(-0.024, abs=2e-4)
    ratio = integrate(spec, 1.176, 1.376) / integrate(spec, 5.776, 5.976)
    assert ratio == pytest.approx(3.0, rel=0.01)
    assert pick_lines(spec, 1.2, 1.35, n=1)[0] == pytest.approx(1.276, abs=5e-4)


def test_baseline_poly_removes_offset_and_slope():
    spec = synthetic([(1.28, 1.0, 3.0)])
    tilted = Spectrum(ppm=spec.ppm, y=spec.y + 0.5 + 0.02 * spec.ppm, sf=SF)
    flat = tilted.baseline_poly()
    assert abs(np.median(flat.y[flat.window(8.5, 9.5)])) < 1e-3


def test_load_bruker_folder(tmp_path):
    si, sf, offset, sw_hz = 1024, 400.13, 12.0, 6402.08
    y = (np.arange(si) % 97).astype("<i4")
    pdata = tmp_path / "sample" / "1" / "pdata" / "1"
    pdata.mkdir(parents=True)
    y.tofile(pdata / "1r")
    (pdata / "procs").write_text(
        f"##$BYTORDP= 0\n##$NC_proc= -2\n##$OFFSET= {offset}\n##$SF= {sf}\n##$SI= {si}\n##$SW_p= {sw_hz}\n"
    )
    (pdata.parent.parent / "acqus").write_text(
        "##$PULPROG= <zg30>\n##$NS= 8\n##$D= (0..63)\n0 1 0 0\n##$TD= 32768\n##$SW_h= 8012.82\n##$SFO1= 400.2324\n"
    )
    spec = load_spectrum(tmp_path / "sample")
    assert spec.ppm[0] == pytest.approx(offset)
    assert spec.ppm[0] - spec.ppm[-1] == pytest.approx(sw_hz / sf * (si - 1) / si)
    assert spec.y[5] == pytest.approx(5 * 2.0**-2)
    assert spec.meta["pulprog"] == "zg30" and spec.meta["ns"] == 8 and spec.meta["d1_s"] == 1.0


def _write_experiment(root, expno, scale=1, acqus=None):
    si, sf = 1024, 400.13
    pdata = root / str(expno) / "pdata" / "1"
    pdata.mkdir(parents=True)
    (scale * (np.arange(si) % 97)).astype("<i4").tofile(pdata / "1r")
    (pdata / "procs").write_text(
        f"##$BYTORDP= 0\n##$NC_proc= 0\n##$OFFSET= 12.0\n##$SF= {sf}\n##$SI= {si}\n##$SW_p= 6402.08\n"
    )
    if acqus:
        (root / str(expno) / "acqus").write_text(acqus)


def test_mass_correction_rejects_unknown_target():
    comps = [Component("a", 1.0, 1, 100.0)]
    with pytest.raises(ValueError, match="target 'b'"):
        mass_correction(comps, 1.0, target="b")
    with pytest.raises(ValueError, match="no components"):
        mass_correction([], 1.0)


def test_jcamp_arrays_and_strings_wrap_over_lines(tmp_path):
    acqus = "##$D= (0..7)\n0 2.5 0 0\n0 0 0 7\n##$TI= <first\nsecond>\n$$ a comment\n##$NS= 16\n"
    _write_experiment(tmp_path / "s", 1, acqus=acqus)
    from nmrtools.bruker import read_jcamp

    params = read_jcamp(tmp_path / "s" / "1" / "acqus")
    assert params["D"].split() == ["0", "2.5", "0", "0", "0", "0", "0", "7"]
    assert params["TI"] == "<first second>"
    assert params["NS"] == "16"
    assert load_spectrum(tmp_path / "s").meta["d1_s"] == 2.5


def test_several_experiments_warn_and_expno_selects(tmp_path):
    _write_experiment(tmp_path / "s", 1)
    _write_experiment(tmp_path / "s", 2, scale=3)
    with pytest.warns(UserWarning, match="2 experiments found"):
        first = load_spectrum(tmp_path / "s")
    assert first.meta["expno"] == "1"
    second = load_spectrum(tmp_path / "s", expno=2)
    assert second.meta["expno"] == "2" and second.y[5] == pytest.approx(15)
    with pytest.raises(FileNotFoundError, match="found 1, 2"):
        load_spectrum(tmp_path / "s", expno=9)


def test_reference_uses_interpolated_line_position():
    # A line centred between two grid points: the tallest point is half a
    # point off, the refined position is not.
    spec = synthetic([(7.0, 1.0, 1.0)])
    step = spec.ppm[0] - spec.ppm[1]
    centre = float(spec.ppm[np.argmin(abs(spec.ppm - 7.28))] - 0.4 * step)
    spec = synthetic([(centre, 1.0, 1.0)])
    ref = spec.referenced()
    assert abs(ref.meta["reference_found"] - centre) < 0.1 * step
    assert pick_lines(ref, 7.2, 7.35, n=1)[0] == pytest.approx(7.26, abs=1e-9)
