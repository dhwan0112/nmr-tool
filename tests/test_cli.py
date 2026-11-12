import io
import json
from pathlib import Path

import numpy as np
import pytest

from nmrtools.cli import NmrShell, main, split_args

EXPT4_PRESET = Path(__file__).resolve().parents[1] / "examples" / "cm3191" / "expt4.json"

SF = 400.13
SI = 65536
OFFSET = 12.0
SW_HZ = 6402.08


def lorentz(x, x0, fwhm, area):
    hw = fwhm / 2
    return area / np.pi * hw / ((x - x0) ** 2 + hw**2)


@pytest.fixture
def bruker_dir(tmp_path):
    """A synthetic processed spectrum on disk: CHCl3 off-reference at 7.284,
    a dd (17.3, 10.7 Hz) and two singlets with areas 1 and 2."""
    ppm = OFFSET - np.arange(SI) * (SW_HZ / SF) / SI
    shift = 0.024
    y = lorentz(ppm, 7.26 + shift, 1.0 / SF, 0.5)
    for s1 in (1, -1):
        for s2 in (1, -1):
            y += lorentz(ppm, 5.91 + shift + (s1 * 17.3 + s2 * 10.7) / 2 / SF, 1.0 / SF, 0.25)
    y += lorentz(ppm, 3.81 + shift, 1.0 / SF, 1.0)
    y += lorentz(ppm, 4.12 + shift, 1.0 / SF, 2.0)
    pdata = tmp_path / "sample" / "1" / "pdata" / "1"
    pdata.mkdir(parents=True)
    (y / y.max() * 2**28).astype("<i4").tofile(pdata / "1r")
    (pdata / "procs").write_text(
        f"##$BYTORDP= 0\n##$NC_proc= 0\n##$OFFSET= {OFFSET}\n##$SF= {SF}\n##$SI= {SI}\n##$SW_p= {SW_HZ}\n"
    )
    (pdata.parent.parent / "acqus").write_text("##$PULPROG= <zg30>\n##$NS= 8\n##$D= (0..63)\n0 1 0 0\n")
    return tmp_path / "sample"


def run(shell, *commands):
    """Run commands and return what they printed."""
    shell.stdout = io.StringIO()
    for command in commands:
        shell.onecmd(command)
    return shell.stdout.getvalue()


def test_split_args():
    assert split_args('5.8 6.0 name="side product" nh=2') == (["5.8", "6.0"], {"name": "side product", "nh": "2"})


def test_load_reference_and_coupling(bruker_dir):
    shell = NmrShell()
    out = run(shell, f"load {bruker_dir} demo", "ref", "j 5.86 5.96")
    assert "demo: 65536 points" in out
    assert "shifted by -0.02" in out
    assert "δ 5.910 dd, J = 17.3, 10.7 Hz" in out
    assert shell.errors == 0


def test_integrals_and_mass(bruker_dir):
    shell = NmrShell()
    out = run(
        shell,
        f"load {bruker_dir}",
        "ref",
        "integrate 5.80 6.00 name=A nh=1 mw=100",
        "integrate 3.70 3.92 name=B nh=1 mw=100",
        "integrate 4.05 4.20 name=C nh=2",
        "mass 1.0 target=A",
    )
    assert shell.errors == 0
    res = shell._mass(1.0, "A")["components"]
    assert res["B"]["mole_ratio"] == pytest.approx(1.0, rel=0.02)
    assert res["A"]["mass_g"] == pytest.approx(0.5, rel=0.02)
    assert "not counted (no mw): C" in out


def test_plot_and_export_write_files(bruker_dir, tmp_path):
    shell = NmrShell()
    svg, png, js = tmp_path / "a.svg", tmp_path / "a.png", tmp_path / "a.json"
    out = run(
        shell,
        f"load {bruker_dir}",
        "ref",
        f"preset {EXPT4_PRESET}",
        f"plot 3.5 6.2 out={svg}",
        f"plot out={png}",
        f"export {js}",
    )
    assert shell.errors == 0, out
    assert svg.read_text().lstrip().startswith("<?xml")
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert [r["name"] for r in json.loads(js.read_text())["s1"]["integrals"]] == ["linalool", "side product", "EtOAc"]


def test_labels_panels_and_stack(bruker_dir, tmp_path):
    shell = NmrShell()
    panels, stack, js = tmp_path / "panels.svg", tmp_path / "stack.svg", tmp_path / "a.json"
    out = run(
        shell,
        f"load {bruker_dir} crude",
        "ref",
        f"load {bruker_dir} pure",
        "ref",
        "label 5.91 Ha",
        "label 3.81 CH(OH)",
        "label 4.12 EtOAc OCH2",
        "labels",
        "unlabel 4.12",
        f"plot 5.7 6.1 3.6 4.3 out={panels}",
        f"stack crude pure 3.4 6.4 norm=5.8:6.0 zoom=3.6:4.3 x=5 out={stack}",
        f"export {js}",
    )
    assert shell.errors == 0, out
    assert shell.labels == {"crude": {}, "pure": {5.91: "Ha", 3.81: "CH(OH)"}}
    drawn = panels.read_text()
    assert "Ha" in drawn and "CH(OH)" in drawn and "EtOAc" not in drawn
    stacked = stack.read_text()
    assert stacked.count("zoom of 3.6–4.3 ppm, ×5 vertical") == 2
    assert "crude" in stacked and "pure" in stacked
    assert json.loads(js.read_text())["pure"]["labels"] == [
        {"ppm": 3.81, "text": "CH(OH)"},
        {"ppm": 5.91, "text": "Ha"},
    ]


def test_figure_argument_errors(bruker_dir, tmp_path):
    shell = NmrShell()
    out = run(
        shell,
        f"load {bruker_dir}",
        "label 99 nowhere",
        "unlabel 5.91",
        "plot 1 2 3",
        "plot 20 30",
        "stack 1 2 zoom=3.6",
        f"stack norm=20:30 out={tmp_path / 'a.png'}",
    )
    assert out.count("error:") == 6
    assert list(tmp_path.glob("*.png")) == []


def test_bad_preset_keeps_the_table(bruker_dir, tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text('{"integrals": [{"name": "x"}]}')
    shell = NmrShell()
    out = run(
        shell, f"load {bruker_dir}", "integrate 3.7 3.9 name=A", f"preset {bad}", "preset /no/such.json", "preset"
    )
    assert out.count("error:") == 3
    assert [r.name for r in shell.integrals["s1"]] == ["A"]


def test_cm3191_example_reads_its_preset(monkeypatch):
    monkeypatch.syspath_prepend(str(EXPT4_PRESET.parent))
    run_expt4 = pytest.importorskip("run_expt4")
    assert run_expt4.QUANT_SIGNALS["EtOAc"] == ("OCH2, δ 4.12", 4.05, 4.20, 2)
    assert run_expt4.MW["linalool"] == 154.25
    assert callable(run_expt4.plot_crude_vs_purified)


def test_errors_are_reported_not_raised(bruker_dir, tmp_path):
    shell = NmrShell()
    out = run(
        shell,
        "peaks 1 2",
        "load /no/such/place",
        f"load {bruker_dir}",
        "frobnicate",
        "integrate a b",
        "j 0.5 0.6",
        f"plot out={tmp_path / 'a.bmp'}",
    )
    assert out.count("error:") == 6
    assert shell.errors == 6
    assert list(tmp_path.glob("*.bmp")) == []


def test_main_runs_commands_and_sets_exit_code(bruker_dir, capsys):
    assert main([str(bruker_dir), "-c", "ref; j 5.86 5.96"]) == 0
    assert "J = 17.3, 10.7 Hz" in capsys.readouterr().out
    assert main([str(bruker_dir), "-c", "nonsense"]) == 1
