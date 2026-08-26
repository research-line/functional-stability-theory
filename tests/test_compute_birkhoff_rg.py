import importlib.util
import sys
from pathlib import Path

import pytest


SCRIPT_PATH = Path(__file__).parents[1] / "scripts" / "yang-mills" / "compute_birkhoff_rg.py"
SPEC = importlib.util.spec_from_file_location("compute_birkhoff_rg_under_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_good_scale_boundary_and_prefix_diagnostics():
    diagnostics = MODULE.analyze_good_scales([0.8, 0.95, 0.9, 0.91, 0.7], epsilon=0.1)

    assert diagnostics.threshold == pytest.approx(0.9)
    assert diagnostics.good_mask == (True, False, True, False, True)
    assert diagnostics.good_count == 3
    assert diagnostics.good_density == pytest.approx(0.6)
    assert diagnostics.bad_count == 2
    assert diagnostics.bad_defect_sum == pytest.approx(0.06)
    assert diagnostics.bad_defects == pytest.approx((0.0, 0.05, 0.0, 0.01, 0.0))
    assert diagnostics.prefix_good_density == pytest.approx((1.0, 0.5, 2 / 3, 0.5, 0.6))
    assert diagnostics.cumulative_bad_defect == pytest.approx((0.0, 0.05, 0.05, 0.06, 0.06))


def test_saturation_zones_are_maximal_contiguous_bad_runs():
    diagnostics = MODULE.analyze_good_scales([0.95, 0.96, 0.8, 0.91, 0.92, 0.93], epsilon=0.1)

    assert diagnostics.good_count == 1
    assert diagnostics.bad_defect_sum == pytest.approx(0.17)
    assert diagnostics.max_bad_run == 3
    assert len(diagnostics.saturation_zones) == 2

    first, second = diagnostics.saturation_zones
    assert (first.start_scale, first.end_scale, first.length) == (0, 1, 2)
    assert first.defect_sum == pytest.approx(0.11)
    assert first.max_tau == pytest.approx(0.96)
    assert first.min_gap_to_one == pytest.approx(0.04)
    assert (second.start_scale, second.end_scale, second.length) == (3, 5, 3)
    assert second.defect_sum == pytest.approx(0.06)


@pytest.mark.parametrize(
    ("taus", "epsilon"),
    [
        ([], 0.1),
        ([0.5], 0.0),
        ([0.5], 1.0),
        ([-0.1], 0.1),
        ([1.1], 0.1),
        ([float("nan")], 0.1),
    ],
)
def test_good_scale_diagnostics_reject_invalid_inputs(taus, epsilon):
    with pytest.raises(ValueError):
        MODULE.analyze_good_scales(taus, epsilon)


def test_small_no_plot_run_reports_claim_neutral_diagnostics(monkeypatch, capsys):
    monkeypatch.setattr(MODULE, "render_plot", lambda *args: pytest.fail("--no-plot wurde ignoriert"))

    exit_code = MODULE.main(
        ["--no-plot", "--epsilon", "0.1", "--n-levels", "3", "--n-bins", "4", "--betas", "1.0"]
    )
    output = capsys.readouterr().out

    assert exit_code == 0
    assert "rho_epsilon=" in output
    assert "bad defects: Summe=" in output
    assert "Saturationszonen:" in output
    assert "kein Kingman-Claim" in output
    assert "beweisen Kingman-Voraussetzungen" in output


def test_diagnostic_plot_can_be_rendered_outside_repository(monkeypatch, tmp_path):
    pytest.importorskip("matplotlib")
    beta = 1.0
    results = MODULE.rg_cascade(beta, n_levels=3, n_bins=4)
    diagnostics = MODULE.analyze_good_scales([result["tau_B"] for result in results], epsilon=0.1)
    monkeypatch.setattr(MODULE, "Path", lambda _: tmp_path / "compute_birkhoff_rg.py")

    output_path = MODULE.render_plot({beta: results}, {beta: diagnostics}, [beta], epsilon=0.1)

    assert output_path == tmp_path / "compute_birkhoff_rg.png"
    assert output_path.is_file()
    assert output_path.stat().st_size > 0
