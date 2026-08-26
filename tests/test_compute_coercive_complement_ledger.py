import importlib.util
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest


SCRIPT_PATH = (
    Path(__file__).parents[1] / "scripts" / "yang-mills" / "compute_coercive_complement_ledger.py"
)
SPEC = importlib.util.spec_from_file_location("compute_coercive_complement_ledger_under_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _results_by_id():
    return {result.case_id: result for result in MODULE.build_results()}


def test_positive_control_computes_exact_outer_gap_and_leakage_bound():
    result = _results_by_id()["finite_predefined_positive_control"]

    assert result.s_lambda == pytest.approx(0.0, abs=1e-14)
    assert result.p_lambda == pytest.approx(0.15)
    assert result.g_star == pytest.approx(3.0)
    assert result.residual_over_gap == pytest.approx(0.05)
    assert result.actual_complement_leakage == pytest.approx(0.05)
    assert result.residual_split_bound_holds
    assert result.leakage_bound_holds
    assert result.diagnostic_status == "finite_control_pass"
    assert result.claim_status == "diagnostic_only_no_yang_mills_claim"


def test_bad_scale_control_rejects_residual_ratio_despite_negative_mean_contraction():
    result = _results_by_id()["negative_bad_scale_residual"]

    assert result.mean_log_tau_B < 0.0
    assert result.g_star == pytest.approx(0.25)
    assert result.residual_over_gap == pytest.approx(0.60)
    assert result.leakage_bound_holds
    assert result.diagnostic_status == "reject_residual_over_gap"


def test_gribov_control_rejects_collapsed_outer_gap():
    result = _results_by_id()["negative_gribov_gap_collapse"]

    assert result.mean_log_tau_B < 0.0
    assert result.g_star == pytest.approx(1e-8, rel=1e-6)
    assert result.residual_over_gap == pytest.approx(0.10, rel=1e-6)
    assert result.leakage_bound_holds
    assert result.diagnostic_status == "reject_outer_gap_collapse"


def test_exactly_closed_complement_gap_reports_gap_collapse_not_a_bound():
    case = MODULE.build_matched_controls()[2]
    case = replace(case, hamiltonian=np.diag([0.0, 1.0, 1.0]))

    result = MODULE.audit_case(case)

    assert result.g_star == pytest.approx(0.0)
    assert np.isinf(result.residual_over_gap)
    assert not result.leakage_bound_holds
    assert result.diagnostic_status == "reject_outer_gap_collapse"


@pytest.mark.parametrize("failure_kind", ["non_self_adjoint", "non_projection", "non_reducing"])
def test_audit_fails_closed_before_reporting_compressed_spectrum(failure_kind):
    case = MODULE.build_matched_controls()[0]

    if failure_kind == "non_self_adjoint":
        hamiltonian = case.hamiltonian.astype(float).copy()
        hamiltonian[0, 2] = 0.2
        case = replace(case, hamiltonian=hamiltonian)
        expected = "nicht selbstadjungiert"
    elif failure_kind == "non_projection":
        projection = case.cluster_projection.astype(float).copy()
        projection[1, 1] = 0.5
        case = replace(case, cluster_projection=projection)
        expected = "kein Orthogonalprojektor"
    else:
        hamiltonian = case.hamiltonian.astype(float).copy()
        hamiltonian[1, 2] = hamiltonian[2, 1] = 0.2
        case = replace(case, hamiltonian=hamiltonian)
        expected = "reduziert hamiltonian nicht"

    with pytest.raises(ValueError, match=expected):
        MODULE.audit_case(case)


def test_cli_writes_claim_neutral_reproducible_ledger(tmp_path, capsys):
    exit_code = MODULE.main(["--output-dir", str(tmp_path)])
    output = capsys.readouterr().out

    assert exit_code == 0
    assert "negative_bad_scale_residual" in output
    assert "negative_gribov_gap_collapse" in output
    assert "claim_pass=0" in output

    json_path = tmp_path / f"{MODULE.STEM}.json"
    csv_path = tmp_path / f"{MODULE.STEM}.csv"
    md_path = tmp_path / f"{MODULE.STEM}.md"
    assert json_path.is_file()
    assert csv_path.is_file()
    assert md_path.is_file()

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["claim_pass"] == 0
    assert payload["scope_status"] == "finite_synthetic_controls_only"
    assert len(payload["rows"]) == 3
    assert {row["diagnostic_status"] for row in payload["rows"]} == {
        "finite_control_pass",
        "reject_outer_gap_collapse",
        "reject_residual_over_gap",
    }
    assert all(row["claim_status"] == "diagnostic_only_no_yang_mills_claim" for row in payload["rows"])
    assert "not a physical Yang--Mills predefinition certificate" in md_path.read_text(encoding="utf-8")
