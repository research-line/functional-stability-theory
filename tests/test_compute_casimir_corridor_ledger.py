import importlib.util
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest


SCRIPT_PATH = (
    Path(__file__).parents[1] / "scripts" / "yang-mills" / "compute_casimir_corridor_ledger.py"
)
SPEC = importlib.util.spec_from_file_location("compute_casimir_corridor_ledger_under_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _results_by_id():
    return {result.case_id: result for result in MODULE.build_results()}


def test_predefined_gauge_casimir_fixture_derives_exact_reducing_corridor():
    result = _results_by_id()["predefined_gauge_casimir_positive_control"]

    assert result.physical_rank == 3
    assert result.corridor_rank == 1
    assert result.physical_exterior_rank == 2
    assert result.casimir_spectral_separation == pytest.approx(0.75)
    assert result.hamiltonian_spectral_separation == pytest.approx(1.0)
    assert result.hamiltonian_gauge_commutator == pytest.approx(0.0)
    assert result.casimir_gauge_commutator == pytest.approx(0.0)
    assert result.corridor_gauge_leakage == pytest.approx(0.0)
    assert result.hamiltonian_corridor_commutator == pytest.approx(0.0)
    assert result.definition_status == "predefined_casimir_window"
    assert result.gauge_casimir_structure_used
    assert not result.rg_coercivity_used
    assert not result.numerical_residual_used
    assert result.symmetry_provenance_status == "gauge_casimir_only_no_rg_or_residual_input"
    assert result.gap_fingerprint_status == "finite_reducing_corridor_gap_fingerprint"
    assert result.transfer_decision == "blocked_physical_predefinition_not_independently_verified"
    assert result.claim_status == "diagnostic_only_no_yang_mills_claim"


def test_posthoc_window_has_identical_operators_but_is_rejected_as_circular():
    results = _results_by_id()
    positive = results["predefined_gauge_casimir_positive_control"]
    posthoc = results["negative_posthoc_window_same_operators"]

    assert posthoc.physical_rank == positive.physical_rank
    assert posthoc.corridor_rank == positive.corridor_rank
    assert posthoc.hamiltonian_gauge_commutator == positive.hamiltonian_gauge_commutator
    assert posthoc.casimir_gauge_commutator == positive.casimir_gauge_commutator
    assert posthoc.hamiltonian_corridor_commutator == positive.hamiltonian_corridor_commutator
    assert posthoc.definition_status == "posthoc_or_adaptive_window"
    assert posthoc.transfer_decision == "reject_posthoc_or_adaptive_casimir_window"


def test_u1_strong_coupling_control_has_unit_flux_gap_fingerprint_without_claim_upgrade():
    result = _results_by_id()["u1_strong_coupling_loop_flux_positive_control"]

    assert result.physical_rank == 5
    assert result.corridor_rank == 2
    assert result.physical_exterior_rank == 3
    assert result.hamiltonian_gauge_commutator == pytest.approx(0.0)
    assert result.casimir_gauge_commutator == pytest.approx(0.0)
    assert result.hamiltonian_corridor_commutator == pytest.approx(0.0)
    assert result.casimir_spectral_separation == pytest.approx(1.0)
    assert result.hamiltonian_spectral_separation == pytest.approx(1.0)
    assert result.gap_fingerprint_status == "finite_reducing_corridor_gap_fingerprint"
    assert result.transfer_decision == "control_pass_u1_strong_coupling_no_yang_mills_claim"
    assert result.claim_status == "diagnostic_only_no_yang_mills_claim"


def test_casimir_crossing_gauge_boundary_is_rejected():
    result = _results_by_id()["negative_casimir_crosses_gauge_boundary"]

    assert result.hamiltonian_gauge_commutator == pytest.approx(0.0)
    assert result.casimir_gauge_commutator == pytest.approx(0.2)
    assert result.transfer_decision == "reject_casimir_not_gauge_compatible"


def test_hamiltonian_breaking_declared_gauge_sector_is_rejected_first():
    result = _results_by_id()["negative_hamiltonian_breaks_gauge_sector"]

    assert result.hamiltonian_gauge_commutator == pytest.approx(0.2)
    assert result.transfer_decision == "reject_hamiltonian_breaks_gauge_sector"


def test_gauge_invariance_alone_does_not_certify_casimir_corridor():
    result = _results_by_id()["negative_gauge_invariant_but_corridor_mixed"]

    assert result.hamiltonian_gauge_commutator == pytest.approx(0.0)
    assert result.casimir_gauge_commutator == pytest.approx(0.0)
    assert result.hamiltonian_corridor_commutator == pytest.approx(0.2)
    assert result.transfer_decision == "reject_hamiltonian_mixes_casimir_corridor"


def test_reducing_corridor_without_spectral_separation_is_rejected():
    result = _results_by_id()["negative_reducing_corridor_gap_collapsed"]

    assert result.hamiltonian_gauge_commutator == pytest.approx(0.0)
    assert result.hamiltonian_corridor_commutator == pytest.approx(0.0)
    assert result.hamiltonian_spectral_separation == pytest.approx(0.0)
    assert result.gap_fingerprint_status == "reject_spectral_or_structural_gate"
    assert result.transfer_decision == "reject_corridor_spectral_separation_collapse"


def test_only_independently_verified_physical_corridor_becomes_review_eligible():
    case = replace(
        MODULE.build_matched_controls()[0],
        window_definition_source="unit_test_preregistered_physical_casimir_sector",
        predefinition_scope="physical_yang_mills",
        predefinition_certificate_id="UNIT_TEST_INDEPENDENT_PHYSICAL_CERTIFICATE",
        predefinition_certificate_status="independently_verified",
    )

    result = MODULE.audit_case(case)

    assert result.transfer_decision == "eligible_for_analytic_review_no_claim"
    assert result.claim_status == "diagnostic_only_no_yang_mills_claim"


def test_independent_certificate_cannot_upgrade_a_synthetic_scope():
    case = replace(
        MODULE.build_matched_controls()[0],
        predefinition_certificate_id="UNIT_TEST_INDEPENDENT_CERTIFICATE",
        predefinition_certificate_status="independently_verified",
    )

    result = MODULE.audit_case(case)

    assert result.transfer_decision == "blocked_nonphysical_predefinition_scope"


@pytest.mark.parametrize(
    ("case", "expected"),
    [
        (
            replace(MODULE.build_matched_controls()[0], gauge_projection=np.zeros((5, 5))),
            "reject_empty_physical_gauge_sector",
        ),
        (
            replace(MODULE.build_matched_controls()[0], window_lower=9.0, window_upper=10.0),
            "reject_empty_casimir_corridor",
        ),
        (
            replace(MODULE.build_matched_controls()[0], window_lower=-1.0, window_upper=3.0),
            "reject_no_physical_exterior_for_corridor_audit",
        ),
    ],
)
def test_degenerate_corridor_definitions_fail_closed(case, expected):
    result = MODULE.audit_case(case)

    assert result.transfer_decision == expected


@pytest.mark.parametrize("failure_kind", ["non_self_adjoint_h", "non_self_adjoint_c", "non_projection"])
def test_audit_fails_closed_on_invalid_operators(failure_kind):
    case = MODULE.build_matched_controls()[0]
    if failure_kind == "non_self_adjoint_h":
        hamiltonian = case.hamiltonian.astype(float).copy()
        hamiltonian[0, 1] = 0.2
        case = replace(case, hamiltonian=hamiltonian)
        expected = "hamiltonian is not self-adjoint"
    elif failure_kind == "non_self_adjoint_c":
        casimir = case.casimir.astype(float).copy()
        casimir[0, 1] = 0.2
        case = replace(case, casimir=casimir)
        expected = "casimir is not self-adjoint"
    else:
        gauge_projection = case.gauge_projection.astype(float).copy()
        gauge_projection[0, 0] = 0.5
        case = replace(case, gauge_projection=gauge_projection)
        expected = "gauge_projection is not an orthogonal projector"

    with pytest.raises(ValueError, match=expected):
        MODULE.audit_case(case)


def test_cli_writes_claim_neutral_reproducible_ledger(tmp_path, capsys):
    exit_code = MODULE.main(["--output-dir", str(tmp_path)])
    output = capsys.readouterr().out

    assert exit_code == 0
    assert "negative_posthoc_window_same_operators" in output
    assert "negative_gauge_invariant_but_corridor_mixed" in output
    assert "claim_pass=0" in output

    json_path = tmp_path / f"{MODULE.STEM}.json"
    csv_path = tmp_path / f"{MODULE.STEM}.csv"
    md_path = tmp_path / f"{MODULE.STEM}.md"
    assert json_path.is_file()
    assert csv_path.is_file()
    assert md_path.is_file()

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 1
    assert payload["claim_pass"] == 0
    assert payload["scope_status"] == "finite_synthetic_and_u1_strong_coupling_controls_only"
    assert len(payload["rows"]) == 7
    assert {row["transfer_decision"] for row in payload["rows"]} == {
        "blocked_physical_predefinition_not_independently_verified",
        "control_pass_u1_strong_coupling_no_yang_mills_claim",
        "reject_posthoc_or_adaptive_casimir_window",
        "reject_casimir_not_gauge_compatible",
        "reject_hamiltonian_breaks_gauge_sector",
        "reject_hamiltonian_mixes_casimir_corridor",
        "reject_corridor_spectral_separation_collapse",
    }
    assert all(row["claim_status"] == "diagnostic_only_no_yang_mills_claim" for row in payload["rows"])
    markdown = md_path.read_text(encoding="utf-8")
    assert "Gauge preservation" in markdown
    assert "compact U(1) loop-flux control" in markdown
    assert "RG coercivity and residual smallness are not inputs" in markdown
