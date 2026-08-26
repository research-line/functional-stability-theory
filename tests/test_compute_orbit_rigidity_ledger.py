import importlib.util
import json
import sys
from dataclasses import replace
from pathlib import Path

import pytest


SCRIPT_PATH = (
    Path(__file__).parents[1] / "scripts" / "yang-mills" / "compute_orbit_rigidity_ledger.py"
)
SPEC = importlib.util.spec_from_file_location("compute_orbit_rigidity_ledger_under_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _results_by_id():
    return {result.case_id: result for result in MODULE.build_results()}


def test_finite_bounded_family_exposes_all_four_selector_fields_without_claim_upgrade():
    result = _results_by_id()["finite_bounded_orbit_family_control"]

    assert result.n_scales == 5
    assert result.bounded_scale_component
    assert result.gauge_orbit_certificate == "FINITE_FIXTURE_GAUSS_ORBIT_A"
    assert result.continuum_residue_gap == pytest.approx(0.27)
    assert result.bad_component_escape == pytest.approx(0.04)
    assert result.good_scale_density == pytest.approx(0.8)
    assert result.finite_gap_floor == pytest.approx(0.44)
    assert result.form_test_status == "finite_form_test_pass"
    assert result.transfer_decision == "form_test_pass_without_physical_origin_or_limit_closure"
    assert result.claim_status == "diagnostic_only_no_yang_mills_claim"


def test_perfect_good_scale_density_cannot_pay_for_unbounded_complexity():
    result = _results_by_id()["negative_high_good_density_unbounded_complexity"]

    assert result.good_scale_density == pytest.approx(1.0)
    assert result.observed_max_component_complexity == 13
    assert result.complexity_overrun == 10
    assert not result.bounded_scale_component
    assert result.transfer_decision == "reject_unbounded_or_incoherent_scale_component"


def test_correct_form_test_without_gauge_origin_remains_diagnostic():
    result = _results_by_id()["negative_stable_gap_missing_gauge_origin"]

    assert result.form_test_status == "finite_form_test_pass"
    assert result.gauge_orbit_certificate == ""
    assert result.transfer_decision == "diagnostic_only_missing_independent_gauge_orbit_certificate"


def test_gauge_certificate_derived_from_gap_signal_is_rejected_as_circular():
    result = _results_by_id()["negative_gauge_certificate_derived_from_gap"]

    assert not result.gauge_orbit_certificate_independent_of_gap_signal
    assert result.transfer_decision == "reject_gauge_orbit_certificate_derived_from_gap_signal"


def test_residue_gap_collapse_blocks_stable_finite_gap_signal():
    result = _results_by_id()["negative_continuum_residue_gap_collapse"]

    assert result.finite_gap_floor >= result.finite_gap_floor_required
    assert result.continuum_residue_gap == pytest.approx(0.0)
    assert result.transfer_decision == "reject_continuum_residue_gap_collapse"


def test_new_remainder_modes_and_bad_escape_fail_closed():
    result = _results_by_id()["negative_bad_component_escape_new_modes"]

    assert result.bad_component_escape == pytest.approx(0.46)
    assert result.new_remainder_modes == 2
    assert result.transfer_decision == "reject_bad_component_escape_or_new_remainder_modes"


def test_posthoc_family_has_same_numbers_but_is_rejected_on_provenance():
    results = _results_by_id()
    positive = results["finite_bounded_orbit_family_control"]
    posthoc = results["negative_posthoc_component_family_same_numbers"]

    assert posthoc.good_scale_density == positive.good_scale_density
    assert posthoc.finite_gap_floor == positive.finite_gap_floor
    assert posthoc.continuum_residue_gap == positive.continuum_residue_gap
    assert posthoc.bad_component_escape == positive.bad_component_escape
    assert posthoc.transfer_decision == "reject_posthoc_or_adaptive_component_family"


def test_physical_label_without_mosco_compact_closure_does_not_upgrade():
    result = _results_by_id()["negative_physical_label_missing_limit_closure"]

    assert result.predefinition_scope == "physical_yang_mills"
    assert result.form_test_status == "finite_form_test_pass"
    assert result.transfer_decision == "diagnostic_only_missing_limit_closure_certificate"


def test_only_independently_verified_physical_origin_and_closure_become_review_eligible():
    case = replace(
        MODULE.build_matched_controls()[0],
        component_family_source="unit_test_physical_predefinition",
        gauge_orbit_certificate="UNIT_TEST_INDEPENDENT_GAUGE_ORBIT",
        gauge_orbit_certificate_status="independently_verified",
        dirichlet_form_status="independently_verified",
        mosco_liminf_status="independently_verified",
        mosco_recovery_status="independently_verified",
        compactness_status="independently_verified",
        limit_closure_certificate_id="UNIT_TEST_INDEPENDENT_LIMIT_CLOSURE",
        limit_closure_certificate_status="independently_verified",
        predefinition_scope="physical_yang_mills",
    )

    result = MODULE.audit_case(case)

    assert result.transfer_decision == "eligible_for_analytic_review_no_claim"
    assert result.claim_status == "diagnostic_only_no_yang_mills_claim"


@pytest.mark.parametrize(
    ("case", "message"),
    [
        (
            replace(MODULE.build_matched_controls()[0], finite_gap_signals=(0.5,)),
            "same nonzero length",
        ),
        (
            replace(
                MODULE.build_matched_controls()[0],
                continuum_residue_gaps=(0.3, 0.2, float("nan"), 0.2, 0.1),
            ),
            "finite nonnegative",
        ),
        (
            replace(MODULE.build_matched_controls()[0], component_complexities=(1, 2, 0, 2, 1)),
            "positive integers",
        ),
        (
            replace(MODULE.build_matched_controls()[0], good_scale_density_floor=1.1),
            "lie in",
        ),
        (
            replace(
                MODULE.build_matched_controls()[0],
                bad_component_escape_values=(0.0, 0.0, -0.1, 0.0, 0.0),
            ),
            "finite nonnegative",
        ),
    ],
)
def test_invalid_scale_ledgers_fail_closed(case, message):
    with pytest.raises(ValueError, match=message):
        MODULE.audit_case(case)


def test_cli_writes_claim_neutral_reproducible_ledger(tmp_path, capsys):
    exit_code = MODULE.main(["--output-dir", str(tmp_path)])
    output = capsys.readouterr().out

    assert exit_code == 0
    assert "negative_high_good_density_unbounded_complexity" in output
    assert "negative_physical_label_missing_limit_closure" in output
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
    assert payload["scope_status"] == "finite_synthetic_orbit_rigidity_controls_only"
    assert len(payload["rows"]) == 8
    assert {row["transfer_decision"] for row in payload["rows"]} == {
        "form_test_pass_without_physical_origin_or_limit_closure",
        "reject_unbounded_or_incoherent_scale_component",
        "diagnostic_only_missing_independent_gauge_orbit_certificate",
        "reject_gauge_orbit_certificate_derived_from_gap_signal",
        "reject_continuum_residue_gap_collapse",
        "reject_bad_component_escape_or_new_remainder_modes",
        "reject_posthoc_or_adaptive_component_family",
        "diagnostic_only_missing_limit_closure_certificate",
    }
    assert all(row["claim_status"] == "diagnostic_only_no_yang_mills_claim" for row in payload["rows"])
    markdown = md_path.read_text(encoding="utf-8")
    for field in (
        "bounded_scale_component",
        "gauge_orbit_certificate",
        "continuum_residue_gap",
        "bad_component_escape",
    ):
        assert field in markdown
    assert "correct finite form test" in markdown
