import csv
import importlib.util
import json
import math
import sys
from dataclasses import replace
from pathlib import Path

import pytest


SCRIPT_PATH = (
    Path(__file__).parents[1]
    / "scripts"
    / "yang-mills"
    / "compute_k41_transfer_split_ledger.py"
)
SPEC = importlib.util.spec_from_file_location(
    "compute_k41_transfer_split_ledger_under_test",
    SCRIPT_PATH,
)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _results_by_id():
    return {result.case_id: result for result in MODULE.build_results()}


def test_local_strong_coupling_result_stays_a_positive_control_only():
    result = _results_by_id()["strong_coupling_local_gap_positive_control"]

    assert result.unconditional_positive_control is True
    assert result.local_statement_scope == "strong_coupling_lattice_local"
    assert result.budget_source_id == (
        "compute_os_capacity_ledger.py:strong_coupling_positive_control:k_max=10"
    )
    assert result.local_gap_lower_bound == pytest.approx(0.19893900145093446)
    assert result.mean_log_rg_contraction < 0.0
    assert result.local_gap_or_hessian_pass is True
    assert result.global_LSI_or_Poincare_status == "local_lattice_only_not_global"
    assert result.continuum_transfer_status == "not_attempted_local_control"
    assert result.penalty_bookkeeping_only is False
    assert result.positive_control_status == "local_positive_control_pass"
    assert result.transfer_hypothesis == "not_invoked_local_result_only"
    assert result.k41_style_split is True
    assert result.transfer_decision == "positive_control_pass_local_lattice_no_continuum_claim"


def test_open_transfer_hypothesis_does_not_erase_the_local_positive_result():
    result = _results_by_id()["open_os_continuum_transfer_hypothesis"]

    assert result.unconditional_positive_control is True
    assert result.transfer_hypothesis == "open_os_continuum_transfer_hypothesis"
    assert result.continuum_transfer_status == "open_missing_os_continuum_bridge"
    assert result.split_status == "split_respected_open_transfer_hypothesis"
    assert result.transfer_decision == "diagnostic_only_no_physical_yang_mills_bridge_evidence"


def test_posthoc_envelope_is_rejected_with_identical_numeric_values():
    results = _results_by_id()
    positive = results["strong_coupling_local_gap_positive_control"]
    posthoc = results["negative_posthoc_envelope_same_values"]

    assert posthoc.ratio_min == positive.ratio_min
    assert posthoc.ratio_max == positive.ratio_max
    assert posthoc.defect_budget_total == positive.defect_budget_total
    assert posthoc.state_correlation_envelope == "posthoc_envelope"
    assert posthoc.transfer_decision == "reject_posthoc_state_correlation_envelope"


def test_warm_corridor_gets_a_dedicated_fail_closed_decision():
    result = _results_by_id()["negative_warm_corridor_escape"]

    assert result.warm_corridor_count == 1
    assert result.warm_corridor_status == "warm_corridor_detected"
    assert result.state_correlation_envelope == "outside_predefined_envelope"
    assert result.transfer_decision == "reject_warm_corridor_outside_local_control"


def test_nonwarm_envelope_escape_is_separate_from_warm_corridor():
    result = _results_by_id()["negative_nonwarm_envelope_escape"]

    assert result.warm_corridor_count == 0
    assert result.state_correlation_envelope == "outside_predefined_envelope"
    assert result.transfer_decision == "reject_state_correlation_envelope_escape"


def test_cross_scale_defect_budget_cannot_be_paid_by_a_local_gap():
    result = _results_by_id()["negative_cross_scale_defect_budget_overrun"]

    assert result.local_gap_lower_bound > 0.0
    assert result.defect_budget_total > result.defect_budget_ceiling
    assert result.defect_budget_status == "predefined_defect_budget_overrun"
    assert result.transfer_decision == "reject_cross_scale_defect_budget_overrun"


def test_bridge_derived_from_local_gap_is_circular():
    result = _results_by_id()["negative_bridge_derived_from_local_gap"]

    assert result.bridge_source_status == "derived_from_local_gap_signal"
    assert result.transfer_decision == "reject_circular_bridge_derived_from_local_gap"


def test_local_result_cannot_be_relabelled_as_continuum_theorem():
    result = _results_by_id()["negative_local_result_relabelled_continuum"]

    assert result.claim_mode == "claimed_continuum_transfer"
    assert result.k41_style_split is False
    assert result.split_status == "split_violated_local_result_relabelled_as_continuum"
    assert result.transfer_decision == "reject_scope_inflation_local_result_relabelled_as_continuum"


def test_assumed_local_gap_does_not_count_as_positive_control():
    result = _results_by_id()["negative_assumed_local_gap"]

    assert result.unconditional_positive_control is False
    assert result.positive_control_status == "missing_or_assumed_local_result"
    assert result.transfer_decision == "reject_missing_local_positive_control"


def test_negative_mean_rg_contraction_does_not_pay_a_harmonic_worst_scale_tail():
    result = _results_by_id()["negative_harmonic_worst_scale_corridor"]

    assert result.mean_log_rg_contraction < 0.0
    assert result.local_gap_or_hessian_pass is True
    assert result.global_LSI_or_Poincare_status == "missing"
    assert result.worst_scale_corridor_sum > 0.0
    assert result.worst_scale_corridor_tail_model == "harmonic_non_summable"
    assert result.worst_scale_corridor_summable is False
    assert result.worst_scale_corridor_summability_status == "non_summable_harmonic"
    assert result.continuum_transfer_status == "blocked_non_summable_worst_scale_corridor"
    assert result.transfer_decision == "reject_non_summable_worst_scale_corridor"


def test_unresolved_gribov_corridor_blocks_an_otherwise_good_control():
    result = _results_by_id()["negative_unresolved_gribov_corridor"]

    assert result.local_gap_or_hessian_pass is True
    assert result.worst_scale_corridor_summable is True
    assert result.gribov_corridor_status == "unresolved"
    assert result.worst_scale_or_gribov_status == "gribov_corridor_not_closed"
    assert result.continuum_transfer_status == "blocked_unresolved_gribov_corridor"
    assert result.transfer_decision == "reject_unresolved_gribov_corridor"


def test_penalty_bookkeeping_never_becomes_a_global_bridge():
    result = _results_by_id()["negative_penalty_bookkeeping_only"]

    assert result.local_gap_or_hessian_pass is True
    assert result.penalty_bookkeeping_only is True
    assert result.global_LSI_or_Poincare_status == "local_lattice_only_not_global"
    assert result.continuum_transfer_status == "blocked_penalty_bookkeeping_only"
    assert result.transfer_decision == "reject_penalty_bookkeeping_only"


def test_independent_physical_bridge_can_only_reach_analytic_review():
    case = MODULE._base_case(
        case_id="hypothetical_independent_physical_bridge",
        evidence_scope="physical_yang_mills",
        claim_mode="transfer_hypothesis",
        local_proof_status="independently_verified",
        rp_status="independently_verified",
        os_reconstruction_status="independently_verified",
        continuum_limit_status="independently_verified",
        physical_normalization_status="independently_verified",
        bridge_source_status="independent",
        global_LSI_or_Poincare_status="independently_verified",
        worst_scale_corridor_tail_model="independently_summable",
        gribov_corridor_status="independently_verified_clear",
    )

    result = MODULE.audit_case(case)

    assert result.k41_style_split is True
    assert result.split_status == "split_respected_bridge_ready_for_analytic_review"
    assert result.continuum_transfer_status == "independently_verified"
    assert result.transfer_decision == "eligible_for_analytic_review_no_claim"
    assert result.claim_status == "diagnostic_only_no_yang_mills_or_turbulence_claim"


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("state_correlation_ratios", (1.0, math.nan), "finite values"),
        ("state_correlation_ratios", (1.0, 0.0), "strictly positive"),
        ("defect_budget_terms", (0.1, -0.01), "nonnegative"),
        ("worst_scale_corridor_costs", (0.1, -0.01), "nonnegative"),
        ("local_gap_lower_bound", math.inf, "finite and nonnegative"),
        ("envelope_min", 0.0, "strictly positive"),
    ],
)
def test_numeric_inputs_fail_closed(field, value, message):
    case = replace(MODULE._base_case(), **{field: value})

    with pytest.raises(ValueError, match=message):
        MODULE.audit_case(case)


def test_warm_corridor_flags_must_align_with_ratio_samples():
    case = replace(MODULE._base_case(), warm_corridor_flags=(False,))

    with pytest.raises(ValueError, match="must match"):
        MODULE.audit_case(case)


def test_bundled_payload_is_claim_neutral_and_has_required_split_fields(tmp_path):
    results = MODULE.build_results()
    paths = MODULE.write_outputs(results, tmp_path)
    payload = json.loads(paths["json"].read_text(encoding="utf-8"))

    required = {
        "unconditional_positive_control",
        "budget_source_id",
        "local_gap_or_hessian_pass",
        "global_LSI_or_Poincare_status",
        "continuum_transfer_status",
        "penalty_bookkeeping_only",
        "transfer_hypothesis",
        "state_correlation_envelope",
        "defect_budget_status",
        "warm_corridor_status",
        "k41_style_split",
        "transfer_decision",
        "claim_status",
    }
    assert payload["schema_version"] == 2
    assert payload["claim_pass"] == 0
    assert payload["scope_status"] == "k41_method_transfer_only_no_turbulence_or_yang_mills_claim"
    assert payload["k41_method_source"] == "https://doi.org/10.5281/zenodo.20131305"
    assert payload["rg_budget_source"] == (
        "compute_os_capacity_ledger.py:strong_coupling_positive_control:k_max=10"
    )
    assert payload["harmonic_negative_source"] == (
        "compute_os_capacity_ledger.py:kingman_false_positive_harmonic:k_max=10"
    )
    assert payload["strict_guardrail_fields"] == [
        "local_gap_or_hessian_pass",
        "global_LSI_or_Poincare_status",
        "continuum_transfer_status",
        "penalty_bookkeeping_only",
    ]
    assert len(payload["rows"]) == 12
    assert all(required <= row.keys() for row in payload["rows"])
    assert all(
        row["claim_status"] == "diagnostic_only_no_yang_mills_or_turbulence_claim"
        for row in payload["rows"]
    )


def test_cli_writes_deterministic_csv_json_and_markdown(tmp_path):
    assert MODULE.main(["--output-dir", str(tmp_path)]) == 0

    csv_path = tmp_path / f"{MODULE.STEM}.csv"
    json_path = tmp_path / f"{MODULE.STEM}.json"
    md_path = tmp_path / f"{MODULE.STEM}.md"
    assert csv_path.is_file()
    assert json_path.is_file()
    assert md_path.is_file()

    with csv_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        assert tuple(reader.fieldnames) == MODULE.FIELDS
        assert len(list(reader)) == 12

    markdown = md_path.read_text(encoding="utf-8")
    assert "no Yang--Mills or turbulence claim" in markdown
    assert "Unconditional positive control" in markdown
    assert "Transfer hypothesis" in markdown
    assert "non-summable worst-scale corridor" in markdown
