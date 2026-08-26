import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest


SCRIPT_PATH = (
    Path(__file__).parents[1]
    / "scripts"
    / "yang-mills"
    / "compute_continuum_gate_audit.py"
)
SPEC = importlib.util.spec_from_file_location(
    "compute_continuum_gate_audit_under_test",
    SCRIPT_PATH,
)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_t0_t3_status_is_complete_and_fail_closed():
    rows = MODULE.build_t0_t3_status()

    assert [row.gate_id for row in rows] == ["T0", "T1", "T2", "T3"]
    assert all(row.gate_pass is False for row in rows)
    assert rows[0].current_status == "open_no_continuum_reconstruction_certificate"
    assert rows[1].current_status == "partial_single_link_restricted_control_only"
    assert rows[2].current_status == "open_gt2_defect_sequence_not_derived"
    assert rows[3].current_status == "partial_finite_lattice_rp_only"


def test_fms_scaling_records_existence_rate_gap_without_overflow():
    row = MODULE.audit_fms_transport(64.0)

    assert row.full_log_lipschitz_scale == 64.0
    assert row.typical_log_lipschitz_scale == 8.0
    assert row.full_log_log_bound_scale == 4096.0
    assert row.typical_log_log_bound_scale == 64.0
    assert row.target_contraction == pytest.approx(0.9375)
    assert row.contraction_certificate_pass is False
    assert row.decision == "open_fms_bound_does_not_imply_required_contraction"


@pytest.mark.parametrize("beta", [0.0, -1.0, float("inf"), float("nan")])
def test_fms_rejects_invalid_beta(beta):
    with pytest.raises(ValueError, match="beta"):
        MODULE.audit_fms_transport(beta)


def test_fms_rejects_target_that_is_not_a_contraction():
    with pytest.raises(ValueError, match="target contraction"):
        MODULE.audit_fms_transport(0.01, target_coefficient=0.5)


def test_existing_kingman_false_positive_is_reused_and_rejected():
    rows = {row.case_id: row for row in MODULE.build_kingman_gt2_audits()}
    row = rows["existing_negative_mean_harmonic_corridor"]

    assert row.mean_log_tau < 0.0
    assert row.defect_budget_summable is False
    assert row.kingman_hypotheses_status == "not_verified_for_physical_rg_cocycle"
    assert row.decision == "reject_negative_mean_without_summable_gt2_budget"
    assert row.claim_pass is False


def test_p2_control_proves_only_the_conditional_product_lemma():
    rows = {row.case_id: row for row in MODULE.build_kingman_gt2_audits()}
    row = rows["analytic_p2_gt2_control"]

    assert row.defect_budget_summable is True
    assert row.defect_budget_upper_bound > 0.0
    assert 0.0 < row.conditional_product_lower_bound < 1.0
    assert row.physical_gt2_certificate == "missing_control_only"
    assert row.decision == "conditional_product_lemma_pass_no_physical_rg_input"
    assert row.claim_pass is False


def test_harmonic_gt2_control_breaks_the_limiting_gap():
    rows = {row.case_id: row for row in MODULE.build_kingman_gt2_audits()}
    row = rows["analytic_harmonic_gt2_negative"]

    assert row.mean_log_tau < 0.0
    assert row.defect_budget_summable is False
    assert row.conditional_product_lower_bound == 0.0
    assert row.decision == "reject_harmonic_defect_even_with_negative_mean"


def test_gribov_audit_separates_unfixed_and_gauge_fixed_routes():
    rows = {row.case_id: row for row in MODULE.build_gribov_dz_audits()}
    unfixed = rows["unfixed_compact_link_product"]
    fixed = rows["global_gauge_fixed_slice"]
    null_shortcut = rows["measure_zero_horizon_shortcut"]

    assert unfixed.global_gauge_section_status == "not_used"
    assert unfixed.gribov_effect.startswith("no_direct_obstruction")
    assert unfixed.claim_pass is False
    assert fixed.global_gauge_section_status == "obstructed_by_singer_gribov_ambiguity"
    assert fixed.dobrushin_zegarlinski_status.startswith("blocked")
    assert null_shortcut.seam_or_horizon_capacity_status == (
        "measure_zero_only_not_a_capacity_bound"
    )


def test_payload_has_no_claim_pass_and_records_blocked_proof_note_writeback():
    payload = MODULE.build_payload()

    assert payload["schema_version"] == 1
    assert payload["claim_pass"] == 0
    assert payload["proof_note_writeback_status"] == (
        "blocked_onedrive_cloud_lock_and_host_artifacts"
    )
    assert all(row["gate_pass"] is False for row in payload["t0_t3"])
    assert all(row["claim_pass"] is False for row in payload["fms_transport"])
    assert all(row["claim_pass"] is False for row in payload["kingman_gt2"])
    assert all(row["claim_pass"] is False for row in payload["gribov_dz"])


def test_outputs_are_deterministic_and_machine_readable(tmp_path):
    payload = MODULE.build_payload()
    paths = MODULE.write_outputs(payload, tmp_path)

    assert set(paths) == {"csv", "json", "md"}
    loaded = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert loaded == payload

    with paths["csv"].open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 15
    assert {row["audit_family"] for row in rows} == {
        "t0_t3",
        "fms_transport",
        "kingman_gt2",
        "gribov_dz",
    }
    assert all(row["claim_pass"] == "False" for row in rows)

    markdown = paths["md"].read_text(encoding="utf-8")
    assert "T0--T3 status" in markdown
    assert "FMS quantitative gap" in markdown
    assert "Kingman versus GT2" in markdown
    assert "Gribov copies and Dobrushin--Zegarlinski" in markdown
    assert "claim_pass` is fixed to zero" in markdown


def test_cli_writes_the_audit(tmp_path):
    assert MODULE.main(["--output-dir", str(tmp_path)]) == 0
    assert (tmp_path / f"{MODULE.STEM}.json").is_file()
    assert (tmp_path / f"{MODULE.STEM}.md").is_file()
    assert (tmp_path / f"{MODULE.STEM}.csv").is_file()
