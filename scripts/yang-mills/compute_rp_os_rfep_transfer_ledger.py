"""
RP/OS RFEP transfer ledger for the Yang-Mills project.

This is a reproducible guardrail, not a Yang-Mills proof. It lifts the
2026-06-23 mirror-sluice ledger with the 2026-06-24 RFEP Tangential-Waterline
fields: target predefinition, source-axis independence, projection
predefinition, outer-gap source, spectral exactness, residual split, and matched
controls.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
RESULT_DIR = ROOT / "_results"
RESULT_DIR.mkdir(exist_ok=True)

DATE = "2026-06-28"
BASE = RESULT_DIR / f"RP_OS_RFEP_TRANSFER_LEDGER_{DATE}"

FIELDS = [
    "case_id",
    "role",
    "source_kind",
    "mirror_leak_class",
    "repair_certificate_id",
    "summable_loss_weight",
    "os_crossing_capacity",
    "area_law_mass_gap_separated",
    "rp_os_cone_status",
    "continuum_limit_status",
    "target_predefined",
    "source_axis_independent",
    "projection_predefined",
    "outer_gap_source",
    "spectral_exactness_status",
    "total_residual_over_gap",
    "projected_residual_over_tangent_gap",
    "normal_residual",
    "normal_paid_by_constraint",
    "matched_control_family",
    "decision",
    "claim_status",
    "notes",
]


CASES: list[dict[str, Any]] = [
    {
        "case_id": "schema_positive_sc_rp_repair",
        "role": "positive_control_schema",
        "source_kind": "finite_lattice_strong_coupling",
        "mirror_leak_class": "internal",
        "repair_certificate_id": "SC_RP_REPAIR_SCHEMA_PREDECLARED",
        "summable_loss_weight": 0.031,
        "os_crossing_capacity": 0.0,
        "area_law_mass_gap_separated": "yes",
        "rp_os_cone_status": "finite_lattice_pass",
        "continuum_limit_status": "not_attempted",
        "target_predefined": True,
        "source_axis_independent": True,
        "projection_predefined": True,
        "outer_gap_source": "local_strong_coupling_transfer_gap",
        "spectral_exactness_status": "finite_lattice_only",
        "total_residual_over_gap": 0.18,
        "projected_residual_over_tangent_gap": 0.07,
        "normal_residual": 0.0,
        "normal_paid_by_constraint": True,
        "matched_control_family": "u1_or_strong_coupling_positive_control",
        "decision": "schema_pass_no_continuum_claim",
        "notes": "Checks that the new columns can represent a clean finite-lattice control; no Clay-level continuum transfer is asserted.",
    },
    {
        "case_id": "proxy_os_capacity_v2_rfep_lift",
        "role": "legacy_proxy_lift",
        "source_kind": "toy_or_proxy_ledger",
        "mirror_leak_class": "paired_boundary",
        "repair_certificate_id": "COFACTOR_FIELD_PRESENT_BUT_SYNTHETIC",
        "summable_loss_weight": 0.124,
        "os_crossing_capacity": 0.019,
        "area_law_mass_gap_separated": "yes",
        "rp_os_cone_status": "schema_pass",
        "continuum_limit_status": "blocked_no_real_rg_data",
        "target_predefined": True,
        "source_axis_independent": False,
        "projection_predefined": True,
        "outer_gap_source": "synthetic_control_input",
        "spectral_exactness_status": "proxy_only",
        "total_residual_over_gap": 0.42,
        "projected_residual_over_tangent_gap": 0.16,
        "normal_residual": 0.04,
        "normal_paid_by_constraint": True,
        "matched_control_family": "legacy_os_capacity_proxy",
        "decision": "diagnostic_only",
        "notes": "Useful regression guardrail, but still proxy data and not source-axis independent.",
    },
    {
        "case_id": "negative_posterior_target_good_local_gap",
        "role": "matched_negative_control",
        "source_kind": "synthetic_bad_protocol",
        "mirror_leak_class": "internal",
        "repair_certificate_id": "POSTERIOR_TARGET_INVALID",
        "summable_loss_weight": 0.02,
        "os_crossing_capacity": 0.0,
        "area_law_mass_gap_separated": "yes",
        "rp_os_cone_status": "would_pass_if_protocol_valid",
        "continuum_limit_status": "blocked_target_not_predefined",
        "target_predefined": False,
        "source_axis_independent": True,
        "projection_predefined": True,
        "outer_gap_source": "local_gap_signal",
        "spectral_exactness_status": "finite_lattice_only",
        "total_residual_over_gap": 0.12,
        "projected_residual_over_tangent_gap": 0.05,
        "normal_residual": 0.0,
        "normal_paid_by_constraint": True,
        "matched_control_family": "posterior_target",
        "decision": "reject",
        "notes": "Good-looking local signal must fail when the target was chosen after seeing the data.",
    },
    {
        "case_id": "negative_os_crossing_normal_leak",
        "role": "matched_negative_control",
        "source_kind": "synthetic_bad_geometry",
        "mirror_leak_class": "os_crossing",
        "repair_certificate_id": "NORMAL_LEAK_UNPAID",
        "summable_loss_weight": 0.04,
        "os_crossing_capacity": 0.31,
        "area_law_mass_gap_separated": "yes",
        "rp_os_cone_status": "failed_os_crossing",
        "continuum_limit_status": "blocked_os_leak",
        "target_predefined": True,
        "source_axis_independent": True,
        "projection_predefined": True,
        "outer_gap_source": "local_gap_signal",
        "spectral_exactness_status": "finite_lattice_only",
        "total_residual_over_gap": 0.78,
        "projected_residual_over_tangent_gap": 0.11,
        "normal_residual": 0.67,
        "normal_paid_by_constraint": False,
        "matched_control_family": "os_crossing_normal_leak",
        "decision": "reject",
        "notes": "Separates tangential residual from an unbudgeted normal component crossing the OS mirror.",
    },
    {
        "case_id": "negative_missing_repair_certificate",
        "role": "matched_negative_control",
        "source_kind": "synthetic_missing_certificate",
        "mirror_leak_class": "paired_boundary",
        "repair_certificate_id": "",
        "summable_loss_weight": 0.03,
        "os_crossing_capacity": 0.0,
        "area_law_mass_gap_separated": "yes",
        "rp_os_cone_status": "unverifiable_without_repair",
        "continuum_limit_status": "blocked_missing_repair_certificate",
        "target_predefined": True,
        "source_axis_independent": True,
        "projection_predefined": True,
        "outer_gap_source": "local_gap_signal",
        "spectral_exactness_status": "finite_lattice_only",
        "total_residual_over_gap": 0.22,
        "projected_residual_over_tangent_gap": 0.08,
        "normal_residual": 0.02,
        "normal_paid_by_constraint": True,
        "matched_control_family": "missing_repair_certificate",
        "decision": "reject",
        "notes": "A candidate with harmless-looking losses still fails without an independent repair or coherence certificate.",
    },
    {
        "case_id": "negative_spectral_pollution_finite_volume",
        "role": "matched_negative_control",
        "source_kind": "synthetic_finite_volume_false_positive",
        "mirror_leak_class": "internal",
        "repair_certificate_id": "FINITE_VOLUME_ONLY",
        "summable_loss_weight": 0.01,
        "os_crossing_capacity": 0.0,
        "area_law_mass_gap_separated": "yes",
        "rp_os_cone_status": "finite_volume_pass_only",
        "continuum_limit_status": "blocked_no_continuum_exactness",
        "target_predefined": True,
        "source_axis_independent": True,
        "projection_predefined": True,
        "outer_gap_source": "finite_volume_gap",
        "spectral_exactness_status": "finite_volume_false_positive",
        "total_residual_over_gap": 0.09,
        "projected_residual_over_tangent_gap": 0.04,
        "normal_residual": 0.0,
        "normal_paid_by_constraint": True,
        "matched_control_family": "spectral_pollution",
        "decision": "reject",
        "notes": "Prevents a finite-volume spectral gap from masquerading as a continuum mass-gap transport result.",
    },
    {
        "case_id": "faizal_shabir_2606_rp_claim_rfep_audit",
        "role": "adversarial_external_audit",
        "source_kind": "external_preprint_claim",
        "mirror_leak_class": "os_crossing",
        "repair_certificate_id": "CLAIMED_TELESCOPING_RP_SLICES_UNAUDITED",
        "summable_loss_weight": 0.08,
        "os_crossing_capacity": 0.21,
        "area_law_mass_gap_separated": "unclear",
        "rp_os_cone_status": "claimed_not_audited",
        "continuum_limit_status": "claimed_not_audited",
        "target_predefined": False,
        "source_axis_independent": False,
        "projection_predefined": False,
        "outer_gap_source": "external_claim",
        "spectral_exactness_status": "not_locally_verified",
        "total_residual_over_gap": 1.0,
        "projected_residual_over_tangent_gap": 0.4,
        "normal_residual": 0.6,
        "normal_paid_by_constraint": False,
        "matched_control_family": "adversarial_complete_claim",
        "decision": "audit_required",
        "notes": "Useful as an audit checklist only; constants, RP-slice stability, OS limit and area-law/mass-gap separation are not accepted locally.",
    },
]


def classify(row: dict[str, Any]) -> str:
    if row["decision"] == "audit_required":
        return "adversarial_audit_required"
    if row["decision"] == "reject":
        return "matched_negative_control_rejected"
    if not row["target_predefined"]:
        return "blocked_target_not_predefined"
    if not row["projection_predefined"]:
        return "blocked_projection_not_predefined"
    if row["mirror_leak_class"] == "os_crossing":
        return "blocked_os_crossing_leak"
    if not row["repair_certificate_id"]:
        return "blocked_missing_repair_certificate"
    if row["area_law_mass_gap_separated"] != "yes":
        return "blocked_conflates_area_law_and_gap"
    if row["spectral_exactness_status"] == "finite_volume_false_positive":
        return "blocked_spectral_pollution"
    if row["continuum_limit_status"] != "continuum_verified":
        return "diagnostic_only_no_claim"
    if not row["source_axis_independent"]:
        return "diagnostic_only_no_claim"
    if not row["normal_paid_by_constraint"]:
        return "blocked_unpaid_normal_residual"
    return "claim_candidate_needs_human_review"


def build_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for case in CASES:
        row = dict(case)
        row["claim_status"] = classify(row)
        rows.append(row)
    return rows


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_md(rows: list[dict[str, Any]], path: Path) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["claim_status"]] = counts.get(row["claim_status"], 0) + 1

    lines = [
        "# RP/OS RFEP Transfer Ledger 2026-06-28",
        "",
        "Status: FORSCHER-Guardrail, kein Yang-Mills-Beweisclaim.",
        "",
        "Dieses Ledger hebt das RP/OS-Mirror-Sluice-Schema auf die RFEP-Tangential-Waterline:",
        "Ziel, Quelle, Projektion, Außenlücke, Spektral-Exaktheit sowie tangentiale und normale",
        "Residualanteile werden getrennt geführt.",
        "",
        "## Summary",
        "",
    ]
    for key in sorted(counts):
        lines.append(f"- `{key}`: {counts[key]}")
    lines.extend(
        [
            "- `claim_pass`: 0",
            "",
            "## Ledger",
            "",
            "| case | role | leak | repair | target predef. | projection predef. | spectral exactness | decision | status |",
            "|---|---|---|---|---:|---:|---|---|---|",
        ]
    )
    for row in rows:
        lines.append(
            "| {case_id} | {role} | {mirror_leak_class} | {repair_certificate_id} | "
            "{target_predefined} | {projection_predefined} | {spectral_exactness_status} | "
            "{decision} | {claim_status} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Der Beweisstand wird nicht erhöht. Der Fortschritt ist ein schärferes Gate:",
            "Ein künftiger RG-/Gauge-Datenlauf muss vorab definierte Ziel- und Projektionsspalten,",
            "unabhängige Source-Achsen, eine externe Außenlücke, ein Reparaturzertifikat und",
            "eine bezahlte Normalresidual-Komponente liefern. Die Negativkontrollen zeigen,",
            "dass lokale Gap-, RP- oder Area-law-Signale ohne diese Felder verworfen bleiben.",
            "",
            "Nächster Schritt: echte RG-/Gauge-Blockdaten oder eine explizite U(1)-/",
            "Strong-Coupling-Positivkontrolle in genau dieses Schema einspeisen.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    write_csv(rows, BASE.with_suffix(".csv"))
    BASE.with_suffix(".json").write_text(
        json.dumps({"rows": rows, "claim_pass": 0}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    write_md(rows, BASE.with_suffix(".md"))
    print(f"wrote {BASE.with_suffix('.md')}")
    print(f"rows={len(rows)} claim_pass=0")


if __name__ == "__main__":
    main()
