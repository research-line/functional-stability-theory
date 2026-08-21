"""
RP/OS mirror-sluice transfer ledger for the Yang-Mills project.

This is a research guardrail, not a Yang-Mills proof. It turns the
2026-06-23 Creative-Innovation idea into a small, reproducible status table:
local lattice gap input, reflection-positive blocking, OS-dangerous leakage,
repair/coherence certificate, and continuum-limit status are evaluated as
separate gates.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULT_DIR = ROOT / "_results"
RESULT_DIR.mkdir(exist_ok=True)

FIELDS = [
    "case_id",
    "source_kind",
    "mirror_leak_class",
    "summable_loss_weight",
    "os_crossing_capacity",
    "repair_certificate_id",
    "rp_os_cone_status",
    "continuum_limit_status",
    "area_law_mass_gap_separated",
    "claim_status",
    "notes",
]


CASES = [
    {
        "case_id": "strong_coupling_local_gap_control",
        "source_kind": "finite_lattice_strong_coupling",
        "mirror_leak_class": "internal",
        "summable_loss_weight": 0.031,
        "os_crossing_capacity": 0.0,
        "repair_certificate_id": "DZ_PI_LOCAL_RP_CONTROL",
        "rp_os_cone_status": "pass_finite_lattice",
        "continuum_limit_status": "not_attempted",
        "area_law_mass_gap_separated": "yes",
        "notes": "Positive local control only: supports the proven strong-coupling/zylinder-sector gap, not Clay-level continuum transfer.",
    },
    {
        "case_id": "existing_os_capacity_v2_proxy",
        "source_kind": "toy_or_proxy_ledger",
        "mirror_leak_class": "paired_boundary",
        "summable_loss_weight": 0.124,
        "os_crossing_capacity": 0.019,
        "repair_certificate_id": "COFACTOR_FIELD_PRESENT_BUT_SYNTHETIC",
        "rp_os_cone_status": "schema_pass",
        "continuum_limit_status": "blocked_no_real_rg_data",
        "area_law_mass_gap_separated": "yes",
        "notes": "Good regression guardrail; still proxy data, so no proof or paper-claim upgrade.",
    },
    {
        "case_id": "faizal_shabir_2606_rp_claim",
        "source_kind": "external_preprint_claim",
        "mirror_leak_class": "os_crossing",
        "summable_loss_weight": 0.08,
        "os_crossing_capacity": 0.21,
        "repair_certificate_id": "CLAIMED_TELESCOPING_RP_SLICES_UNAUDITED",
        "rp_os_cone_status": "claimed_not_audited",
        "continuum_limit_status": "claimed_not_audited",
        "area_law_mass_gap_separated": "unclear",
        "notes": "New arXiv 2606.19362 is useful as an adversarial checklist, but needs line-by-line constants, RP-slice stability, and universality audit before use.",
    },
    {
        "case_id": "wilson_loop_confinement_guardrail_2605",
        "source_kind": "survey_guardrail",
        "mirror_leak_class": "internal",
        "summable_loss_weight": 0.0,
        "os_crossing_capacity": 0.0,
        "repair_certificate_id": "AREA_LAW_NOT_MASS_GAP_GATE",
        "rp_os_cone_status": "not_a_transfer_claim",
        "continuum_limit_status": "not_a_transfer_claim",
        "area_law_mass_gap_separated": "yes",
        "notes": "Useful negative boundary: Wilson-loop area law and Wightman/OS mass gap must remain separate ledger columns.",
    },
    {
        "case_id": "su3_positivity_bootstrap_2502",
        "source_kind": "lattice_bootstrap_bounds",
        "mirror_leak_class": "paired_boundary",
        "summable_loss_weight": 0.0,
        "os_crossing_capacity": 0.0,
        "repair_certificate_id": "PSD_RP_WILSON_LOOP_BOUNDS_ONLY",
        "rp_os_cone_status": "bounds_pass",
        "continuum_limit_status": "blocked_no_gap_transfer",
        "area_law_mass_gap_separated": "yes",
        "notes": "Promising source of RP/PSD witness constraints; does not by itself give continuum mass-gap transport.",
    },
    {
        "case_id": "github_coq_mass_gap_claim",
        "source_kind": "github_unreviewed_claim",
        "mirror_leak_class": "os_crossing",
        "summable_loss_weight": 0.5,
        "os_crossing_capacity": 1.0,
        "repair_certificate_id": "UNTRUSTED_COMPLETE_PROOF_CLAIM",
        "rp_os_cone_status": "untrusted",
        "continuum_limit_status": "untrusted",
        "area_law_mass_gap_separated": "unclear",
        "notes": "Repository claim is community signal only; no local proof input without independent mathematical audit.",
    },
]


def classify(row: dict) -> str:
    if row["source_kind"] in {"external_preprint_claim", "github_unreviewed_claim"}:
        return "adversarial_audit_required"
    if row["area_law_mass_gap_separated"] != "yes":
        return "blocked_conflates_area_law_and_gap"
    if row["mirror_leak_class"] == "os_crossing":
        return "blocked_os_crossing_leak"
    if row["continuum_limit_status"] in {"not_attempted", "blocked_no_real_rg_data", "blocked_no_gap_transfer"}:
        return "control_or_schema_only_no_claim"
    if row["rp_os_cone_status"].startswith("pass") and float(row["os_crossing_capacity"]) == 0.0:
        return "local_control_pass_no_clay_claim"
    return "not_validated"


def build_rows() -> list[dict]:
    rows = []
    for case in CASES:
        row = dict(case)
        row["claim_status"] = classify(row)
        rows.append(row)
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_md(rows: list[dict], path: Path) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["claim_status"]] = counts.get(row["claim_status"], 0) + 1

    lines = [
        "# RP/OS Mirror-Sluice Transfer Ledger 2026-06-23",
        "",
        "Status: Creative-Innovation-Guardrail, kein Yang-Mills-Beweisclaim.",
        "",
        "Die Tabelle trennt lokale Gap-Quellen, Spiegel-/OS-Leckklassen, Reparaturzertifikate,",
        "Kontinuumsstatus und die Area-law-vs.-Mass-gap-Wasserlinie.",
        "",
        "## Summary",
        "",
    ]
    for key in sorted(counts):
        lines.append(f"- `{key}`: {counts[key]}")
    lines.extend(
        [
            "",
            "## Ledger",
            "",
            "| case | source | leak | repair | continuum | status |",
            "|---|---|---|---|---|---|",
        ]
    )
    for row in rows:
        lines.append(
            "| {case_id} | {source_kind} | {mirror_leak_class} | {repair_certificate_id} | "
            "{continuum_limit_status} | {claim_status} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Die überlebende neue Idee ist nicht ein stärkerer externer Proof-Claim, sondern ein",
            "engeres Audit-Schema: Ein Transferkandidat muss `mirror_leak_class`,",
            "`repair_certificate_id`, `summable_loss_weight`, `os_crossing_capacity`,",
            "`rp_os_cone_status` und `continuum_limit_status` gleichzeitig ausweisen.",
            "Die neuen 2026-Quellen liefern gute Checklisten, aber keine lokale Claim-Erhöhung.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    base = RESULT_DIR / "RP_OS_MIRROR_SLUICE_LEDGER_2026-06-23"
    write_csv(rows, base.with_suffix(".csv"))
    base.with_suffix(".json").write_text(
        json.dumps({"rows": rows}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    write_md(rows, base.with_suffix(".md"))
    print(f"wrote {base.with_suffix('.md')}")
    print(f"rows={len(rows)} claim_pass=0")


if __name__ == "__main__":
    main()
