#!/usr/bin/env python3
"""Build the source-bound FST-II homochirality mechanism taxonomy gate."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-07-12"
RESULTS = ROOT / "_results"
SOURCE_DIR = ROOT / "_sources" / "2026-07-12_forschungsstand"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


SOURCES = {
    "soai": SOURCE_DIR / "trapp_et_al_2025_soai_mechanism.pdf",
    "crystal": SOURCE_DIR / "deck_et_al_2025_homochiral_suspensions.pdf",
    "template": SOURCE_DIR / "li_et_al_2026_crystallization_template_autocatalysis.pdf",
}


ROWS = [
    {
        "mechanism_id": "symmetric_open_drive_null",
        "evidence_class": "exact_minimal_model",
        "bias_source": "initial condition only",
        "amplification_channel": "none; d/dt log(L/D)=0",
        "nonequilibrium_condition": "symmetric feed and dilution",
        "diagnostic_observable": "ee(t) and log(L/D)",
        "matched_negative_control": "h=0 open chemostat",
        "control_status": "executed_2026-07-11",
        "source_anchor": "HOMOCHIRALITY_OPEN_DRIVE_GATE_2026-07-11",
        "claim_ceiling": "Openness alone does not amplify ee in this symmetric minimal model.",
    },
    {
        "mechanism_id": "frank_cross_antagonism",
        "evidence_class": "exact_minimal_model_plus_numeric_controls",
        "bias_source": "nonzero perturbation or noise around the racemic branch",
        "amplification_channel": "mutual antagonism; lambda_chi=h*x_star>0",
        "nonequilibrium_condition": "open autocatalytic growth with h>0",
        "diagnostic_observable": "chiral eigenrate and ee(t)",
        "matched_negative_control": "same parameters with h=0",
        "control_status": "executed_2026-07-11",
        "source_anchor": "HOMOCHIRALITY_OPEN_DRIVE_GATE_2026-07-11",
        "claim_ceiling": "Sufficient branch-instability mechanism in the stated two-enantiomer model.",
    },
    {
        "mechanism_id": "soai_hemiacetalate_autocatalysis",
        "evidence_class": "source_bound_experimental_mechanism",
        "bias_source": "small initial product ee or asymmetric seed",
        "amplification_channel": "product autocatalysis through transient diastereomeric zinc hemiacetalates and dual stereocontrol",
        "nonequilibrium_condition": "diisopropylzinc alkylation under substrate- and solvent-specific homogeneous conditions",
        "diagnostic_observable": "HRMS intermediate profiles, chiral HPLC reaction profiles, substrate consumption, ee(t)",
        "matched_negative_control": "same substrate class without the autocatalytic product/intermediate channel",
        "control_status": "specified_not_executed_in_fst_ii",
        "source_anchor": "Mohler et al. 2025, doi:10.1038/s41467-025-62591-3",
        "claim_ceiling": "Mechanistic and parameter-fit benchmark; not a prebiotic universality claim.",
    },
    {
        "mechanism_id": "stochastic_nucleation_deracemization",
        "evidence_class": "source_bound_experiment_plus_stochastic_model",
        "bias_source": "finite-volume stochastic primary nucleation",
        "amplification_channel": "temperature cycling plus solution racemization and growth/dissolution turnover",
        "nonequilibrium_condition": "supersaturation, cycling, and a racemizing solution",
        "diagnostic_observable": "nucleation-time distribution, crystal counts, solid-phase ee, time to homochirality",
        "matched_negative_control": "remove cycling or solution racemization while keeping the initial racemic solution",
        "control_status": "specified_not_executed_in_fst_ii",
        "source_anchor": "Deck et al. 2025, doi:10.1021/jacs.5c02651",
        "claim_ceiling": "Solid-state stochastic pathway; local handedness is not deterministic global selection.",
    },
    {
        "mechanism_id": "crystallization_template_autocatalysis",
        "evidence_class": "source_bound_experimental_mechanism",
        "bias_source": "kinetically controlled helical bias during nucleation or a chiral seed",
        "amplification_channel": "secondary nucleation, fiber elongation, and template-assisted replication",
        "nonequilibrium_condition": "crystallization-driven reaction with persistent crystalline template fibers",
        "diagnostic_observable": "CD signal, P/M fiber counts, nucleation time, template/precursor conversion",
        "matched_negative_control": "conventional solution crystallization under thermodynamic control",
        "control_status": "source_executed_not_reproduced_in_fst_ii",
        "source_anchor": "Wu et al. 2026, doi:10.1038/s41467-026-70105-y",
        "claim_ceiling": "Material-specific template mechanism; not a universal origin-of-life pathway.",
    },
]


REQUIRED_FIELDS = [
    "mechanism_id",
    "evidence_class",
    "bias_source",
    "amplification_channel",
    "nonequilibrium_condition",
    "diagnostic_observable",
    "matched_negative_control",
    "control_status",
    "source_anchor",
    "claim_ceiling",
]


def main() -> None:
    for path in SOURCES.values():
        if not path.is_file() or path.read_bytes()[:5] != b"%PDF-":
            raise RuntimeError(f"Missing or invalid source PDF: {path}")

    for row in ROWS:
        missing = [field for field in REQUIRED_FIELDS if not row.get(field)]
        if missing:
            raise AssertionError(f"{row['mechanism_id']} missing {missing}")

    ids = {row["mechanism_id"] for row in ROWS}
    if len(ids) != len(ROWS):
        raise AssertionError("Mechanism identifiers must be unique")

    source_hashes = {key: sha256(path) for key, path in SOURCES.items()}
    summary = {
        "date": DATE,
        "rows": len(ROWS),
        "complete_rows": sum(all(row[field] for field in REQUIRED_FIELDS) for row in ROWS),
        "executed_or_source_executed_controls": sum(
            row["control_status"].startswith("executed") or row["control_status"].startswith("source_executed")
            for row in ROWS
        ),
        "fst_ii_unreproduced_literature_rows": sum(
            row["evidence_class"].startswith("source_bound")
            and row["control_status"].endswith("in_fst_ii")
            for row in ROWS
        ),
        "mechanism_collapse_rejected": True,
        "claim_upgrade_allowed": False,
        "decision": "taxonomy_gate_pass_claim_upgrade_blocked",
        "source_sha256": source_hashes,
    }

    RESULTS.mkdir(exist_ok=True)
    stem = RESULTS / f"HOMOCHIRALITY_MECHANISM_TAXONOMY_GATE_{DATE}"

    with stem.with_suffix(".csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REQUIRED_FIELDS)
        writer.writeheader()
        writer.writerows(ROWS)

    with stem.with_suffix(".json").open("w", encoding="utf-8") as handle:
        json.dump({"summary": summary, "rows": ROWS}, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    lines = [
        "# FST-II Homochirality Mechanism Taxonomy Gate",
        "",
        f"**Stand:** {DATE}",
        f"**Entscheidung:** `{summary['decision']}`",
        "",
        "The gate rejects the collapse of distinct homochirality mechanisms into a single Frank-style channel. "
        "It does not upgrade P2: three literature-bound arms still lack an FST-II reproduction or fit.",
        "",
        "| Mechanism | Bias source | Amplification channel | Control status | Claim ceiling |",
        "|---|---|---|---|---|",
    ]
    for row in ROWS:
        lines.append(
            f"| `{row['mechanism_id']}` | {row['bias_source']} | {row['amplification_channel']} | "
            f"`{row['control_status']}` | {row['claim_ceiling']} |"
        )
    lines.extend(
        [
            "",
            "## Machine checks",
            "",
            f"- Rows / complete rows: `{summary['rows']}` / `{summary['complete_rows']}`",
            f"- Executed or source-executed controls: `{summary['executed_or_source_executed_controls']}`",
            f"- Literature arms not executed or reproduced in FST-II: `{summary['fst_ii_unreproduced_literature_rows']}`",
            f"- Mechanism collapse rejected: `{str(summary['mechanism_collapse_rejected']).lower()}`",
            f"- Claim upgrade allowed: `{str(summary['claim_upgrade_allowed']).lower()}`",
            "",
            "## Source hashes",
            "",
        ]
    )
    for key, value in source_hashes.items():
        lines.append(f"- `{key}`: `{value}`")
    stem.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    if summary["rows"] != 5 or summary["complete_rows"] != 5:
        raise AssertionError("Expected five complete mechanism rows")
    if summary["fst_ii_unreproduced_literature_rows"] != 3:
        raise AssertionError("Expected three not-yet-reproduced FST-II literature arms")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
