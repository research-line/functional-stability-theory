#!/usr/bin/env python3
"""Audit FST-II claim candidates against a dissipation/topology waterline.

This script is deliberately conservative: it does not try to prove P3, P5, or
homochirality claims. It records which independent channels would be required
before a stronger claim could pass.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any


RUN_DATE = "2026-06-19"
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "_results"
OUT_BASE = f"DISSIPATION_TOPOLOGY_WATERLINE_LEDGER_{RUN_DATE}"


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def stable_count(result: dict[str, Any]) -> int:
    return sum(1 for point in result.get("fixed_points", []) if point.get("stable"))


def best_stability_status(result: dict[str, Any]) -> str:
    statuses = Counter(point.get("stability_status", "unknown") for point in result.get("fixed_points", []))
    if not statuses:
        return "unknown"
    return ", ".join(f"{name}:{count}" for name, count in sorted(statuses.items()))


def classify(row: dict[str, Any]) -> str:
    required = [
        "independent_dissipation",
        "independent_selection",
        "topology_varied",
        "stability_margin_measured",
        "negative_control_passed",
    ]
    if all(row[key] for key in required):
        return "strong_pass"
    if row["case_type"] == "current_proxy":
        return "reject_strong_pass"
    if row["case_type"] == "source_guardrail":
        return "literature_guardrail"
    if row["case_type"] == "required_control":
        return "required_control_missing"
    return "not_validated"


def blockers(row: dict[str, Any]) -> list[str]:
    checks = [
        ("independent_dissipation", "dissipation_not_independent"),
        ("independent_selection", "selection_not_independent"),
        ("topology_varied", "topology_not_varied"),
        ("stability_margin_measured", "stability_margin_missing"),
        ("negative_control_passed", "negative_control_missing"),
    ]
    return [label for key, label in checks if not row[key]]


def build_rows() -> list[dict[str, Any]]:
    ii_calc = load_json(ROOT / "scripts" / "output" / "ii_calc_01_results.json")
    rep3 = load_json(ROOT / "results" / "replicator_azoarcus_3species.json")
    rep4 = load_json(ROOT / "results" / "replicator_azoarcus_4species.json")

    tangent_alpha = ii_calc["stability"].get("tangent_max_real_eigenvalue")
    s_dot = ii_calc["entropy_production"].get("S_dot")

    rows: list[dict[str, Any]] = [
        {
            "case_id": "current_ii_calc_01_three_species",
            "case_type": "current_proxy",
            "claim_area": "P3/MEPP",
            "evidence": f"internal ODE flux proxy; S_dot={s_dot:.6g}; alpha_T={tangent_alpha:.6g}",
            "independent_dissipation": False,
            "independent_selection": False,
            "topology_varied": False,
            "stability_margin_measured": True,
            "negative_control_passed": False,
            "next_step": "repeat with independent calorimetry/flux observable and topology ensemble",
        },
        {
            "case_id": "current_replicator_azoarcus_3species",
            "case_type": "current_proxy",
            "claim_area": "P3/MEPP",
            "evidence": (
                f"{stable_count(rep3)} stable fixed point(s); "
                f"statuses={best_stability_status(rep3)}; sigma_proxy guardrail present"
            ),
            "independent_dissipation": False,
            "independent_selection": False,
            "topology_varied": False,
            "stability_margin_measured": True,
            "negative_control_passed": False,
            "next_step": "archive as toy throughput proxy, not thermodynamic entropy evidence",
        },
        {
            "case_id": "current_replicator_azoarcus_4species",
            "case_type": "current_proxy",
            "claim_area": "P3/P4",
            "evidence": (
                f"{stable_count(rep4)} stable fixed point(s); "
                f"statuses={best_stability_status(rep4)}; neutral boundary channel recorded"
            ),
            "independent_dissipation": False,
            "independent_selection": False,
            "topology_varied": False,
            "stability_margin_measured": True,
            "negative_control_passed": False,
            "next_step": "use only as invasion/boundary diagnostic until independent dissipation exists",
        },
        {
            "case_id": "proposed_p3_calorimetry_sequence_design",
            "case_type": "design",
            "claim_area": "P3/MEPP",
            "evidence": "design separates population/sequence channel from calorimetry or heat-flux channel",
            "independent_dissipation": True,
            "independent_selection": True,
            "topology_varied": False,
            "stability_margin_measured": False,
            "negative_control_passed": False,
            "next_step": "define matched networks and calibrate calorimetry before claim testing",
        },
        {
            "case_id": "proposed_raf_topology_ensemble",
            "case_type": "design",
            "claim_area": "topology/P3",
            "evidence": "CatReNet/RAF tool path can provide RAF/CAF/pseudo-RAF signatures and shuffled controls",
            "independent_dissipation": False,
            "independent_selection": False,
            "topology_varied": True,
            "stability_margin_measured": False,
            "negative_control_passed": False,
            "next_step": "build RAF ensemble with degree-profile and stoichiometry-shuffle controls",
        },
        {
            "case_id": "proposed_p5_hopf_floquet_window",
            "case_type": "design",
            "claim_area": "P5/spatial patterns",
            "evidence": "Hopf/Floquet observables can replace purely visual spiral-wave language",
            "independent_dissipation": False,
            "independent_selection": False,
            "topology_varied": True,
            "stability_margin_measured": False,
            "negative_control_passed": False,
            "next_step": "add reaction-diffusion model with Hopf distance and Floquet multipliers",
        },
        {
            "case_id": "required_homochirality_no_cross_inhibition_control",
            "case_type": "required_control",
            "claim_area": "homochirality",
            "evidence": "open driven model without Frank/Viedma-style cross-inhibition must not auto-pass",
            "independent_dissipation": False,
            "independent_selection": False,
            "topology_varied": False,
            "stability_margin_measured": False,
            "negative_control_passed": False,
            "next_step": "run negative control before any stronger homochirality wording",
        },
        {
            "case_id": "source_tamano_saito_entropy_criticality",
            "case_type": "source_guardrail",
            "claim_area": "CRN entropy production",
            "evidence": "arXiv:2605.03192 suggests entropy-production criticality is a separate observable layer",
            "independent_dissipation": False,
            "independent_selection": False,
            "topology_varied": False,
            "stability_margin_measured": False,
            "negative_control_passed": False,
            "next_step": "add entropy-production criticality field; do not infer selection from it alone",
        },
        {
            "case_id": "source_gagrani_path_ranking",
            "case_type": "source_guardrail",
            "claim_area": "thermodynamic pathway rank",
            "evidence": "arXiv:2506.23496 supports path ranking as its own channel, not selection proof",
            "independent_dissipation": False,
            "independent_selection": False,
            "topology_varied": False,
            "stability_margin_measured": False,
            "negative_control_passed": False,
            "next_step": "keep path_rank_channel separate from selection_channel",
        },
    ]

    for row in rows:
        row["decision"] = classify(row)
        row["blockers"] = ";".join(blockers(row))
    return rows


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id",
        "case_type",
        "claim_area",
        "decision",
        "blockers",
        "evidence",
        "independent_dissipation",
        "independent_selection",
        "topology_varied",
        "stability_margin_measured",
        "negative_control_passed",
        "next_step",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_markdown(rows: list[dict[str, Any]], path: Path) -> None:
    counts = Counter(row["decision"] for row in rows)
    lines = [
        f"# Dissipation-Topology-Waterline-Ledger - {RUN_DATE}",
        "",
        "Status: konservatives Audit, kein Claim-Upgrade.",
        "",
        "## Summary",
        "",
        f"- Fälle gesamt: {len(rows)}",
        f"- Strong passes: {counts.get('strong_pass', 0)}",
        f"- Reject strong pass: {counts.get('reject_strong_pass', 0)}",
        f"- Not validated: {counts.get('not_validated', 0)}",
        f"- Required controls missing: {counts.get('required_control_missing', 0)}",
        f"- Literature guardrails: {counts.get('literature_guardrail', 0)}",
        "",
        "Interpretation: Der aktuelle Projektstand bleibt P3-/P5-/Homochiralitäts-heuristisch. "
        "Die neue Wasserlinie ist nützlich, weil sie exakt zeigt, welche Kanäle fehlen.",
        "",
        "## Ledger",
        "",
        "| case | area | decision | blockers | next step |",
        "|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            "| {case_id} | {claim_area} | {decision} | {blockers} | {next_step} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Wichtigster Befund",
            "",
            "Kein aktueller Toy-Fall passiert die starke Wasserlinie, weil Dissipation und "
            "Spiel-/Selektionsstatus nicht unabhängig gemessen werden und die Topologie nicht "
            "als Ensemble variiert wird. Das ist kein Rückschritt, sondern eine präzisere "
            "Arbeitsdefinition der nächsten mathematisch-experimentellen Hürde.",
            "",
            "## Nächste Arbeitsarme",
            "",
            "1. `P3`: Sequenzierung/Populationen plus Kalorimetrie oder Heat-Flux-Kanal als unabhängige Messung.",
            "2. `Topologie`: RAF-/CAF-/pseudo-RAF-Ensemble mit Gradprofil- und Stöchiometriekontrollen.",
            "3. `P5`: Reaktions-Diffusionsmodell mit Hopf-Abstand und Floquet-Multiplikatoren.",
            "4. `Homochiralität`: Negativkontrolle ohne Frank-/Viedma-artige Kreuzinhibition.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    summary = {
        "run_date": RUN_DATE,
        "project": "FST-II_Chemical",
        "claim_upgrade": False,
        "counts": dict(Counter(row["decision"] for row in rows)),
        "rows": rows,
    }

    (RESULTS / f"{OUT_BASE}.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_csv(rows, RESULTS / f"{OUT_BASE}.csv")
    write_markdown(rows, RESULTS / f"{OUT_BASE}.md")

    print(json.dumps(summary["counts"], sort_keys=True))
    print(f"claim_upgrade={summary['claim_upgrade']}")


if __name__ == "__main__":
    main()
