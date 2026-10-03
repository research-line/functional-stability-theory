#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S-Unit / Frey shadow audit for Beal root-survivor diagnostics.

This is a guardrail runner. It classifies the current local-pass survivors
against external generalized-Fermat source classes and matched controls. It does
not prove Beal's conjecture or close HD3'.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


RUN_DATE = "2026-06-22"


@dataclass(frozen=True)
class SourceClass:
    source_id: str
    title: str
    arxiv: str
    method_axis: str
    signature_scope: str
    beal_relevance: str
    local_path: str
    source_url: str


@dataclass
class AuditRow:
    row_id: str
    row_type: str
    rank: str
    signature: str
    tuple_data: str
    best_pass_Q: str
    first_empty_after_best: str
    low_sections_at_best: str
    duplicate_key_count: str
    source_fit_sahoo_2p2qr: str
    source_fit_cazorla_2rp: str
    source_fit_sahoo_zl: str
    source_fit_billerey_13_13_7: str
    source_fit_abc_quality: str
    q97_status: str
    anchor_status: str
    matched_control_status: str
    verdict: str
    claim_level: str
    next_data_needed: str


SOURCES: Sequence[SourceClass] = (
    SourceClass(
        source_id="sahoo_2509_2p2qr",
        title="Effective Generalized Fermat equation of signature (2p, 2q, r) with odd narrow class number",
        arxiv="2509.21083",
        method_axis="S-unit + Hilbert modular / modularity",
        signature_scope="x^(2p)+y^(2q)=z^r, p,q>=5, r prime, class-number hypotheses",
        beal_relevance="direct only for two even exponent channels with large prime halves; otherwise method anchor",
        local_path="_sources/2026-05-19_creative_innovation/2025_Sahoo_Effective_GFE_2p_2q_r_arxiv2509.21083.pdf",
        source_url="https://arxiv.org/abs/2509.21083",
    ),
    SourceClass(
        source_id="cazorla_2605_2rp",
        title="The generalized Fermat equation Ax^2 + By^r = Cz^p and applications",
        arxiv="2605.02632",
        method_axis="Frey hyperelliptic curves / Darmon program",
        signature_scope="A*x^2 + B*y^r = C*z^p and signature transfers involving (2,r,p)",
        beal_relevance="shadow fit when one Beal exponent can be folded to 2; not a primitive HD3' bridge by itself",
        local_path="_sources/2026-05-19_creative_innovation/2026_CazorlaGarcia_Koutsianas_VillagraTorcomian_Ax2_Byr_Czp_arxiv2605.02632.pdf",
        source_url="https://arxiv.org/abs/2605.02632",
    ),
    SourceClass(
        source_id="moriwaki_2601_arithmetic_dynamics",
        title="Arithmetic Dynamics and Generalized Fermat's Conjecture",
        arxiv="2601.08207",
        method_axis="arithmetic dynamics / adelic heights / generalized Fermat over arithmetic function fields",
        signature_scope="height-zero Fermat property for dynamical pullback families, not integer Beal signatures",
        beal_relevance="remote height analogy only; no direct S-unit, Frey, Beal, or HD3' anchor",
        local_path="_sources/2026-05-19_creative_innovation/2026_Moriwaki_Arithmetic_Dynamics_Generalized_Fermat_arxiv2601.08207.pdf",
        source_url="https://arxiv.org/abs/2601.08207",
    ),
    SourceClass(
        source_id="sahoo_2605_zl",
        title="Generalized Fermat equation over cyclotomic Z_l-extensions of totally real fields",
        arxiv="2605.20860",
        method_axis="asymptotic FLT / unit coefficient classes in Z_l-towers",
        signature_scope="Ax^p+By^p+Cz^p=0 and restricted unit / power coefficient classes",
        beal_relevance="unit vocabulary and asymptotic context, not a heterogeneous Beal-survivor anchor",
        local_path="_sources/2026-06-02_forscher_check/2026_Sahoo_GFE_cyclotomic_Zl_extensions_arxiv2605.20860.pdf",
        source_url="https://arxiv.org/abs/2605.20860",
    ),
    SourceClass(
        source_id="billerey_2510_x13_y13_3z7",
        title="Revisiting the Fermat-type equation x^13 + y^13 = 3z^7",
        arxiv="2510.13773",
        method_axis="unit sieve + multi-Frey modular method + level raising",
        signature_scope="specific equation x^13+y^13=3z^7",
        beal_relevance="strong method warning: unit sieve alone had an obstruction; no direct fit to current survivors",
        local_path="_sources/2026-06-22_forscher_check/2025_Billerey_Chen_Dembele_Dieulefait_Freitas_x13_y13_3z7_arxiv2510.13773.pdf",
        source_url="https://arxiv.org/abs/2510.13773",
    ),
    SourceClass(
        source_id="sankaran_2606_abc_quality",
        title="Variants on the abc-Conjecture using Alternative Quality Metrics",
        arxiv="2606.08416",
        method_axis="abc quality metrics / smoothness phase heuristics / Frey-curve context",
        signature_scope="abc triples, not a Beal signature family",
        beal_relevance="context for radical-quality accounting only; no HD3' or S-unit closure",
        local_path="_sources/2026-06-22_forscher_check/2026_Sankaran_abc_alternative_quality_metrics_arxiv2606.08416.pdf",
        source_url="https://arxiv.org/abs/2606.08416",
    ),
    SourceClass(
        source_id="hanners_zenodo_beal_claim_low_reliability",
        title="Resolution of Beal's Conjecture",
        arxiv="Zenodo 10.5281/zenodo.15338623",
        method_axis="entropy / harmonic-coherence claim",
        signature_scope="purported global Beal claim without source-bound root, S-unit, or Frey anchors",
        beal_relevance="low-reliability negative control; not used as evidence",
        local_path="_sources/2026-05-19_creative_innovation/2026_Hanners_Resolution_of_Beals_Conjecture_Zenodo19401153_claim_LOW_RELIABILITY.pdf",
        source_url="https://zenodo.org/records/15338623",
    ),
)


def load_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def signature(row: Dict[str, str]) -> str:
    return f"({row['x']},{row['y']},{row['z']})"


def tuple_data(row: Dict[str, str]) -> str:
    return f"A={row['A']},x={row['x']},B={row['B']},y={row['y']},C={row['C']},z={row['z']}"


def has_even_side(row: Dict[str, str]) -> bool:
    return int(row["x"]) % 2 == 0 or int(row["y"]) % 2 == 0


def fits_sahoo_2p2qr(row: Dict[str, str]) -> str:
    x = int(row["x"])
    y = int(row["y"])
    z = int(row["z"])
    if x % 2 == 0 and y % 2 == 0 and x // 2 >= 5 and y // 2 >= 5 and z >= 5:
        return "direct_signature_scope"
    if x % 2 == 0 or y % 2 == 0:
        return "weak_shadow_even_exponent_but_scope_mismatch"
    return "no_signature_fit"


def fits_cazorla_2rp(row: Dict[str, str]) -> str:
    if has_even_side(row):
        return "shadow_fit_after_even_exponent_folding_only"
    return "no_signature_fit"


def fits_sahoo_zl(row: Dict[str, str]) -> str:
    if row["x"] == row["y"] == row["z"]:
        return "possible_equal_exponent_context_only"
    return "no_heterogeneous_signature_fit"


def fits_billerey(row: Dict[str, str]) -> str:
    if row["x"] == "13" and row["y"] == "13" and row["z"] == "7":
        return "direct_specific_signature"
    return "no_specific_signature_fit"


def q97_status(row: Dict[str, str]) -> str:
    first_empty = int(row.get("first_empty_after_best") or 0)
    if int(row.get("best_pass_Q") or 0) >= 97 and first_empty == 0:
        return "q97_plus_open"
    if first_empty and first_empty <= 97:
        return f"fails_at_q{first_empty}"
    return "not_tested_to_q97"


def classify_row(row: Dict[str, str]) -> AuditRow:
    sahoo_2p2qr = fits_sahoo_2p2qr(row)
    cazorla_2rp = fits_cazorla_2rp(row)
    sahoo_zl = fits_sahoo_zl(row)
    billerey = fits_billerey(row)
    q97 = q97_status(row)

    duplicate_count = int(row.get("duplicate_key_count") or 1)
    direct_anchor = sahoo_2p2qr == "direct_signature_scope" or billerey == "direct_specific_signature"
    weak_shadow = "shadow" in sahoo_2p2qr or "shadow" in cazorla_2rp

    if duplicate_count > 1:
        anchor_status = "duplicate_power_anchor_not_independent"
        matched_control_status = "duplicate_channel_control_rejects_independence"
        verdict = "reject_as_independent_survivor"
        next_data_needed = "Replace duplicate power presentation by one canonical family row before any source matching."
    elif direct_anchor and q97 == "q97_plus_open":
        anchor_status = "source_bound_candidate_needs_full_local_global_audit"
        matched_control_status = "controls_missing"
        verdict = "blocked_unfilled"
        next_data_needed = "Run source-bound S-unit data and matched negative controls, then re-check Q97-plus survival."
    elif weak_shadow:
        anchor_status = "method_shadow_only"
        matched_control_status = "negative_controls_required"
        verdict = "local_transition_no_s_unit_anchor"
        next_data_needed = "Either construct a synthetic Q97-plus survivor with pre-registered root section or obtain direct S-unit data for this signature."
    else:
        anchor_status = "no_external_signature_anchor"
        matched_control_status = "negative_controls_required"
        verdict = "local_transition_no_s_unit_anchor"
        next_data_needed = "Search in known signature families or build synthetic exact-section controls; do not upgrade from local passes."

    if q97.startswith("fails_at_q"):
        verdict = "local_transition_fails_before_q97"

    return AuditRow(
        row_id=f"survivor_rank_{row['rank']}",
        row_type="natural_q61_survivor",
        rank=row["rank"],
        signature=signature(row),
        tuple_data=tuple_data(row),
        best_pass_Q=row["best_pass_Q"],
        first_empty_after_best=row["first_empty_after_best"],
        low_sections_at_best=row["low_sections_at_best"],
        duplicate_key_count=row["duplicate_key_count"],
        source_fit_sahoo_2p2qr=sahoo_2p2qr,
        source_fit_cazorla_2rp=cazorla_2rp,
        source_fit_sahoo_zl=sahoo_zl,
        source_fit_billerey_13_13_7=billerey,
        source_fit_abc_quality="context_only_not_signature_anchor",
        q97_status=q97,
        anchor_status=anchor_status,
        matched_control_status=matched_control_status,
        verdict=verdict,
        claim_level="no_claim_upgrade",
        next_data_needed=next_data_needed,
    )


def source_matrix_lines() -> List[str]:
    lines = [
        "| source_id | arXiv | method_axis | signature_scope | Beal-Relevanz |",
        "|---|---|---|---|---|",
    ]
    for source in SOURCES:
        lines.append(
            f"| `{source.source_id}` | `{source.arxiv}` | {source.method_axis} | "
            f"{source.signature_scope} | {source.beal_relevance} |"
        )
    return lines


def summarize(rows: Sequence[AuditRow]) -> Dict[str, int]:
    unique_keys = {(row.signature, row.tuple_data) for row in rows}
    return {
        "survivor_rows": len(rows),
        "unique_tuple_rows": len(unique_keys),
        "q97_plus_survivors": sum(1 for row in rows if row.q97_status == "q97_plus_open"),
        "rows_with_direct_s_unit_anchor": sum(
            1
            for row in rows
            if row.source_fit_sahoo_2p2qr == "direct_signature_scope"
            or row.source_fit_billerey_13_13_7 == "direct_specific_signature"
        ),
        "rows_with_method_shadow_only": sum(1 for row in rows if row.anchor_status == "method_shadow_only"),
        "duplicate_reject_rows": sum(1 for row in rows if row.anchor_status == "duplicate_power_anchor_not_independent"),
        "claim_upgrades": 0,
    }


def write_csv(rows: Sequence[AuditRow], path: Path) -> None:
    fields = list(AuditRow.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(rows: Sequence[AuditRow], summary: Dict[str, int], path: Path) -> None:
    payload = {
        "run_date": RUN_DATE,
        "status": "guardrail_audit_no_claim_upgrade",
        "summary": summary,
        "sources": [asdict(source) for source in SOURCES],
        "rows": [asdict(row) for row in rows],
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_markdown(rows: Sequence[AuditRow], summary: Dict[str, int], path: Path, input_path: Path) -> None:
    lines: List[str] = [
        "# S-Unit-Shadow-Audit -- Beal Height Dominance",
        "",
        f"Datum: {RUN_DATE}",
        "",
        "## Scope",
        "",
        "Dieser FORSCHER-Lauf prüft, ob die stärksten lokalen `Q=61`-Survivors aus der Beal-Root-Coherence-Diagnostik durch aktuelle S-Unit-, Frey- oder `abc`-Quellen gestützt werden. Das Audit ist ein Guardrail-Runner: kein Beal-Beweis, keine HD3'-Schließung, kein TeX-/PDF-Build und kein Upload.",
        "",
        "## Aktueller Quellenstand",
        "",
        *source_matrix_lines(),
        "",
        "## Aggregate",
        "",
        f"- Input-Survivors: `{input_path}`",
        f"- Survivor-Zeilen: `{summary['survivor_rows']}`",
        f"- Q97-plus-Survivors: `{summary['q97_plus_survivors']}`",
        f"- Direkte S-Unit-/Frey-Signaturanker: `{summary['rows_with_direct_s_unit_anchor']}`",
        f"- Nur methodische Shadow-Fits: `{summary['rows_with_method_shadow_only']}`",
        f"- Duplikat-Reject-Zeilen: `{summary['duplicate_reject_rows']}`",
        f"- Claim-Upgrades: `{summary['claim_upgrades']}`",
        "",
        "## Survivor-Audit",
        "",
        "| row_id | Signatur | best Q | erste leere q danach | Duplikate | Sahoo (2p,2q,r) | Cazorla (2,r,p) | Q97 | Verdict |",
        "|---|---|---:|---:|---:|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| `{row.row_id}` | `{row.signature}` | {row.best_pass_Q} | {row.first_empty_after_best} | "
            f"{row.duplicate_key_count} | `{row.source_fit_sahoo_2p2qr}` | `{row.source_fit_cazorla_2rp}` | "
            f"`{row.q97_status}` | `{row.verdict}` |"
        )
    lines.extend(
        [
            "",
            "## Befund",
            "",
            "Die stärksten natürlichen Survivor bleiben lokale Übergangsfälle: keine Zeile überlebt als `Q97-plus`-Kandidat, und keine Zeile besitzt einen direkten externen S-Unit-/Frey-Signaturanker. Die drei Zeilen mit geradem Exponenten erhalten höchstens einen methodischen Shadow-Fit zur `(2,r,p)`- bzw. `(2p,2q,r)`-Literatur; das ist nützlich für die Sprache des nächsten Audits, aber kein Strukturbeleg.",
            "",
            "Der neue Quellenfund `x^13+y^13=3z^7` ist als Warnanker wichtig: Selbst in einer engen Spezialsignatur war ein reiner Unit-Sieve-Pfad fehleranfällig und musste mit multi-Frey-/Level-Raising-Information kombiniert werden. Für `beal_height_dominance` bedeutet das: S-Unit-Shadow-Felder dürfen nur als source-bound gelten, wenn die tatsächliche Signatur, die lokalen Bedingungen und die Negativkontrollen gemeinsam vorliegen.",
            "",
            "Der `abc`-Metriken-Preprint stützt nur das Radical-/Qualitätsvokabular. Er liefert keinen Beal- oder HD3'-Hebel und wird deshalb als Kontextanker, nicht als Claim-Evidenz geführt.",
            "",
            "## Nächster Schritt",
            "",
            "Das nächste prüfbare Gate ist ein synthetischer `Q97-plus`-Kontrolllauf mit vorregistrierter Root-Sektion und geschuffelter Branch-Map. Alternativ müsste eine bekannte Signaturfamilie direkt in eine der aktuellen Survivor-Signaturen passen und die S-Unit-/Frey-Daten source-bound liefern.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    project = here.parent
    default_input = here / "ABOVE_THRESHOLD_SURVIVOR_SEARCH_2026-06-02.csv"
    default_prefix = project / "_results" / f"S_UNIT_SHADOW_AUDIT_{RUN_DATE}"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=default_input)
    parser.add_argument("--csv", type=Path, default=default_prefix.with_suffix(".csv"))
    parser.add_argument("--json", type=Path, default=default_prefix.with_suffix(".json"))
    parser.add_argument("--report", type=Path, default=default_prefix.with_suffix(".md"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_rows = load_rows(args.input)
    rows = [classify_row(row) for row in input_rows]
    summary = summarize(rows)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    write_csv(rows, args.csv)
    write_json(rows, summary, args.json)
    write_markdown(rows, summary, args.report, args.input)
    print(f"wrote {args.csv}")
    print(f"wrote {args.json}")
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()
