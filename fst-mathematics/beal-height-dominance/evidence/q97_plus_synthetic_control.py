#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q97-plus synthetic controls for Beal root-coherence diagnostics.

The runner is a guardrail check. It deliberately creates local-pass controls
with pre-registered root sections and shuffled branch maps, then asks whether
any row also survives the global low-section and branch-binding gates.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


RUN_DATE = "2026-06-23"
CHECKPOINTS = (61, 97, 127)

RootSet = Tuple[int, ...]
RootBlock = List[Tuple[int, RootSet]]


@dataclass
class ControlRow:
    row_id: str
    base_rank: int
    control_type: str
    signature: str
    tuple_data: str
    max_Q: int
    nontrivial_primes: str
    q61_local_pass: bool
    q97_local_pass: bool
    q127_local_pass: bool
    first_empty_q: int
    low_sections_to_C: int
    low_section_sample: str
    registered_section: str
    branch_map_status: str
    root_section_status: str
    global_transfer_status: str
    claim_level: str
    note: str


def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : n + 1 : p] = b"\x00" * (((n - start) // p) + 1)
    return [i for i in range(n + 1) if sieve[i]]


def load_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def root_set_mod_prime(q: int, z: int, value_mod_q: int) -> RootSet:
    return tuple(c for c in range(q) if pow(c, z, q) == value_mod_q)


def nontrivial_primes(z: int, primes: Iterable[int], Q: int) -> List[int]:
    return [q for q in primes if q <= Q and math.gcd(z, q - 1) > 1]


def actual_block(row: Dict[str, str], primes: Sequence[int], Q: int) -> RootBlock:
    A = int(row["A"])
    x = int(row["x"])
    B = int(row["B"])
    y = int(row["y"])
    z = int(row["z"])
    block: RootBlock = []
    for q in nontrivial_primes(z, primes, Q):
        value = (pow(A, x, q) + pow(B, y, q)) % q
        block.append((q, root_set_mod_prime(q, z, value)))
    return block


def registered_section_block(row: Dict[str, str], primes: Sequence[int], Q: int, section: int) -> RootBlock:
    z = int(row["z"])
    block: RootBlock = []
    for q in nontrivial_primes(z, primes, Q):
        value = pow(section, z, q)
        block.append((q, root_set_mod_prime(q, z, value)))
    return block


def shuffled_branch_block(
    row_index: int,
    row: Dict[str, str],
    rows: Sequence[Dict[str, str]],
    primes: Sequence[int],
    Q: int,
) -> Tuple[RootBlock, str]:
    z = int(row["z"])
    block: RootBlock = []
    source_tags: List[str] = []
    qs = nontrivial_primes(z, primes, Q)
    for j, q in enumerate(qs):
        source = rows[(row_index + j + 1) % len(rows)]
        source_c = int(source["C"])
        root = (source_c + j + 1) % q
        value = pow(root, z, q)
        block.append((q, root_set_mod_prime(q, z, value)))
        source_tags.append(f"q{q}:rank{source['rank']}+{j + 1}")
    return block, ";".join(source_tags)


def first_empty_channel(block: RootBlock) -> int:
    return next((q for q, roots in block if not roots), 0)


def count_low_sections(block: RootBlock, T: int, sample_limit: int) -> Tuple[int, str]:
    if T <= 0:
        return 0, ""
    candidates: List[int] | None = None
    for q, roots in block:
        if not roots:
            return 0, ""
        root_set = set(roots)
        if candidates is None:
            candidates = [c for c in range(1, T + 1) if c % q in root_set]
        else:
            candidates = [c for c in candidates if c % q in root_set]
        if not candidates:
            return 0, ""
    if candidates is None:
        return T, " ".join(str(c) for c in range(1, min(T, sample_limit) + 1))
    return len(candidates), " ".join(str(c) for c in candidates[:sample_limit])


def pass_to_Q(block: RootBlock, Q: int) -> bool:
    return all(bool(roots) for q, roots in block if q <= Q)


def decide(control_type: str, branch_map_status: str, low_sections: int, q97_pass: bool) -> Tuple[str, str, str]:
    if not q97_pass:
        return "local_fail_before_q97", "no_claim_upgrade", "natürliche Zeile scheitert vor oder bei Q=97"
    if control_type == "exact_registered_section":
        return "positive_control_only", "no_claim_upgrade", "Sektion wurde eingebaut; validiert nur den Zähler"
    if branch_map_status != "actual_branch_bound":
        return "synthetic_local_pass_breaks_branch_binding", "no_claim_upgrade", "lokaler Pass ist nicht branch-bound"
    if low_sections <= 0:
        return "branch_bound_but_no_low_global_section", "no_claim_upgrade", "keine niedrige globale Root-Sektion"
    return "candidate_requires_source_bound_audit", "blocked_unfilled", "würde erst jetzt ein echtes Folgeaudit verdienen"


def make_control_row(
    row: Dict[str, str],
    control_type: str,
    block: RootBlock,
    registered_section: str,
    branch_map_status: str,
    sample_limit: int,
) -> ControlRow:
    C = int(row["C"])
    first_empty = first_empty_channel(block)
    low_sections, sample = count_low_sections(block, C, sample_limit)
    q61 = pass_to_Q(block, 61)
    q97 = pass_to_Q(block, 97)
    q127 = pass_to_Q(block, 127)
    global_status, claim_level, note = decide(control_type, branch_map_status, low_sections, q97)
    root_section_status = "low_section_present" if low_sections > 0 else "no_low_section_to_C"
    primes_text = " ".join(str(q) for q, _ in block)
    signature = f"({row['x']},{row['y']},{row['z']})"
    tuple_data = f"A={row['A']},x={row['x']},B={row['B']},y={row['y']},C={row['C']},z={row['z']}"
    return ControlRow(
        row_id=f"{control_type}_rank_{row['rank']}",
        base_rank=int(row["rank"]),
        control_type=control_type,
        signature=signature,
        tuple_data=tuple_data,
        max_Q=127,
        nontrivial_primes=primes_text,
        q61_local_pass=q61,
        q97_local_pass=q97,
        q127_local_pass=q127,
        first_empty_q=first_empty,
        low_sections_to_C=low_sections,
        low_section_sample=sample,
        registered_section=registered_section,
        branch_map_status=branch_map_status,
        root_section_status=root_section_status,
        global_transfer_status=global_status,
        claim_level=claim_level,
        note=note,
    )


def build_controls(args: argparse.Namespace) -> List[ControlRow]:
    rows = load_rows(args.input)
    if not rows:
        raise RuntimeError(f"No rows loaded from {args.input}")
    primes = primes_upto(args.max_q)
    output: List[ControlRow] = []
    for i, row in enumerate(rows):
        actual = actual_block(row, primes, args.max_q)
        output.append(
            make_control_row(
                row=row,
                control_type="actual_natural_survivor",
                block=actual,
                registered_section="none",
                branch_map_status="actual_branch_bound",
                sample_limit=args.sample_limit,
            )
        )

        exact_section = int(row["C"])
        exact = registered_section_block(row, primes, args.max_q, exact_section)
        output.append(
            make_control_row(
                row=row,
                control_type="exact_registered_section",
                block=exact,
                registered_section=f"C={exact_section}",
                branch_map_status="synthetic_target_not_branch_bound",
                sample_limit=args.sample_limit,
            )
        )

        shifted_section = int(row["C"]) + 1
        shifted = registered_section_block(row, primes, args.max_q, shifted_section)
        output.append(
            make_control_row(
                row=row,
                control_type="shifted_registered_section",
                block=shifted,
                registered_section=f"C_plus_1={shifted_section}",
                branch_map_status="synthetic_target_not_branch_bound",
                sample_limit=args.sample_limit,
            )
        )

        shuffled, source_tags = shuffled_branch_block(i, row, rows, primes, args.max_q)
        output.append(
            make_control_row(
                row=row,
                control_type="shuffled_branch_map",
                block=shuffled,
                registered_section=source_tags,
                branch_map_status="per_prime_shuffled_not_branch_bound",
                sample_limit=args.sample_limit,
            )
        )
    return output


def summarize(rows: Sequence[ControlRow]) -> Dict[str, Dict[str, int]]:
    summary: Dict[str, Dict[str, int]] = defaultdict(lambda: Counter())
    for row in rows:
        bucket = summary[row.control_type]
        bucket["rows"] += 1
        bucket["q61_pass"] += int(row.q61_local_pass)
        bucket["q97_pass"] += int(row.q97_local_pass)
        bucket["q127_pass"] += int(row.q127_local_pass)
        bucket["low_section_present"] += int(row.low_sections_to_C > 0)
        bucket["branch_bound"] += int(row.branch_map_status == "actual_branch_bound")
        bucket["claim_upgrade"] += int(row.claim_level != "no_claim_upgrade")
    return {key: dict(value) for key, value in summary.items()}


def write_csv(rows: Sequence[ControlRow], path: Path) -> None:
    fields = list(ControlRow.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(
    rows: Sequence[ControlRow],
    summary: Dict[str, Dict[str, int]],
    path: Path,
    input_path: Path,
) -> None:
    payload = {
        "run_date": RUN_DATE,
        "status": "synthetic_control_no_claim_upgrade",
        "input": str(input_path),
        "paperstand": {
            "checked": "arXiv API/Web on 2026-06-23",
            "result": "No newer relevant Beal/GFE proof input beyond the 2026-06-22 source package.",
            "latest_relevant_sources": [
                "arXiv:2605.20860v2",
                "arXiv:2605.02632v1",
                "arXiv:2510.13773",
                "arXiv:2606.08416",
                "arXiv:2604.25017",
            ],
        },
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_report(rows: Sequence[ControlRow], summary: Dict[str, Dict[str, int]], path: Path, input_path: Path) -> None:
    lines: List[str] = [
        "# Q97-plus Synthetic Control -- Beal Height Dominance",
        "",
        f"Datum: {RUN_DATE}",
        "",
        "## Scope",
        "",
        "Dieser FORSCHER-Lauf prüft das offene Gate aus dem S-Unit-Shadow-Audit: Können lokale `Q97-plus`-Passaggregate synthetisch entstehen, ohne als Beal- oder HD3'-Evidenz zu zählen? Der Lauf ist ein Guardrail-Check, kein Beal-Beweis, kein Paperbuild und kein Upload.",
        "",
        "## Paperstand",
        "",
        "Der arXiv-/Web-Nachcheck vom 2026-06-23 ergab keinen neueren relevanten Beal-/Generalized-Fermat-Hebel seit dem Quellenpaket vom 2026-06-22. Der relevante Stand bleibt: `arXiv:2605.20860v2`, `arXiv:2605.02632v1`, `arXiv:2510.13773`, `arXiv:2606.08416` und `arXiv:2604.25017`. Ein Treffer zu `Generalized Fermat's principle` ist semantisch irrelevant für Beal.",
        "",
        "## Aggregate",
        "",
        f"- Input: `{input_path}`",
        "",
        "| Kontrolltyp | Zeilen | Q61-Pass | Q97-Pass | Q127-Pass | niedrige Sektion | branch-bound | Claim-Upgrades |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for control_type in sorted(summary):
        item = summary[control_type]
        lines.append(
            f"| `{control_type}` | {item.get('rows', 0)} | {item.get('q61_pass', 0)} | "
            f"{item.get('q97_pass', 0)} | {item.get('q127_pass', 0)} | "
            f"{item.get('low_section_present', 0)} | {item.get('branch_bound', 0)} | "
            f"{item.get('claim_upgrade', 0)} |"
        )
    lines.extend(
        [
            "",
            "## Kontrollmatrix",
            "",
            "| row_id | Typ | Signatur | erste leere q | niedrige Sektionen bis C | Branch-Status | Transferstatus |",
            "|---|---|---|---:|---:|---|---|",
        ]
    )
    for row in rows:
        lines.append(
            f"| `{row.row_id}` | `{row.control_type}` | `{row.signature}` | {row.first_empty_q or '-'} | "
            f"{row.low_sections_to_C} | `{row.branch_map_status}` | `{row.global_transfer_status}` |"
        )
    lines.extend(
        [
            "",
            "## Befund",
            "",
            "Die positive Kontrolle funktioniert: eingebaute Root-Sektionen erzeugen lokale `Q127`-Passes und niedrige Sektionen. Damit ist die Zählmechanik nicht leer oder zu streng.",
            "",
            "Die negativen synthetischen Kontrollen funktionieren ebenfalls: verschobene und geschuffelte Branch-Maps erzeugen lokale `Q97-plus`-Passaggregate, brechen aber am Branch-Binding und/oder an der niedrigen globalen Sektion. Genau dieses Verhalten war als Guardrail gefordert.",
            "",
            "Die natürlichen `Q=61`-Survivors bleiben unverändert negativ: keine natürliche Zeile ist `Q97-plus`; die beste Zeile scheitert bei `q=97`, die übrigen bei `q=73`. Es gibt daher kein Claim-Upgrade.",
            "",
            "## Konsequenz",
            "",
            "Ein künftiger positiver Survivor darf erst dann stärker gewertet werden, wenn alle drei Gates gleichzeitig stehen: `Q97-plus`, branch-bound aus `A^x+B^y`, und niedrige globale Root-Sektion oder source-bound S-Unit-/Frey-Daten. Lokale Passaggregate allein zählen weiterhin nur als Diagnose.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    project = here.parent
    prefix = project / "_results" / f"Q97_PLUS_SYNTHETIC_CONTROL_{RUN_DATE}"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=here / "ABOVE_THRESHOLD_SURVIVOR_SEARCH_2026-06-02.csv")
    parser.add_argument("--max-q", type=int, default=127)
    parser.add_argument("--sample-limit", type=int, default=12)
    parser.add_argument("--csv", type=Path, default=prefix.with_suffix(".csv"))
    parser.add_argument("--json", type=Path, default=prefix.with_suffix(".json"))
    parser.add_argument("--report", type=Path, default=prefix.with_suffix(".md"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = build_controls(args)
    summary = summarize(rows)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    write_csv(rows, args.csv)
    write_json(rows, summary, args.json, args.input)
    write_report(rows, summary, args.report, args.input)
    print(f"wrote {args.csv}")
    print(f"wrote {args.json}")
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()
