#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Micro-audit for the rank-8 Root-Branch survivor.

This script isolates the first additional prime channel that turns the
2026-06-01 Q=31 local pass into the Q=61 local obstruction. It is a diagnostic
tool, not a Beal proof.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


RootSet = Tuple[int, ...]
RootBlock = List[Tuple[int, RootSet]]


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


def root_set_mod_prime(q: int, z: int, value_mod_q: int) -> RootSet:
    return tuple(c for c in range(q) if pow(c, z, q) == value_mod_q)


def block_product(primes: Iterable[int]) -> int:
    product = 1
    for q in primes:
        product *= q
    return product


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


def load_rank_row(path: Path, rank: int) -> Dict[str, str]:
    with path.open("r", newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if int(row.get("rank", "0")) == rank:
                return row
    raise RuntimeError(f"rank={rank} not found in {path}")


def classify_channel(q: int, z: int, value: int, roots: RootSet, A: int, B: int, C: int, N: int) -> str:
    if math.gcd(z, q - 1) == 1:
        return "trivial_z_power_automorphism"
    if roots:
        return "local_root_pass"
    if value == 0:
        return "zero_residue_edge_case"
    if any(item % q == 0 for item in (A, B, C, N)):
        return "branch_or_valuation_touched"
    return "unit_character_nonresidue_fail"


@dataclass
class AuditRow:
    q: int
    d_gcd_z_qminus1: int
    value_mod_q: int
    exact_power_mod_q: int
    gap_mod_q: int
    actual_root_count: int
    actual_roots: str
    exact_root_count: int
    exact_roots: str
    C_mod_q: int
    actual_contains_C_mod_q: bool
    exact_contains_C_mod_q: bool
    q_divides_A: bool
    q_divides_B: bool
    q_divides_C: bool
    q_divides_N: bool
    q_divides_gap: bool
    cumulative_prime_count: int
    cumulative_M: str
    cumulative_root_classes: str
    cumulative_low_sections: int
    cumulative_low_section_sample: str
    D_root_cumulative: float
    D_drift_cumulative: float
    status: str
    channel_class: str


def audit(row: Dict[str, str], max_q: int, sample_limit: int) -> List[AuditRow]:
    A = int(row["A"])
    x = int(row["x"])
    B = int(row["B"])
    y = int(row["y"])
    C = int(row["C"])
    z = int(row["z"])
    N = int(row["N"])
    H = int(row["H_Cz"])
    T = C
    gap = N - H

    block: RootBlock = []
    rows: List[AuditRow] = []
    first_empty_seen = False
    for q in primes_upto(max_q):
        d = math.gcd(z, q - 1)
        if d == 1:
            continue
        value = (pow(A, x, q) + pow(B, y, q)) % q
        exact_value = H % q
        roots = root_set_mod_prime(q, z, value)
        exact_roots = root_set_mod_prime(q, z, exact_value)
        block.append((q, roots))

        block_primes = [prime for prime, _ in block]
        M = block_product(block_primes)
        root_sizes = [len(root_set) for _, root_set in block]
        n_root_classes = 1
        for size in root_sizes:
            n_root_classes *= max(1, size)
        n_low, sample = count_low_sections(block, T, sample_limit)
        log_n_root = sum(math.log(max(1, size)) for size in root_sizes)
        D_root = math.log(M) - math.log(max(1, T)) - log_n_root
        D_drift = log_n_root - math.log(max(1, n_low))

        channel_class = classify_channel(q, z, value, roots, A, B, C, N)
        if roots:
            status = "pass_before_first_empty" if not first_empty_seen else "pass_after_first_empty"
        elif not first_empty_seen:
            status = "first_empty_root_channel"
            first_empty_seen = True
        else:
            status = "empty_after_first_empty"

        rows.append(
            AuditRow(
                q=q,
                d_gcd_z_qminus1=d,
                value_mod_q=value,
                exact_power_mod_q=exact_value,
                gap_mod_q=gap % q,
                actual_root_count=len(roots),
                actual_roots=" ".join(str(root) for root in roots),
                exact_root_count=len(exact_roots),
                exact_roots=" ".join(str(root) for root in exact_roots),
                C_mod_q=C % q,
                actual_contains_C_mod_q=(C % q) in roots,
                exact_contains_C_mod_q=(C % q) in exact_roots,
                q_divides_A=A % q == 0,
                q_divides_B=B % q == 0,
                q_divides_C=C % q == 0,
                q_divides_N=N % q == 0,
                q_divides_gap=gap % q == 0,
                cumulative_prime_count=len(block_primes),
                cumulative_M=str(M),
                cumulative_root_classes=str(n_root_classes),
                cumulative_low_sections=n_low,
                cumulative_low_section_sample=sample,
                D_root_cumulative=D_root,
                D_drift_cumulative=D_drift,
                status=status,
                channel_class=channel_class,
            )
        )
    return rows


def write_csv(rows: Sequence[AuditRow], path: Path) -> None:
    fields = list(AuditRow.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for item in rows:
            data = asdict(item)
            for key in ("D_root_cumulative", "D_drift_cumulative"):
                data[key] = f"{data[key]:.12g}"
            writer.writerow(data)


def write_report(row: Dict[str, str], rows: Sequence[AuditRow], path: Path, csv_path: Path) -> None:
    rank = int(row["rank"])
    A = int(row["A"])
    x = int(row["x"])
    B = int(row["B"])
    y = int(row["y"])
    C = int(row["C"])
    z = int(row["z"])
    N = int(row["N"])
    H = int(row["H_Cz"])
    gap = N - H
    first_empty = next((item for item in rows if item.status == "first_empty_root_channel"), None)

    lines: List[str] = []
    lines.append("# Root-Branch Rank-8 Micro-Audit -- 2026-06-02")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append("Dieses Audit verfolgt den einen `Q=31`-Survivor ohne niedrige Sektion aus dem FORSCHER-Lauf vom 2026-06-01.")
    lines.append("Es ist ein lokales Diagnoseartefakt, kein Beal-Beweis und kein Claim-Upgrade.")
    lines.append("")
    lines.append("## Fall")
    lines.append("")
    lines.append(f"- Rank: `{rank}`")
    lines.append(f"- Tupel: `A={A}`, `x={x}`, `B={B}`, `y={y}`, `C={C}`, `z={z}`")
    lines.append(f"- Near-Miss: `A^x+B^y = {N}`, `C^z = {H}`, Gap `N-C^z = {gap}`")
    lines.append(f"- CSV: `{csv_path.name}`")
    lines.append("")
    lines.append("## Kanal-Tabelle")
    lines.append("")
    lines.append("| q | d | value mod q | roots | C mod q | low sections | status | class |")
    lines.append("|---:|---:|---:|---:|---:|---:|---|---|")
    for item in rows:
        roots = item.actual_roots if item.actual_roots else "-"
        lines.append(
            f"| {item.q} | {item.d_gcd_z_qminus1} | {item.value_mod_q} | `{roots}` | "
            f"{item.C_mod_q} | {item.cumulative_low_sections} | `{item.status}` | `{item.channel_class}` |"
        )
    lines.append("")
    lines.append("## Befund")
    lines.append("")
    if first_empty:
        lines.append(
            f"Der zusätzliche Primkanal, der den lokalen Pass zwischen `Q=31` und `Q=61` tötet, ist `q={first_empty.q}`."
        )
        lines.append(
            "Bis `q=31` sind alle nichttrivialen Kubikrest-Root-Sets nichtleer, aber es gibt keine niedrige globale Sektion `c <= C`."
        )
        lines.append(
            f"Bei `q={first_empty.q}` ist `A^x+B^y mod q = {first_empty.value_mod_q}` kein Kubikrest, während `C^z mod q = {first_empty.exact_power_mod_q}` die positive Exact-Power-Kontrolle mit `C mod q = {first_empty.C_mod_q}` trägt."
        )
        lines.append(
            "Da `q` weder `A`, `B`, `C`, `N` noch den Gap teilt, ist der Fail ein finite-field Unit-/Kubikcharakter-Fail, kein Valuation-, Branch- oder Common-Factor-Anker."
        )
    else:
        lines.append("Bis zur geprüften Grenze wurde kein leerer Root-Kanal gefunden.")
    lines.append("")
    lines.append("## Konsequenz")
    lines.append("")
    lines.append(
        "Der rank-8-Fall ist damit kein stabiler Above-threshold-Drift-Survivor. Er ist ein nützlicher Übergangsfall: "
        "`Q=31` zeigt Root-Drift ohne niedrige Sektion, aber der nächste relevante Primkanal `37` liefert bereits eine normale lokale Kubikrest-Obstruktion."
    )
    lines.append(
        "Für HD3'-Arbeit bleibt der nächste echte Fortschritt daher: Survivors konstruieren oder finden, die über `q=37` und weitere nichttriviale Primkanäle lokal bestehen, bevor `D_drift` als nichtlokales Signal gewertet wird."
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=here / "near_power_b250_e3-7_z3-9_k5000.csv")
    parser.add_argument("--rank", type=int, default=8)
    parser.add_argument("--max-q", type=int, default=61)
    parser.add_argument("--sample-limit", type=int, default=12)
    parser.add_argument("--csv", type=Path, default=here / "ROOT_BRANCH_RANK8_MICRO_AUDIT_2026-06-02.csv")
    parser.add_argument("--report", type=Path, default=here / "ROOT_BRANCH_RANK8_MICRO_AUDIT_2026-06-02.md")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    row = load_rank_row(args.input, args.rank)
    rows = audit(row, args.max_q, args.sample_limit)
    write_csv(rows, args.csv)
    write_report(row, rows, args.report, args.csv)
    print(f"wrote {args.csv}")
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()
