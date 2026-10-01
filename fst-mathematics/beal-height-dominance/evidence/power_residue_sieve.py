#!/usr/bin/env python3
"""
Power-residue sieve for Beal Height-Dominance evidence.

Given near-power records A^x+B^y ~= C^z, test the necessary local condition

    A^x + B^y mod q is a z-th power residue mod q

for many small primes q.  An exact solution C^z passes every q.  Near-misses
usually fail many nontrivial verifier primes even when the real-valued gap to
C^z is tiny.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set, Tuple


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


def z_power_residues(q: int, z: int) -> Set[int]:
    return {pow(a, z, q) for a in range(q)}


def local_density(q: int, x: int, y: int, z: int, residues: Set[int]) -> float:
    accepted = 0
    q2 = q * q
    ax = [pow(a, x, q) for a in range(q)]
    by = [pow(b, y, q) for b in range(q)]
    for av in ax:
        for bv in by:
            if (av + bv) % q in residues:
                accepted += 1
    return accepted / q2


@dataclass
class RowScore:
    rank: int
    A: int
    x: int
    B: int
    y: int
    C: int
    z: int
    rel_gap: float
    nontrivial_tests: int
    passed: int
    failed: int
    fail_rate: float
    verifier_bits: float
    first_fail_q: int
    failed_qs: str
    passed_qs: str


def load_rows(path: Path, limit: int) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
            if limit and len(rows) >= limit:
                break
    return rows


def score_rows(rows: List[Dict[str, str]], args: argparse.Namespace) -> Tuple[List[RowScore], Dict[Tuple[int, int, int, int], float]]:
    primes = primes_upto(args.prime_bound)
    residue_cache: Dict[Tuple[int, int], Set[int]] = {}
    density_cache: Dict[Tuple[int, int, int, int], float] = {}
    scores: List[RowScore] = []

    for idx, row in enumerate(rows, start=1):
        A = int(row["A"])
        x = int(row["x"])
        B = int(row["B"])
        y = int(row["y"])
        C = int(row["C"])
        z = int(row["z"])
        rel_gap = float(row.get("rel_gap", "nan"))

        passed_qs: List[int] = []
        failed_qs: List[int] = []
        verifier_bits = 0.0

        for q in primes:
            d = math.gcd(z, q - 1)
            if args.only_nontrivial and d == 1:
                continue
            key = (q, z)
            if key not in residue_cache:
                residue_cache[key] = z_power_residues(q, z)
            residues = residue_cache[key]

            density_key = (q, x, y, z)
            if density_key not in density_cache:
                density_cache[density_key] = local_density(q, x, y, z, residues)
            delta = density_cache[density_key]
            if delta > 0:
                verifier_bits += -math.log(delta, 2)

            value = (pow(A, x, q) + pow(B, y, q)) % q
            if value in residues:
                passed_qs.append(q)
            else:
                failed_qs.append(q)

        tests = len(passed_qs) + len(failed_qs)
        failed = len(failed_qs)
        passed = len(passed_qs)
        scores.append(
            RowScore(
                rank=int(row.get("rank", idx)),
                A=A,
                x=x,
                B=B,
                y=y,
                C=C,
                z=z,
                rel_gap=rel_gap,
                nontrivial_tests=tests,
                passed=passed,
                failed=failed,
                fail_rate=failed / tests if tests else 0.0,
                verifier_bits=verifier_bits,
                first_fail_q=failed_qs[0] if failed_qs else 0,
                failed_qs=" ".join(map(str, failed_qs[: args.list_qs])),
                passed_qs=" ".join(map(str, passed_qs[: args.list_qs])),
            )
        )

    return scores, density_cache


def write_scores(scores: List[RowScore], path: Path) -> None:
    fields = [
        "rank",
        "A",
        "x",
        "B",
        "y",
        "C",
        "z",
        "rel_gap",
        "nontrivial_tests",
        "passed",
        "failed",
        "fail_rate",
        "verifier_bits",
        "first_fail_q",
        "failed_qs",
        "passed_qs",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for s in scores:
            writer.writerow(
                {
                    "rank": s.rank,
                    "A": s.A,
                    "x": s.x,
                    "B": s.B,
                    "y": s.y,
                    "C": s.C,
                    "z": s.z,
                    "rel_gap": f"{s.rel_gap:.12e}",
                    "nontrivial_tests": s.nontrivial_tests,
                    "passed": s.passed,
                    "failed": s.failed,
                    "fail_rate": f"{s.fail_rate:.8f}",
                    "verifier_bits": f"{s.verifier_bits:.4f}",
                    "first_fail_q": s.first_fail_q,
                    "failed_qs": s.failed_qs,
                    "passed_qs": s.passed_qs,
                }
            )


def median(values: List[float]) -> float:
    vals = sorted(values)
    n = len(vals)
    if n == 0:
        return float("nan")
    mid = n // 2
    if n % 2:
        return vals[mid]
    return 0.5 * (vals[mid - 1] + vals[mid])


def bar(count: int, total: int, width: int = 36) -> str:
    if total <= 0:
        return ""
    filled = round(width * count / total)
    return "#" * filled + "." * (width - filled)


def write_summary(scores: List[RowScore], density_cache: Dict[Tuple[int, int, int, int], float], path: Path, score_path: Path, args: argparse.Namespace) -> None:
    total = len(scores)
    zero_fail = sum(1 for s in scores if s.failed == 0)
    all_fail_rates = [s.fail_rate for s in scores]
    first_fail_counts = Counter(s.first_fail_q for s in scores if s.first_fail_q)
    z_counts = Counter(s.z for s in scores)
    fail_bins = [
        ("0", 0.0, 1e-12),
        ("(0,0.25)", 1e-12, 0.25),
        ("[0.25,0.50)", 0.25, 0.50),
        ("[0.50,0.75)", 0.50, 0.75),
        ("[0.75,1.00]", 0.75, 1.0000001),
    ]

    lines: List[str] = []
    lines.append("# Power-Residue-Sieve Summary")
    lines.append("")
    lines.append(f"- Input CSV: `{args.input.name}`")
    lines.append(f"- Rows scored: `{total}`")
    lines.append(f"- Prime bound: `{args.prime_bound}`")
    lines.append(f"- Only nontrivial verifier primes: `{args.only_nontrivial}`")
    lines.append(f"- Score CSV: `{score_path.name}`")
    lines.append("")
    lines.append("## Kerndeutung")
    lines.append("")
    lines.append("Ein exakter Treffer `A^x+B^y=C^z` muss für jede Primzahl `q` bestehen:")
    lines.append("")
    lines.append("```text")
    lines.append("A^x+B^y mod q in {u^z mod q}")
    lines.append("```")
    lines.append("")
    lines.append("Die Near-Misses sind metrisch nah an `C^z`, scheitern aber lokal an vielen kleinen Verifiern.")
    lines.append("")
    lines.append("## Hauptbefund")
    lines.append("")
    lines.append(f"- Null lokale Failures: {zero_fail}/{total} ({zero_fail/total:.2%}) `{bar(zero_fail,total)}`")
    lines.append(f"- Mindestens eine lokale Failure: {total-zero_fail}/{total} ({(total-zero_fail)/total:.2%})")
    lines.append(f"- Median Fail-Rate: `{median(all_fail_rates):.4f}`")
    lines.append(f"- Minimum Fail-Rate: `{min(all_fail_rates):.4f}`")
    lines.append(f"- Maximum Fail-Rate: `{max(all_fail_rates):.4f}`")
    lines.append(f"- Median verifier entropy: `{median([s.verifier_bits for s in scores]):.2f}` bits")
    lines.append("")
    lines.append("## Fail-Rate-Verteilung")
    lines.append("")
    for label, lo, hi in fail_bins:
        count = sum(1 for s in scores if lo <= s.fail_rate < hi)
        lines.append(f"- {label}: {count} ({count/total:.2%}) `{bar(count,total)}`")
    lines.append("")
    lines.append("## Ziel-Exponenten")
    lines.append("")
    for z, count in sorted(z_counts.items()):
        lines.append(f"- z={z}: {count} ({count/total:.2%})")
    lines.append("")
    lines.append("## Häufigste erste lokale Failure")
    lines.append("")
    for q, count in first_fail_counts.most_common(15):
        lines.append(f"- q={q}: {count} ({count/total:.2%})")
    lines.append("")
    lines.append("## 20 engste Near-Misses mit Sieve-Score")
    lines.append("")
    lines.append("| rank | A^x | B^y | C^z | rel_gap | tests | failed | fail_rate | first_fail_q | failed_qs |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for s in scores[:20]:
        lines.append(
            f"| {s.rank} | {s.A}^{s.x} | {s.B}^{s.y} | {s.C}^{s.z} | "
            f"{s.rel_gap:.3e} | {s.nontrivial_tests} | {s.failed} | "
            f"{s.fail_rate:.3f} | {s.first_fail_q} | `{s.failed_qs}` |"
        )
    lines.append("")
    lines.append("## Lokale Acceptance-Dichten")
    lines.append("")
    lines.append("Die Dichte `delta_q(x,y,z)` misst, welcher Anteil aller Paare `(a,b) mod q` den lokalen Verifier besteht.")
    lines.append("Hier sind die 20 stärksten lokalen Filter in der beobachteten Parametermenge:")
    lines.append("")
    lines.append("| q | x | y | z | delta | bits |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    strongest = sorted(density_cache.items(), key=lambda kv: kv[1])[:20]
    for (q, x, y, z), delta in strongest:
        bits = -math.log(delta, 2) if delta > 0 else float("inf")
        lines.append(f"| {q} | {x} | {y} | {z} | {delta:.5f} | {bits:.3f} |")
    lines.append("")
    lines.append("## Deutung")
    lines.append("")
    lines.append(
        "Der Power-Residue-Sieve ist die lokale kryptographische Sicht: "
        "Ein echter perfect-power Hash muss jeden kleinen Verifier bestehen. "
        "Fast-Potenzen sehen reell nah aus, fallen aber meist schon modulo kleinen "
        "Primzahlen auseinander. Eine Beal-Lösung wäre daher nicht nur ein "
        "metrischer Near-Hit, sondern ein globaler Verifier-Pass über alle q."
    )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parent / "near_power_b250_e3-7_z3-9_k5000.csv")
    parser.add_argument("--limit", type=int, default=5000)
    parser.add_argument("--prime-bound", type=int, default=251)
    parser.add_argument("--list-qs", type=int, default=20)
    parser.add_argument("--include-trivial", action="store_true")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.only_nontrivial = not args.include_trivial
    return args


def main() -> None:
    args = parse_args()
    rows = load_rows(args.input, args.limit)
    scores, density_cache = score_rows(rows, args)
    stem = f"power_residue_sieve_q{args.prime_bound}_n{len(scores)}"
    if not args.only_nontrivial:
        stem += "_allq"
    score_path = args.out_dir / f"{stem}.csv"
    summary_path = args.out_dir / f"{stem}_summary.md"
    write_scores(scores, score_path)
    write_summary(scores, density_cache, summary_path, score_path, args)
    print(f"Scored {len(scores)} records.")
    print(f"Scores: {score_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
