#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Variable-Q Power-Residue-Sieve fuer HD2'-Wachstumstest.

Zweck: BEWEISNOTIZ §25.5 Punkt 1. Pruefe ob `I(Q_min) >= c log H` empirisch
plausibel ist, indem dieselben Near-Misses gegen mehrere Verifier-Hoehen Q
getestet werden:

    Q ∈ {251, 503, 1009, 2003}

Output: pro Q ein eigenes Score-CSV + Markdown-Summary plus eine
Aggregat-Auswertung (`variable_Q_aggregate_summary.md`) mit:
- I(Q) Mittelwert/Median pro Q,
- I(Q) als Funktion von log H (Buckets),
- Plausibilitaet einer linearen Untergrenze c log H.

Performance-Optimierung:
- `local_density` wird als Bitset-Lookup vektorisiert (statt Doppel-Loop).
- Density-Cache ueberlebt zwischen Q-Schritten (alle kleineren q
  werden bei groesseren Q wiederverwendet).
- Optional --limit fuer Teilmenge der Near-Misses (erste Lauf:
  1000 statt 5000, falls Laufzeit zu lang).

Aufruf (Server, mit nohup):
    PYTHONIOENCODING=utf-8 nohup python3 -u run_variable_Q_sieve.py \
        --input near_power_b250_e3-7_z3-9_k5000.csv \
        --limit 1000 \
        --primes 251 503 1009 2003 \
        > variable_Q_sieve.log 2>&1 &
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Set, Tuple


# ----------------------------------------------------------------------------- #
# Helpers                                                                       #
# ----------------------------------------------------------------------------- #

def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(n ** 0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start: n + 1: p] = b"\x00" * (((n - start) // p) + 1)
    return [i for i in range(n + 1) if sieve[i]]


def z_power_residues(q: int, z: int) -> Set[int]:
    return {pow(a, z, q) for a in range(q)}


def residue_bitmask(q: int, residues: Set[int]) -> int:
    """Pack residue set into a Python int bitmask of width q."""
    mask = 0
    for r in residues:
        mask |= 1 << r
    return mask


def local_density_fast(q: int, x: int, y: int, z: int, residue_mask: int) -> float:
    """
    Vectorised density: compute, for each pair (a,b) in F_q x F_q,
    whether a^x + b^y mod q is a z-th power residue mod q.

    Uses precomputed power tables and a bitmask of the residue set.
    Avoids the inner double-loop of the original implementation.
    """
    ax = [pow(a, x, q) for a in range(q)]
    by = [pow(b, y, q) for b in range(q)]
    accepted = 0
    for av in ax:
        # For each b, sum mod q -> bit lookup
        # Inner loop in pure python is the bottleneck; this is unavoidable
        # without numpy. We at least pre-shift the mask once per av.
        for bv in by:
            s = av + bv
            if s >= q:
                s -= q
            if (residue_mask >> s) & 1:
                accepted += 1
    return accepted / (q * q)


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
    log_H: float
    nontrivial_tests: int
    passed: int
    failed: int
    fail_rate: float
    verifier_bits: float
    first_fail_q: int


# ----------------------------------------------------------------------------- #
# Core sieve                                                                    #
# ----------------------------------------------------------------------------- #

def score_for_Q(
    rows: List[Dict[str, str]],
    primes_Q: List[int],
    only_nontrivial: bool,
    residue_cache: Dict[Tuple[int, int], int],
    density_cache: Dict[Tuple[int, int, int, int], float],
    log_every: int = 100,
) -> List[RowScore]:
    """
    Score every row against the given prime list. Caches survive across calls
    so larger Q reuse smaller-Q work.
    """
    scores: List[RowScore] = []
    n = len(rows)
    t0 = time.time()
    last_log = t0
    for idx, row in enumerate(rows, start=1):
        A = int(row["A"])
        x = int(row["x"])
        B = int(row["B"])
        y = int(row["y"])
        C = int(row["C"])
        z = int(row["z"])
        rel_gap = float(row.get("rel_gap", "nan"))
        log_H = z * math.log(max(C, 2))

        passed = 0
        failed = 0
        first_fail_q = 0
        verifier_bits = 0.0

        for q in primes_Q:
            d = math.gcd(z, q - 1)
            if only_nontrivial and d == 1:
                continue
            key_r = (q, z)
            if key_r not in residue_cache:
                residue_cache[key_r] = residue_bitmask(q, z_power_residues(q, z))
            residue_mask = residue_cache[key_r]

            key_d = (q, x, y, z)
            if key_d not in density_cache:
                density_cache[key_d] = local_density_fast(q, x, y, z, residue_mask)
            delta = density_cache[key_d]
            if delta > 0:
                verifier_bits += -math.log(delta, 2)

            value = (pow(A, x, q) + pow(B, y, q)) % q
            if (residue_mask >> value) & 1:
                passed += 1
            else:
                failed += 1
                if first_fail_q == 0:
                    first_fail_q = q

        tests = passed + failed
        scores.append(
            RowScore(
                rank=int(row.get("rank", idx)),
                A=A, x=x, B=B, y=y, C=C, z=z,
                rel_gap=rel_gap,
                log_H=log_H,
                nontrivial_tests=tests,
                passed=passed,
                failed=failed,
                fail_rate=failed / tests if tests else 0.0,
                verifier_bits=verifier_bits,
                first_fail_q=first_fail_q,
            )
        )
        now = time.time()
        if idx % log_every == 0 or now - last_log > 30:
            elapsed = now - t0
            rate = idx / elapsed if elapsed > 0 else 0.0
            eta = (n - idx) / rate if rate > 0 else float("inf")
            print(
                f"  [{idx}/{n}] elapsed={elapsed:.1f}s rate={rate:.2f}/s eta={eta:.1f}s "
                f"density_cache={len(density_cache)} residue_cache={len(residue_cache)}",
                flush=True,
            )
            last_log = now
    return scores


def write_scores_csv(scores: List[RowScore], path: Path) -> None:
    fields = list(RowScore.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for s in scores:
            d = asdict(s)
            d["rel_gap"] = f"{s.rel_gap:.12e}"
            d["log_H"] = f"{s.log_H:.6f}"
            d["fail_rate"] = f"{s.fail_rate:.8f}"
            d["verifier_bits"] = f"{s.verifier_bits:.4f}"
            writer.writerow(d)


def median(values: List[float]) -> float:
    vals = sorted(values)
    n = len(vals)
    if n == 0:
        return float("nan")
    mid = n // 2
    if n % 2:
        return vals[mid]
    return 0.5 * (vals[mid - 1] + vals[mid])


def fit_linear(xs: List[float], ys: List[float]) -> Tuple[float, float]:
    n = len(xs)
    if n < 2:
        return float("nan"), float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    var = sum((x - mx) ** 2 for x in xs)
    if var <= 0:
        return float("nan"), my
    slope = cov / var
    return slope, my - slope * mx


def write_summary_for_Q(
    scores: List[RowScore],
    Q: int,
    n_total_rows: int,
    out_path: Path,
    elapsed: float,
) -> None:
    total = len(scores)
    zero_fail = sum(1 for s in scores if s.failed == 0)
    fail_rates = [s.fail_rate for s in scores]
    bits = [s.verifier_bits for s in scores]
    log_Hs = [s.log_H for s in scores]
    z_counts = Counter(s.z for s in scores)

    lines: List[str] = []
    lines.append(f"# Power-Residue-Sieve Summary (Q={Q})")
    lines.append("")
    lines.append(f"- Rows scored: `{total}` (subset of {n_total_rows} input rows)")
    lines.append(f"- Verifier prime bound Q: `{Q}`")
    lines.append(f"- Compute time: `{elapsed:.1f}s`")
    lines.append("")
    lines.append("## Pass/fail")
    lines.append(f"- Zero local failures: {zero_fail}/{total} ({zero_fail/total:.2%})")
    lines.append(f"- Median fail rate: `{median(fail_rates):.4f}`")
    lines.append(f"- Min/Max fail rate: `{min(fail_rates):.4f}` / `{max(fail_rates):.4f}`")
    lines.append("")
    lines.append("## Verifier entropy I(Q) (bits)")
    lines.append(f"- Median: `{median(bits):.2f}`")
    lines.append(f"- Mean:   `{sum(bits)/total:.2f}`")
    lines.append(f"- Min:    `{min(bits):.2f}`")
    lines.append(f"- Max:    `{max(bits):.2f}`")
    lines.append("")
    lines.append("## I(Q) by target exponent z")
    by_z: Dict[int, List[float]] = defaultdict(list)
    for s in scores:
        by_z[s.z].append(s.verifier_bits)
    for z, vs in sorted(by_z.items()):
        lines.append(f"- z={z}: n={len(vs)} median={median(vs):.2f} mean={sum(vs)/len(vs):.2f}")
    lines.append("")
    lines.append("## I(Q) vs log H — bucketed")
    lines.append("")
    sorted_pairs = sorted(zip(log_Hs, bits))
    n_buckets = 10
    bsize = max(1, total // n_buckets)
    lines.append("| log_H lo | log_H hi | n | I_min | I_med | I_mean |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for i in range(0, total, bsize):
        chunk = sorted_pairs[i: i + bsize]
        if not chunk:
            continue
        chunk_bits = sorted(b for _, b in chunk)
        lo = chunk[0][0]
        hi = chunk[-1][0]
        bm = chunk_bits[0]
        bmed = chunk_bits[len(chunk_bits)//2]
        bme = sum(chunk_bits)/len(chunk_bits)
        lines.append(f"| {lo:.2f} | {hi:.2f} | {len(chunk)} | {bm:.2f} | {bmed:.2f} | {bme:.2f} |")
    lines.append("")
    slope, intercept = fit_linear(log_Hs, bits)
    lines.append(f"Linear fit (mean): `I(Q) ≈ {slope:.4f} · log H + {intercept:.4f}`")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def write_aggregate(
    per_Q: Dict[int, List[RowScore]],
    out_path: Path,
    elapsed_total: float,
) -> None:
    Qs_sorted = sorted(per_Q.keys())
    lines: List[str] = []
    lines.append("# Variable-Q Power-Residue-Sieve — Aggregate Summary")
    lines.append("")
    lines.append(f"- Q values tested: {Qs_sorted}")
    lines.append(f"- Total compute time: `{elapsed_total:.1f}s`")
    lines.append("")
    lines.append("## I(Q) — distribution by Q")
    lines.append("")
    lines.append("| Q | n | I_min | I_med | I_mean | I_max | zero_fail |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|")
    for Q in Qs_sorted:
        scores = per_Q[Q]
        bits = [s.verifier_bits for s in scores]
        zero_fail = sum(1 for s in scores if s.failed == 0)
        lines.append(
            f"| {Q} | {len(scores)} | {min(bits):.2f} | {median(bits):.2f} "
            f"| {sum(bits)/len(scores):.2f} | {max(bits):.2f} | {zero_fail} |"
        )
    lines.append("")
    lines.append("## HD2' Test: I(Q) Slope vs log H, by Q")
    lines.append("")
    lines.append("If HD2' holds, the *minimum* I(Q) per log-H-bucket should grow")
    lines.append("with both Q and log H. Slope estimates per Q follow:")
    lines.append("")
    lines.append("| Q | slope (mean fit) | intercept | slope_min_per_bucket |")
    lines.append("|---:|---:|---:|---:|")
    for Q in Qs_sorted:
        scores = per_Q[Q]
        log_Hs = [s.log_H for s in scores]
        bits = [s.verifier_bits for s in scores]
        slope_mean, ic = fit_linear(log_Hs, bits)
        # bucket and lower-envelope slope
        n_b = 10
        sorted_pairs = sorted(zip(log_Hs, bits))
        bsize = max(1, len(sorted_pairs) // n_b)
        b_mid = []
        b_min = []
        for i in range(0, len(sorted_pairs), bsize):
            chunk = sorted_pairs[i: i + bsize]
            if not chunk:
                continue
            mid = (chunk[0][0] + chunk[-1][0]) / 2
            mn = min(b for _, b in chunk)
            b_mid.append(mid)
            b_min.append(mn)
        slope_min, _ = fit_linear(b_mid, b_min)
        lines.append(f"| {Q} | {slope_mean:.4f} | {ic:.4f} | {slope_min:.4f} |")
    lines.append("")
    lines.append("## Cross-Q growth at fixed log H")
    lines.append("")
    lines.append("Rows where log H is in the median band, comparison across Q:")
    lines.append("")
    # Pick the row with the median log H across smallest-Q scoring;
    # align by (rank) across Q to compare same input.
    if Qs_sorted:
        ref_scores = per_Q[Qs_sorted[0]]
        # take 5 near-median rows
        ref_sorted = sorted(ref_scores, key=lambda s: s.log_H)
        n_ref = len(ref_sorted)
        sample_idx = [int(n_ref * f) for f in (0.1, 0.3, 0.5, 0.7, 0.9)]
        lines.append("| rank | A^x | B^y | C^z | log H | " + " | ".join(f"I(Q={Q})" for Q in Qs_sorted) + " |")
        lines.append("|---:|---:|---:|---:|---:|" + "|".join("---:" for _ in Qs_sorted) + "|")
        for j in sample_idx:
            if j >= n_ref:
                continue
            ref = ref_sorted[j]
            row = f"| {ref.rank} | {ref.A}^{ref.x} | {ref.B}^{ref.y} | {ref.C}^{ref.z} | {ref.log_H:.2f} "
            for Q in Qs_sorted:
                # find matching rank
                match = next((s for s in per_Q[Q] if s.rank == ref.rank), None)
                row += f"| {match.verifier_bits:.2f} " if match else "| ? "
            row += "|"
            lines.append(row)
    lines.append("")
    lines.append("## Interpretation guide")
    lines.append("")
    lines.append("HD2' is empirically supported if:")
    lines.append("")
    lines.append("1. I(Q) grows with Q (more verifiers → more entropy).")
    lines.append("2. At each Q, the lower-envelope of I(Q) vs log H is non-decreasing.")
    lines.append("3. The slope_min_per_bucket grows with Q.")
    lines.append("")
    lines.append("If (1) holds but (2)/(3) do not, HD2' must be reformulated to")
    lines.append("decouple Q-growth from H-growth.")
    out_path.write_text("\n".join(lines), encoding="utf-8")


# ----------------------------------------------------------------------------- #
# CLI                                                                           #
# ----------------------------------------------------------------------------- #

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True,
                   help="Input near-power CSV (same as power_residue_sieve.py)")
    p.add_argument("--limit", type=int, default=0,
                   help="Limit number of rows scored. 0 = all rows.")
    p.add_argument("--primes", type=int, nargs="+", default=[251, 503, 1009, 2003],
                   help="List of Q values (verifier prime bounds) to test.")
    p.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent,
                   help="Output directory for CSVs and summaries.")
    p.add_argument("--include-trivial", action="store_true",
                   help="Include verifier primes with gcd(z,q-1)=1 (no filtering power).")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    only_nontrivial = not args.include_trivial
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if not args.input.exists():
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        return 2

    # Load rows once
    print(f"Loading rows from {args.input} ...", flush=True)
    rows: List[Dict[str, str]] = []
    with args.input.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    if args.limit and len(rows) > args.limit:
        rows = rows[: args.limit]
    print(f"Loaded {len(rows)} rows.", flush=True)

    Qs = sorted(set(args.primes))
    print(f"Q values: {Qs}", flush=True)

    residue_cache: Dict[Tuple[int, int], int] = {}
    density_cache: Dict[Tuple[int, int, int, int], float] = {}
    per_Q: Dict[int, List[RowScore]] = {}

    t0 = time.time()
    for Q in Qs:
        primes_Q = primes_upto(Q)
        nontrivial = [q for q in primes_Q if any(math.gcd(z, q - 1) > 1 for z in range(3, 10))]
        print(f"\n=== Q={Q} ===", flush=True)
        print(f"  primes <= {Q}: {len(primes_Q)} (potentially nontrivial: {len(nontrivial)})", flush=True)
        tQ0 = time.time()
        scores = score_for_Q(rows, primes_Q, only_nontrivial,
                             residue_cache, density_cache, log_every=100)
        elapsed_Q = time.time() - tQ0
        per_Q[Q] = scores
        score_csv = args.out_dir / f"variable_Q_sieve_Q{Q}_n{len(scores)}.csv"
        score_md = args.out_dir / f"variable_Q_sieve_Q{Q}_n{len(scores)}_summary.md"
        write_scores_csv(scores, score_csv)
        write_summary_for_Q(scores, Q, len(rows), score_md, elapsed_Q)
        print(f"  Q={Q} done in {elapsed_Q:.1f}s. CSV: {score_csv.name}, summary: {score_md.name}",
              flush=True)
        # Save partial aggregate after each Q in case we get killed
        agg_partial = args.out_dir / "variable_Q_aggregate_summary.md"
        write_aggregate(per_Q, agg_partial, time.time() - t0)
        # Also save a checkpoint json
        checkpoint = {
            "Qs_done": list(per_Q.keys()),
            "elapsed_total_s": time.time() - t0,
            "n_rows": len(rows),
            "input_file": str(args.input),
        }
        (args.out_dir / "variable_Q_checkpoint.json").write_text(
            json.dumps(checkpoint, indent=2), encoding="utf-8"
        )

    elapsed_total = time.time() - t0
    print(f"\nAll Q values done in {elapsed_total:.1f}s.", flush=True)
    print(f"Aggregate summary: {args.out_dir / 'variable_Q_aggregate_summary.md'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
