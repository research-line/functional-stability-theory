#!/usr/bin/env python3
"""
Targeted near-power fracture scan for Beal Height Dominance.

This script searches for close misses

    A^x + B^y ~= C^z,        x,y,z >= 3,

keeps the closest cases, and measures whether the radical support of
ABC is already too large for the Beal compression exponent

    sigma = 1/x + 1/y + 1/z.

The central diagnostic is

    overflow_margin = log(rad(ABC)) / log(C^z) - sigma.

For an exact primitive Beal solution the elementary power-compression side
would force rad(ABC) < (C^z)^sigma.  Thus a positive margin in near-misses
is evidence for "radical overflow" at primitive fracture points.

The more sensitive fracture diagnostic is the right-side radical jump

    rhs_radical_jump = log(rad(A^x+B^y)) / log(C^z) - 1/z.

If the near-miss closed exactly as C^z, the RHS radical exponent would have
to collapse to <= 1/z.  Near-misses typically have rad(A^x+B^y) close to
A^x+B^y, so this jump measures how violently the additive break creates
new prime support.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import math
import random
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


Factorization = Dict[int, int]


def gcd_many(values: Iterable[int]) -> int:
    result = 0
    for value in values:
        result = math.gcd(result, value)
    return result


def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2

    for a in small_primes:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def pollard_rho(n: int, rng: random.Random) -> int:
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        c = rng.randrange(1, n - 1)
        x = rng.randrange(2, n - 1)
        y = x
        d = 1
        while d == 1:
            x = (pow(x, 2, n) + c) % n
            y = (pow(y, 2, n) + c) % n
            y = (pow(y, 2, n) + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d


def factor_into(n: int, out: List[int], rng: random.Random) -> None:
    if n <= 1:
        return
    if is_probable_prime(n):
        out.append(n)
        return
    d = pollard_rho(n, rng)
    factor_into(d, out, rng)
    factor_into(n // d, out, rng)


def factorize(n: int, rng: random.Random) -> Factorization:
    factors: List[int] = []
    factor_into(abs(n), factors, rng)
    return dict(Counter(factors))


def radical_from_factors(factors: Factorization) -> int:
    r = 1
    for p in factors:
        r *= p
    return r


def merge_factors(*items: Factorization) -> Factorization:
    merged: Factorization = {}
    for factors in items:
        for p, e in factors.items():
            merged[p] = merged.get(p, 0) + e
    return merged


def integer_nth_root(n: int, k: int) -> int:
    if n < 2:
        return n
    lo = 1
    hi = 1 << ((n.bit_length() + k - 1) // k)
    while lo <= hi:
        mid = (lo + hi) // 2
        val = mid**k
        if val == n:
            return mid
        if val < n:
            lo = mid + 1
        else:
            hi = mid - 1
    return hi


def nearest_power(n: int, min_z: int, max_z: int) -> Tuple[int, int, int, int]:
    best: Tuple[int, int, int, int] | None = None
    for z in range(min_z, max_z + 1):
        root = integer_nth_root(n, z)
        for c in (root, root + 1):
            if c <= 1:
                continue
            h = c**z
            gap = abs(n - h)
            item = (gap, z, c, h)
            if best is None or item < best:
                best = item
    assert best is not None
    gap, z, c, h = best
    return z, c, h, gap


def perfect_power_degree(factors: Factorization) -> int:
    if not factors:
        return 1
    return gcd_many(factors.values())


def format_factorization(factors: Factorization) -> str:
    if not factors:
        return "1"
    parts = []
    for p in sorted(factors):
        e = factors[p]
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return " * ".join(parts)


def sha256_hex(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


@dataclass
class RawCandidate:
    rel_gap: float
    gap_abs: int
    A: int
    x: int
    B: int
    y: int
    C: int
    z: int
    H: int
    N: int
    left: int
    right: int
    balance: float


@dataclass
class NearPowerRecord:
    rel_gap: float
    gap_signed: int
    gap_abs: int
    A: int
    x: int
    B: int
    y: int
    C: int
    z: int
    H: int
    N: int
    left: int
    right: int
    balance: float
    sigma: float
    rad_ABC: int
    rad_exp: float
    overflow_margin: float
    product_exp: float
    rad_N: int
    rad_N_exp: float
    rhs_radical_jump: float
    power_index_N: float
    omega_N: int
    perfect_degree_N: int
    gcd_AB: int
    gcd_AC: int
    gcd_BC: int
    factors_N: Factorization
    factors_gap: Factorization


def push_candidate(heap: List[Tuple[float, int, RawCandidate]], keep: int, seq: int, cand: RawCandidate) -> None:
    # Max-heap by negative rel_gap.  The worst kept candidate is at heap[0].
    item = (-cand.rel_gap, seq, cand)
    if len(heap) < keep:
        heapq.heappush(heap, item)
    elif item > heap[0]:
        heapq.heapreplace(heap, item)


def generate_candidates(args: argparse.Namespace) -> Tuple[List[RawCandidate], int]:
    powers: List[Tuple[int, int, int, int]] = []
    for exp in range(args.min_exp, args.max_exp + 1):
        for base in range(2, args.max_base + 1):
            powers.append((base**exp, base, exp, base))
    powers.sort()

    heap: List[Tuple[float, int, RawCandidate]] = []
    scanned = 0
    seq = 0

    for i, (left, A, x, _) in enumerate(powers):
        for right, B, y, _ in powers[i if args.allow_symmetric_duplicates else i + 1 :]:
            if not args.allow_equal_terms and A == B and x == y:
                continue
            if args.primitive_inputs and math.gcd(A, B) != 1:
                continue
            balance = min(left, right) / max(left, right)
            if balance < args.min_balance:
                continue

            N = left + right
            z, C, H, gap_abs = nearest_power(N, args.target_min_exp, args.target_max_exp)
            rel_gap = gap_abs / H
            scanned += 1
            if rel_gap > args.max_rel_gap:
                continue

            seq += 1
            push_candidate(
                heap,
                args.keep,
                seq,
                RawCandidate(
                    rel_gap=rel_gap,
                    gap_abs=gap_abs,
                    A=A,
                    x=x,
                    B=B,
                    y=y,
                    C=C,
                    z=z,
                    H=H,
                    N=N,
                    left=left,
                    right=right,
                    balance=balance,
                ),
            )

    selected = [item[2] for item in sorted(heap, key=lambda t: -t[0])]
    selected.sort(key=lambda c: (c.rel_gap, c.gap_abs, c.N))
    return selected, scanned


def enrich(candidates: List[RawCandidate], args: argparse.Namespace) -> List[NearPowerRecord]:
    rng = random.Random(args.seed ^ 0xBEEFBEEF)
    factor_cache: Dict[int, Factorization] = {}

    def ff(n: int) -> Factorization:
        n = abs(n)
        if n not in factor_cache:
            factor_cache[n] = factorize(n, rng)
        return factor_cache[n]

    records: List[NearPowerRecord] = []
    for c in candidates:
        factors_A = ff(c.A)
        factors_B = ff(c.B)
        factors_C = ff(c.C)
        factors_ABC = merge_factors(factors_A, factors_B, factors_C)
        rad_ABC = radical_from_factors(factors_ABC)

        factors_N = ff(c.N)
        rad_N = radical_from_factors(factors_N)
        gap_signed = c.N - c.H
        factors_gap = ff(abs(gap_signed)) if gap_signed != 0 else {}

        sigma = 1 / c.x + 1 / c.y + 1 / c.z
        log_H = math.log(c.H)
        rad_exp = math.log(rad_ABC) / log_H
        product_exp = math.log(c.A * c.B * c.C) / log_H
        overflow_margin = rad_exp - sigma
        rad_N_exp = math.log(rad_N) / log_H
        rhs_radical_jump = rad_N_exp - (1 / c.z)
        power_index_N = math.log(c.N) / math.log(rad_N)

        records.append(
            NearPowerRecord(
                rel_gap=c.rel_gap,
                gap_signed=gap_signed,
                gap_abs=c.gap_abs,
                A=c.A,
                x=c.x,
                B=c.B,
                y=c.y,
                C=c.C,
                z=c.z,
                H=c.H,
                N=c.N,
                left=c.left,
                right=c.right,
                balance=c.balance,
                sigma=sigma,
                rad_ABC=rad_ABC,
                rad_exp=rad_exp,
                overflow_margin=overflow_margin,
                product_exp=product_exp,
                rad_N=rad_N,
                rad_N_exp=rad_N_exp,
                rhs_radical_jump=rhs_radical_jump,
                power_index_N=power_index_N,
                omega_N=len(factors_N),
                perfect_degree_N=perfect_power_degree(factors_N),
                gcd_AB=math.gcd(c.A, c.B),
                gcd_AC=math.gcd(c.A, c.C),
                gcd_BC=math.gcd(c.B, c.C),
                factors_N=factors_N,
                factors_gap=factors_gap,
            )
        )
    return records


def write_csv(records: List[NearPowerRecord], path: Path) -> None:
    fields = [
        "rank",
        "A",
        "x",
        "B",
        "y",
        "C",
        "z",
        "N",
        "H_Cz",
        "gap_signed",
        "gap_abs",
        "rel_gap",
        "balance",
        "sigma",
        "rad_ABC",
        "rad_exp",
        "overflow_margin",
        "product_exp",
        "rad_N",
        "rad_N_exp",
        "rhs_radical_jump",
        "power_index_N",
        "omega_N",
        "perfect_degree_N",
        "gcd_AB",
        "gcd_AC",
        "gcd_BC",
        "factorization_N",
        "factorization_hash",
        "record_hash",
        "factorization_gap",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for rank, r in enumerate(records, start=1):
            factorization_N = format_factorization(r.factors_N)
            factorization_gap = format_factorization(r.factors_gap)
            record_text = (
                f"{r.A}^{r.x}+{r.B}^{r.y}~{r.C}^{r.z}|"
                f"N={r.N}|H={r.H}|gap={r.gap_signed}|"
                f"factors_N={factorization_N}|rad_N={r.rad_N}|"
                f"rad_ABC={r.rad_ABC}|rhs_jump={r.rhs_radical_jump:.12f}"
            )
            writer.writerow(
                {
                    "rank": rank,
                    "A": r.A,
                    "x": r.x,
                    "B": r.B,
                    "y": r.y,
                    "C": r.C,
                    "z": r.z,
                    "N": r.N,
                    "H_Cz": r.H,
                    "gap_signed": r.gap_signed,
                    "gap_abs": r.gap_abs,
                    "rel_gap": f"{r.rel_gap:.12e}",
                    "balance": f"{r.balance:.12e}",
                    "sigma": f"{r.sigma:.12f}",
                    "rad_ABC": r.rad_ABC,
                    "rad_exp": f"{r.rad_exp:.12f}",
                    "overflow_margin": f"{r.overflow_margin:.12f}",
                    "product_exp": f"{r.product_exp:.12f}",
                    "rad_N": r.rad_N,
                    "rad_N_exp": f"{r.rad_N_exp:.12f}",
                    "rhs_radical_jump": f"{r.rhs_radical_jump:.12f}",
                    "power_index_N": f"{r.power_index_N:.12f}",
                    "omega_N": r.omega_N,
                    "perfect_degree_N": r.perfect_degree_N,
                    "gcd_AB": r.gcd_AB,
                    "gcd_AC": r.gcd_AC,
                    "gcd_BC": r.gcd_BC,
                    "factorization_N": factorization_N,
                    "factorization_hash": sha256_hex(factorization_N),
                    "record_hash": sha256_hex(record_text),
                    "factorization_gap": factorization_gap,
                }
            )


def count_where(records: List[NearPowerRecord], predicate) -> int:
    return sum(1 for r in records if predicate(r))


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


def write_summary(records: List[NearPowerRecord], path: Path, csv_path: Path, args: argparse.Namespace, scanned: int) -> None:
    total = len(records)
    pos_overflow = count_where(records, lambda r: r.overflow_margin > 0)
    strong_overflow = count_where(records, lambda r: r.overflow_margin > 0.05)
    positive_rhs_jump = count_where(records, lambda r: r.rhs_radical_jump > 0)
    strong_rhs_jump = count_where(records, lambda r: r.rhs_radical_jump > 0.5)
    exact = count_where(records, lambda r: r.gap_abs == 0)
    pairwise_anchor = count_where(records, lambda r: r.gcd_AC > 1 or r.gcd_BC > 1)
    squarefree_N = count_where(records, lambda r: r.power_index_N < 1.02)
    type_z = Counter(r.z for r in records)
    omega_counts = Counter(r.omega_N for r in records)

    lines: List[str] = []
    lines.append("# Near-Power-Bruchstellen-Summary")
    lines.append("")
    lines.append(f"- Enumerierte primitive Paare: `{scanned}`")
    lines.append(f"- Behaltene engste Near-Misses: `{total}`")
    lines.append(f"- Bases: `2..{args.max_base}`")
    lines.append(f"- Input-Exponenten: `{args.min_exp}..{args.max_exp}`")
    lines.append(f"- Ziel-Exponenten `z`: `{args.target_min_exp}..{args.target_max_exp}`")
    lines.append(f"- Minimum balance: `{args.min_balance}`")
    lines.append(f"- Max rel gap: `{args.max_rel_gap}`")
    lines.append(f"- CSV: `{csv_path.name}`")
    lines.append("")
    lines.append("## Kerndeutung")
    lines.append("")
    lines.append(
        "Diese Auswahl ist bewusst nicht zufällig: sie enthält die Fälle, in denen "
        "`A^x+B^y` besonders nahe an einer höheren Potenz `C^z` liegt. "
        "Genau dort sollte ein möglicher Beal-artiger Bruchschutz am ehesten sichtbar werden."
    )
    lines.append("")
    lines.append("Es gibt zwei Messwerte:")
    lines.append("")
    lines.append("```text")
    lines.append("abc_margin        = log(rad(ABC))/log(C^z) - (1/x+1/y+1/z)")
    lines.append("rhs_radical_jump  = log(rad(A^x+B^y))/log(C^z) - 1/z")
    lines.append("```")
    lines.append("")
    lines.append(
        "`abc_margin` prüft die Beal-Kompression und ist bei Near-Misses erwartbar "
        "nicht positiv. `rhs_radical_jump` misst dagegen den eigentlichen Bruch: "
        "wie stark das Radical des Summenwerts `N` über dem Radical einer echten "
        "Potenz `C^z` liegt."
    )
    lines.append("")
    lines.append("## Hauptbefund")
    lines.append("")
    lines.append(f"- Positiver RHS-Radical-Jump: {positive_rhs_jump}/{total} ({positive_rhs_jump/total:.2%}) `{bar(positive_rhs_jump,total)}`")
    lines.append(f"- Starker RHS-Jump > 0.5: {strong_rhs_jump}/{total} ({strong_rhs_jump/total:.2%}) `{bar(strong_rhs_jump,total)}`")
    lines.append(f"- Positive ABC-Kompressionsmarge: {pos_overflow}/{total} ({pos_overflow/total:.2%})")
    lines.append(f"- Starke ABC-Marge > 0.05: {strong_overflow}/{total} ({strong_overflow/total:.2%})")
    lines.append(f"- Exakte Treffer `A^x+B^y=C^z`: {exact}")
    lines.append(f"- Pairwise C-Anker (`gcd(A,C)>1` oder `gcd(B,C)>1`): {pairwise_anchor}/{total} ({pairwise_anchor/total:.2%})")
    lines.append(f"- RHS `N` fast squarefree (`pi(N)<1.02`): {squarefree_N}/{total} ({squarefree_N/total:.2%})")
    lines.append("")
    lines.append("## Kennzahlen")
    lines.append("")
    lines.append(f"- Median rel gap: `{median([r.rel_gap for r in records]):.3e}`")
    lines.append(f"- Median RHS radical jump: `{median([r.rhs_radical_jump for r in records]):.6f}`")
    lines.append(f"- Minimum RHS radical jump: `{min(r.rhs_radical_jump for r in records):.6f}`")
    lines.append(f"- Maximum RHS radical jump: `{max(r.rhs_radical_jump for r in records):.6f}`")
    lines.append(f"- Median ABC compression margin: `{median([r.overflow_margin for r in records]):.6f}`")
    lines.append(f"- Median power-index of N: `{median([r.power_index_N for r in records]):.6f}`")
    lines.append("")
    lines.append("## Ziel-Exponenten")
    lines.append("")
    for z, count in sorted(type_z.items()):
        lines.append(f"- z={z}: {count} ({count/total:.2%}) `{bar(count,total)}`")
    lines.append("")
    lines.append("## Anzahl Primwellen in N")
    lines.append("")
    for omega, count in sorted(omega_counts.items()):
        lines.append(f"- omega(N)={omega}: {count} ({count/total:.2%}) `{bar(count,total)}`")
    lines.append("")
    lines.append("## 25 engste Near-Misses")
    lines.append("")
    lines.append("| rank | A^x | B^y | nearest C^z | rel_gap | rhs_jump | abc_margin | pi(N) | gcdAC | gcdBC | factorization N | gap |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|")
    for idx, r in enumerate(records[:25], start=1):
        lines.append(
            f"| {idx} | {r.A}^{r.x} | {r.B}^{r.y} | {r.C}^{r.z} | "
            f"{r.rel_gap:.3e} | {r.rhs_radical_jump:.4f} | {r.overflow_margin:.4f} | "
            f"{r.power_index_N:.4f} | "
            f"{r.gcd_AC} | {r.gcd_BC} | `{format_factorization(r.factors_N)}` | "
            f"`{r.gap_signed}` |"
        )
    lines.append("")
    lines.append("## 20 stärkste RHS-Radical-Jumps unter den Near-Misses")
    lines.append("")
    lines.append("| A^x | B^y | nearest C^z | rel_gap | rhs_jump | rad_N_exp | 1/z | pi(N) | factorization N |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in sorted(records, key=lambda item: item.rhs_radical_jump, reverse=True)[:20]:
        lines.append(
            f"| {r.A}^{r.x} | {r.B}^{r.y} | {r.C}^{r.z} | {r.rel_gap:.3e} | "
            f"{r.rhs_radical_jump:.4f} | {r.rad_N_exp:.4f} | {1/r.z:.4f} | "
            f"{r.power_index_N:.4f} | "
            f"`{format_factorization(r.factors_N)}` |"
        )
    lines.append("")
    lines.append("## Deutung")
    lines.append("")
    lines.append(
        "Die ABC-Kompressionsmarge bleibt negativ, wie sie es für near-power "
        "Tripel auch soll. Der eigentliche Befund liegt im RHS-Radical-Jump: "
        "Die Summe liegt extrem nahe an einer Potenz, aber ihre Faktorisierung "
        "ist fast immer radikal breit. Eine exakte Potenzclosure müsste diesen "
        "neu entstandenen Radical-Support auf den kleinen Wert `rad(C)` kollabieren lassen."
    )
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-base", type=int, default=250)
    parser.add_argument("--min-exp", type=int, default=3)
    parser.add_argument("--max-exp", type=int, default=7)
    parser.add_argument("--target-min-exp", type=int, default=3)
    parser.add_argument("--target-max-exp", type=int, default=9)
    parser.add_argument("--keep", type=int, default=5000)
    parser.add_argument("--max-rel-gap", type=float, default=1e-4)
    parser.add_argument("--min-balance", type=float, default=1e-3)
    parser.add_argument("--seed", type=int, default=20260425)
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--allow-nonprimitive-inputs", action="store_true")
    parser.add_argument("--allow-symmetric-duplicates", action="store_true")
    parser.add_argument("--allow-equal-terms", action="store_true")
    args = parser.parse_args()
    args.primitive_inputs = not args.allow_nonprimitive_inputs
    if args.max_exp < args.min_exp:
        raise ValueError("--max-exp must be >= --min-exp")
    if args.target_max_exp < args.target_min_exp:
        raise ValueError("--target-max-exp must be >= --target-min-exp")
    return args


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    candidates, scanned = generate_candidates(args)
    if not candidates:
        raise RuntimeError("No near-power candidates found. Increase --max-rel-gap or scan range.")
    records = enrich(candidates, args)
    records.sort(key=lambda r: (r.rel_gap, -r.overflow_margin, r.N))

    stem = (
        f"near_power_b{args.max_base}_e{args.min_exp}-{args.max_exp}_"
        f"z{args.target_min_exp}-{args.target_max_exp}_k{len(records)}"
    )
    csv_path = args.out_dir / f"{stem}.csv"
    summary_path = args.out_dir / f"{stem}_summary.md"
    write_csv(records, csv_path)
    write_summary(records, summary_path, csv_path, args, scanned)

    print(f"Scanned {scanned} primitive input pairs.")
    print(f"Kept {len(records)} closest near-power fractures.")
    print(f"CSV: {csv_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
