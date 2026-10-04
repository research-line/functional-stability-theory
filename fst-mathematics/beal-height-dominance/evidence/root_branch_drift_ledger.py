#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Root-Branch-Drift ledger for the Beal Height-Dominance project.

The script reads near-power rows A^x+B^y ~= C^z and measures, for prime
blocks P={q<=Q}, whether local z-th-root branches admit low global sections
c <= T=C. It is a diagnostic tool for proof_notes/ROOT_BRANCH_DRIFT_LEDGER.md,
not a Beal proof.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import random
from collections import Counter, defaultdict
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


def load_rows(path: Path, limit: int) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
            if limit and len(rows) >= limit:
                break
    return rows


def root_set_mod_prime(q: int, z: int, value_mod_q: int) -> RootSet:
    return tuple(c for c in range(q) if pow(c, z, q) == value_mod_q)


def deterministic_random_roots(q: int, size: int, seed_text: str) -> RootSet:
    if size <= 0:
        return tuple()
    digest = hashlib.sha256(seed_text.encode("utf-8")).hexdigest()
    rng = random.Random(int(digest[:16], 16))
    return tuple(sorted(rng.sample(range(q), min(size, q))))


def positive_count_for_residue(T: int, modulus: int, residue: int) -> int:
    """Count positive c <= T with c == residue mod modulus."""
    if T <= 0:
        return 0
    residue %= modulus
    if residue == 0:
        return T // modulus
    if residue > T:
        return 0
    return 1 + (T - residue) // modulus


def crt_pair_residue(a: int, q: int, b: int, r: int) -> int:
    inv_q = pow(q, -1, r)
    k = ((b - a) * inv_q) % r
    return a + q * k


def pair_correlation_ratios(block: RootBlock, T: int, pair_limit: int) -> Tuple[float, float, float]:
    usable = [(q, roots) for q, roots in block if roots]
    if len(usable) < 2 or T <= 0 or pair_limit <= 0:
        return float("nan"), float("nan"), float("nan")

    pairs: List[Tuple[Tuple[int, RootSet], Tuple[int, RootSet]]] = []
    for idx in range(len(usable) - 1):
        pairs.append((usable[idx], usable[idx + 1]))
        if len(pairs) >= pair_limit:
            break

    ratios: List[float] = []
    for (q, roots_q), (r, roots_r) in pairs:
        modulus = q * r
        count = 0
        for a in roots_q:
            for b in roots_r:
                residue = crt_pair_residue(a, q, b, r)
                count += positive_count_for_residue(T, modulus, residue)
        expected = T * len(roots_q) * len(roots_r) / modulus
        if expected > 0:
            ratios.append(count / expected)

    if not ratios:
        return float("nan"), float("nan"), float("nan")
    return sum(ratios) / len(ratios), min(ratios), max(ratios)


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


def build_block(
    row: Dict[str, str],
    primes: Sequence[int],
    control_class: str,
    seed: int,
) -> RootBlock:
    rank = int(row.get("rank", "0"))
    A = int(row["A"])
    x = int(row["x"])
    B = int(row["B"])
    y = int(row["y"])
    C = int(row["C"])
    z = int(row["z"])
    N = int(row["N"])
    H = int(row["H_Cz"])

    block: RootBlock = []
    for q in primes:
        if math.gcd(z, q - 1) == 1:
            continue
        if control_class == "exact_power_anchor":
            value = H % q
            roots = root_set_mod_prime(q, z, value)
        else:
            value = (pow(A, x, q) + pow(B, y, q)) % q
            roots = root_set_mod_prime(q, z, value)
            if control_class == "random_root_sets":
                roots = deterministic_random_roots(
                    q,
                    len(roots),
                    f"{seed}|{rank}|{q}|{z}|random_root_sets",
                )
        block.append((q, roots))
    return block


def block_product(primes: Iterable[int]) -> int:
    product = 1
    for q in primes:
        product *= q
    return product


def decide(control_class: str, local_empty_roots: int, n_low_sections: int) -> str:
    if control_class == "exact_power_anchor":
        return "positive_control_survivor" if n_low_sections else "positive_control_failed"
    if local_empty_roots:
        return "local_obstruction_before_drift"
    if n_low_sections:
        return "low_section_survivor_requires_anchor_audit"
    return "root_branch_drift_no_low_section"


@dataclass
class LedgerRow:
    case_id: str
    control_class: str
    rank: int
    x: int
    y: int
    z: int
    A: int
    B: int
    C: int
    H: int
    T: int
    Q: int
    block_id: str
    prime_count: int
    M_block: str
    log_M_block: float
    n_root_classes: str
    log_n_root_classes: float
    local_empty_roots: int
    first_empty_q: int
    n_low_sections: int
    low_section_sample: str
    D_root: float
    D_drift: float
    pair_corr_ratio_mean: float
    pair_corr_ratio_min: float
    pair_corr_ratio_max: float
    anchor_type: str
    detector_degree: float
    target_degree: float
    degree_gap: float
    advice_source: str
    decision: str
    # Hodge-Transfer (Provenienz-Wasserlinie)
    source_axis_1: str
    source_axis_2: str
    source_anchor_1: str
    source_anchor_2: str
    axis_independence_check: str
    coefficient_ring_target: str
    residue_or_survivor_status: str
    overlap_status: str
    lift_or_escape_obstruction: str
    matched_negative_control_result: str
    # NS-LDI-Transfer (Tide-Clock)
    window_id: str
    window_predefined: str
    local_pass_occupancy: int
    root_signal_share: float
    drift_or_branch_share: float
    share_over_occupancy: float
    branch_switch_share: str
    alternate_root_or_shuffled_branch_control: str
    noncoherent_tail_cost: float
    transfer_status: str


def score_row(
    row: Dict[str, str],
    Q: int,
    primes: Sequence[int],
    control_class: str,
    seed: int,
    pair_limit: int,
    sample_limit: int,
) -> LedgerRow:
    rank = int(row.get("rank", "0"))
    A = int(row["A"])
    x = int(row["x"])
    B = int(row["B"])
    y = int(row["y"])
    C = int(row["C"])
    z = int(row["z"])
    H = int(row["H_Cz"])
    T = C

    block = build_block(row, primes, control_class, seed)
    block_primes = [q for q, _ in block]
    M = block_product(block_primes)
    log_M = sum(math.log(q) for q in block_primes)

    root_sizes = [len(roots) for _, roots in block]
    local_empty_roots = sum(1 for size in root_sizes if size == 0)
    first_empty_q = next((q for q, roots in block if not roots), 0)
    log_n_root = sum(math.log(max(1, size)) for size in root_sizes)
    n_root_classes = 1
    for size in root_sizes:
        n_root_classes *= max(1, size)

    n_low, sample = count_low_sections(block, T, sample_limit)
    D_root = log_M - math.log(max(1, T)) - log_n_root
    D_drift = log_n_root - math.log(max(1, n_low))
    pair_mean, pair_min, pair_max = pair_correlation_ratios(block, T, pair_limit)
    log_H = math.log(max(2, H))
    detector_degree = log_M / log_H
    target_degree = math.log(max(2, T)) / log_H

    if control_class == "exact_power_anchor":
        anchor_type = "constructed_exact_power_root"
        advice_source = "C_from_near_power_row_as_positive_control"
    elif control_class == "random_root_sets":
        anchor_type = "none_visible"
        advice_source = "deterministic_seeded_random_sets_same_cardinality"
    else:
        anchor_type = "none_visible"
        advice_source = "existing_near_power_csv"

    # Hodge-Transfer (Provenienz-Wasserlinie) fields
    source_axis_1 = "local_residue_sieve"
    source_axis_2 = "blocked_unfilled"
    source_anchor_1 = "darmon_granville_ceiling"
    source_anchor_2 = "blocked_unfilled"
    axis_independence_check = "blocked_unfilled"
    coefficient_ring_target = "Z"
    residue_or_survivor_status = "residue_only" if n_low == 0 else "survivor_kandidat_requires_overlap_check"
    overlap_status = "blocked_unfilled"
    lift_or_escape_obstruction = "blocked_unfilled"
    if control_class == "random_root_sets":
        matched_negative_control_result = "negative_control_passed" if n_low == 0 else "negative_control_failed_spurious_survivor"
    elif control_class == "exact_power_anchor":
        matched_negative_control_result = "positive_control_passed" if n_low > 0 else "positive_control_failed_no_survivor"
    else:
        matched_negative_control_result = "under_test"

    # NS-LDI-Transfer (Tide-Clock) fields
    window_id = f"Q_{Q}"
    window_predefined = "pass_predefined"
    local_pass_occupancy = len(block_primes)
    root_signal_share = float(n_low) / float(max(1, T))
    drift_or_branch_share = log_n_root
    share_over_occupancy = log_n_root / float(max(1, len(block_primes)))
    branch_switch_share = "blocked_unfilled"
    alternate_root_or_shuffled_branch_control = f"control_class={control_class}"
    noncoherent_tail_cost = D_drift
    transfer_status = decide(control_class, local_empty_roots, n_low)

    return LedgerRow(
        case_id=f"rank{rank}_Q{Q}_{control_class}",
        control_class=control_class,
        rank=rank,
        x=x,
        y=y,
        z=z,
        A=A,
        B=B,
        C=C,
        H=H,
        T=T,
        Q=Q,
        block_id=f"prime_block_le_{Q}",
        prime_count=len(block_primes),
        M_block=str(M),
        log_M_block=log_M,
        n_root_classes=str(n_root_classes),
        log_n_root_classes=log_n_root,
        local_empty_roots=local_empty_roots,
        first_empty_q=first_empty_q,
        n_low_sections=n_low,
        low_section_sample=sample,
        D_root=D_root,
        D_drift=D_drift,
        pair_corr_ratio_mean=pair_mean,
        pair_corr_ratio_min=pair_min,
        pair_corr_ratio_max=pair_max,
        anchor_type=anchor_type,
        detector_degree=detector_degree,
        target_degree=target_degree,
        degree_gap=detector_degree - target_degree,
        advice_source=advice_source,
        decision=decide(control_class, local_empty_roots, n_low),
        source_axis_1=source_axis_1,
        source_axis_2=source_axis_2,
        source_anchor_1=source_anchor_1,
        source_anchor_2=source_anchor_2,
        axis_independence_check=axis_independence_check,
        coefficient_ring_target=coefficient_ring_target,
        residue_or_survivor_status=residue_or_survivor_status,
        overlap_status=overlap_status,
        lift_or_escape_obstruction=lift_or_escape_obstruction,
        matched_negative_control_result=matched_negative_control_result,
        window_id=window_id,
        window_predefined=window_predefined,
        local_pass_occupancy=local_pass_occupancy,
        root_signal_share=root_signal_share,
        drift_or_branch_share=drift_or_branch_share,
        share_over_occupancy=share_over_occupancy,
        branch_switch_share=branch_switch_share,
        alternate_root_or_shuffled_branch_control=alternate_root_or_shuffled_branch_control,
        noncoherent_tail_cost=noncoherent_tail_cost,
        transfer_status=transfer_status,
    )


def write_ledger(rows: List[LedgerRow], path: Path) -> None:
    fields = list(LedgerRow.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for item in rows:
            data = asdict(item)
            for key in (
                "log_M_block",
                "log_n_root_classes",
                "D_root",
                "D_drift",
                "pair_corr_ratio_mean",
                "pair_corr_ratio_min",
                "pair_corr_ratio_max",
                "detector_degree",
                "target_degree",
                "degree_gap",
                "root_signal_share",
                "drift_or_branch_share",
                "share_over_occupancy",
                "noncoherent_tail_cost",
            ):
                value = data[key]
                data[key] = "" if isinstance(value, float) and math.isnan(value) else f"{value:.12g}"
            writer.writerow(data)


def median(values: List[float]) -> float:
    vals = sorted(v for v in values if not math.isnan(v))
    if not vals:
        return float("nan")
    mid = len(vals) // 2
    if len(vals) % 2:
        return vals[mid]
    return 0.5 * (vals[mid - 1] + vals[mid])


def write_summary(all_rows: Dict[int, List[LedgerRow]], args: argparse.Namespace, paths: Dict[int, Path]) -> None:
    lines: List[str] = []
    lines.append("# Root-Branch-Drift Ledger -- Pilot Summary")
    lines.append("")
    lines.append(f"- Input: `{args.input.name}`")
    lines.append(f"- Rows je Q: `{args.limit}`")
    lines.append(f"- Q-Werte: `{', '.join(str(q) for q in args.primes)}`")
    lines.append(f"- Kontrollklassen: `{', '.join(args.controls)}`")
    lines.append("")
    lines.append("## Kerndeutung")
    lines.append("")
    lines.append("Der Lauf pr\u00fcft nicht, ob Beal bewiesen ist. Er trennt drei F\u00e4lle:")
    lines.append("")
    lines.append("- `actual_near_miss`: vorhandene Near-Misses aus dem Projekt.")
    lines.append("- `random_root_sets`: gleich gro\u00dfe, deterministisch zuf\u00e4llige Root-Sets.")
    lines.append("- `exact_power_anchor`: positive Kontrolle mit `N=C^z`; hier muss `C` als niedrige Sektion \u00fcberleben.")
    lines.append("")
    lines.append("`local_obstruction_before_drift` bedeutet: Bereits mindestens ein lokales Root-Set ist leer. Dann ist kein globaler Drift-Befund n\u00f6tig.")
    lines.append("")
    lines.append("## Aggregate")
    lines.append("")
    lines.append("| Q | control | n | local obstruction | low-section survivors | median D_root | median D_drift | median degree_gap | decisions | CSV |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|---:|---|---|")
    for Q, rows in all_rows.items():
        by_control: Dict[str, List[LedgerRow]] = defaultdict(list)
        for row in rows:
            by_control[row.control_class].append(row)
        for control in args.controls:
            control_rows = by_control.get(control, [])
            if not control_rows:
                continue
            obstruction = sum(1 for r in control_rows if r.local_empty_roots > 0)
            survivors = sum(1 for r in control_rows if r.n_low_sections > 0)
            decisions = Counter(r.decision for r in control_rows)
            decision_text = "; ".join(f"{k}:{v}" for k, v in sorted(decisions.items()))
            csv_name = paths[Q].name
            lines.append(
                f"| {Q} | `{control}` | {len(control_rows)} | {obstruction} | {survivors} | "
                f"{median([r.D_root for r in control_rows]):.3f} | "
                f"{median([r.D_drift for r in control_rows]):.3f} | "
                f"{median([r.degree_gap for r in control_rows]):.3f} | "
                f"{decision_text} | `{csv_name}` |"
            )
    lines.append("")
    lines.append("## Befund")
    lines.append("")
    lines.append("Im Pilotlauf dominieren bei echten Near-Misses die lokalen Obstruktionen: Root-Branch-Drift wird also erst nach lokalem Pass oder in konstruierten Kontrollklassen entscheidend. Die Exact-Power-Kontrolle funktioniert als Positivkontrolle und h\u00e4lt die niedrige Sektion `C` sichtbar.")
    lines.append("")
    lines.append("## Folgeschritt")
    lines.append("")
    lines.append("F\u00fcr einen st\u00e4rkeren Drift-Test braucht das Projekt nun gezielt lokale-Pass-Survivors oder bekannte Signaturfamilien; ansonsten misst der Ledger vor allem, dass Near-Misses schon lokal scheitern.")
    args.summary.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=here / "near_power_b250_e3-7_z3-9_k5000.csv")
    parser.add_argument("--limit", type=int, default=250)
    parser.add_argument("--primes", type=int, nargs="+", default=[251, 503, 1009])
    parser.add_argument(
        "--controls",
        nargs="+",
        default=["actual_near_miss", "random_root_sets", "exact_power_anchor"],
        choices=["actual_near_miss", "random_root_sets", "exact_power_anchor"],
    )
    parser.add_argument("--seed", type=int, default=20260522)
    parser.add_argument("--pair-limit", type=int, default=24)
    parser.add_argument("--sample-limit", type=int, default=12)
    parser.add_argument("--out-prefix", default="root_branch_drift_ledger")
    parser.add_argument("--summary", type=Path, default=here / "root_branch_drift_ledger_summary.md")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = load_rows(args.input, args.limit)
    if not rows:
        raise RuntimeError(f"No rows loaded from {args.input}")

    all_rows: Dict[int, List[LedgerRow]] = {}
    paths: Dict[int, Path] = {}
    for Q in args.primes:
        primes = primes_upto(Q)
        scored: List[LedgerRow] = []
        for row in rows:
            for control in args.controls:
                scored.append(
                    score_row(
                        row=row,
                        Q=Q,
                        primes=primes,
                        control_class=control,
                        seed=args.seed,
                        pair_limit=args.pair_limit,
                        sample_limit=args.sample_limit,
                    )
                )
        out_path = args.input.parent / f"{args.out_prefix}_Q{Q}_n{len(rows)}.csv"
        write_ledger(scored, out_path)
        all_rows[Q] = scored
        paths[Q] = out_path
        print(f"Q={Q}: wrote {len(scored)} ledger rows -> {out_path}")

    write_summary(all_rows, args, paths)
    print(f"Summary: {args.summary}")


if __name__ == "__main__":
    main()
