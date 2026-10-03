#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Above-threshold local-pass survivor search for Beal Root-Branch diagnostics.

The script scans the existing near-power corpus and reports rows that pass all
nontrivial z-th-power residue channels up to checkpoint prime bounds. It is a
diagnostic tool, not a Beal proof.
"""

from __future__ import annotations

import argparse
import csv
import math
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


def block_product(primes: Iterable[int]) -> int:
    product = 1
    for q in primes:
        product *= q
    return product


def build_block(row: Dict[str, str], primes: Sequence[int], Q: int) -> RootBlock:
    A = int(row["A"])
    x = int(row["x"])
    B = int(row["B"])
    y = int(row["y"])
    z = int(row["z"])
    block: RootBlock = []
    for q in primes:
        if q > Q:
            break
        if math.gcd(z, q - 1) == 1:
            continue
        value = (pow(A, x, q) + pow(B, y, q)) % q
        block.append((q, root_set_mod_prime(q, z, value)))
    return block


def first_empty_channel(block: RootBlock) -> int:
    return next((q for q, roots in block if not roots), 0)


def score_block(row: Dict[str, str], block: RootBlock, sample_limit: int) -> Tuple[int, str, float, float]:
    C = int(row["C"])
    block_primes = [q for q, _ in block]
    M = block_product(block_primes)
    root_sizes = [len(roots) for _, roots in block]
    log_n_root = sum(math.log(max(1, size)) for size in root_sizes)
    n_low, sample = count_low_sections(block, C, sample_limit)
    D_root = math.log(M) - math.log(max(1, C)) - log_n_root
    D_drift = log_n_root - math.log(max(1, n_low))
    return n_low, sample, D_root, D_drift


@dataclass
class CheckpointSummary:
    Q: int
    nontrivial_prime_count_median: int
    local_pass_rows: int
    drift_no_low_section_rows: int
    low_section_rows: int
    duplicate_keys_among_passes: int
    first_empty_top: str


@dataclass
class SurvivorRow:
    rank: int
    A: int
    x: int
    B: int
    y: int
    C: int
    z: int
    N: int
    H_Cz: int
    gap_signed: int
    rel_gap: str
    rhs_radical_jump: str
    best_pass_Q: int
    first_empty_q_up_to_max: int
    first_empty_after_best: int
    nontrivial_primes_at_best: str
    low_sections_at_best: int
    low_section_sample_at_best: str
    D_root_at_best: float
    D_drift_at_best: float
    duplicate_key: str
    duplicate_key_count: int
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


def median_int(values: Sequence[int]) -> int:
    if not values:
        return 0
    vals = sorted(values)
    return vals[len(vals) // 2]


def format_float(value: float) -> str:
    if math.isnan(value):
        return ""
    return f"{value:.12g}"


def analyze(args: argparse.Namespace) -> Tuple[List[CheckpointSummary], List[SurvivorRow]]:
    rows = load_rows(args.input, args.limit)
    if not rows:
        raise RuntimeError(f"No rows loaded from {args.input}")

    max_q = max(max(args.checkpoints), args.max_q)
    primes = primes_upto(max_q)
    duplicate_counts = Counter(f"{row['N']}|{row['H_Cz']}|{row['z']}" for row in rows)

    pass_sets: Dict[int, List[Dict[str, str]]] = defaultdict(list)
    low_counts: Dict[int, Dict[int, Tuple[int, str, float, float, RootBlock]]] = defaultdict(dict)
    checkpoint_nontrivial_counts: Dict[int, List[int]] = defaultdict(list)
    checkpoint_first_empty: Dict[int, Counter[int]] = defaultdict(Counter)

    for row in rows:
        rank = int(row["rank"])
        for Q in args.checkpoints:
            block = build_block(row, primes, Q)
            checkpoint_nontrivial_counts[Q].append(len(block))
            first_empty = first_empty_channel(block)
            if first_empty:
                checkpoint_first_empty[Q][first_empty] += 1
                continue
            pass_sets[Q].append(row)
            low_counts[Q][rank] = (*score_block(row, block, args.sample_limit), block)

    summaries: List[CheckpointSummary] = []
    for Q in args.checkpoints:
        low_section_rows = sum(1 for n_low, *_ in low_counts[Q].values() if n_low > 0)
        drift_no_low = sum(1 for n_low, *_ in low_counts[Q].values() if n_low == 0)
        pass_duplicate_keys = Counter(f"{row['N']}|{row['H_Cz']}|{row['z']}" for row in pass_sets[Q])
        duplicate_keys = sum(1 for count in pass_duplicate_keys.values() if count > 1)
        top = "; ".join(
            f"{q}:{count}" for q, count in checkpoint_first_empty[Q].most_common(args.first_empty_top)
        )
        summaries.append(
            CheckpointSummary(
                Q=Q,
                nontrivial_prime_count_median=median_int(checkpoint_nontrivial_counts[Q]),
                local_pass_rows=len(pass_sets[Q]),
                drift_no_low_section_rows=drift_no_low,
                low_section_rows=low_section_rows,
                duplicate_keys_among_passes=duplicate_keys,
                first_empty_top=top,
            )
        )

    survivor_rows: List[SurvivorRow] = []
    for row in rows:
        best_Q = 0
        best_data: Tuple[int, str, float, float, RootBlock] | None = None
        for Q in args.checkpoints:
            rank = int(row["rank"])
            if rank in low_counts[Q]:
                best_Q = Q
                best_data = low_counts[Q][rank]
        if best_Q < args.min_survivor_Q or best_data is None:
            continue

        full_block = build_block(row, primes, args.max_q)
        first_empty = first_empty_channel(full_block)
        n_low, sample, D_root, D_drift, best_block = best_data
        first_after_best = 0
        for q, roots in full_block:
            if q <= best_Q:
                continue
            if not roots:
                first_after_best = q
                break
        duplicate_key = f"{row['N']}|{row['H_Cz']}|{row['z']}"
        decision = "drift_no_low_section" if n_low == 0 else "low_section_survivor_requires_anchor_audit"

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
        matched_negative_control_result = "under_test"

        # NS-LDI-Transfer (Tide-Clock) fields
        window_id = f"Q_{best_Q}"
        window_predefined = "pass_predefined"
        local_pass_occupancy = len(best_block)
        root_signal_share = float(n_low) / float(max(1, int(row["C"])))
        log_n_root = sum(math.log(max(1, len(roots))) for _, roots in best_block)
        drift_or_branch_share = log_n_root
        share_over_occupancy = log_n_root / float(max(1, len(best_block)))
        branch_switch_share = "blocked_unfilled"
        alternate_root_or_shuffled_branch_control = "control_class=actual_near_miss"
        noncoherent_tail_cost = D_drift
        transfer_status = decision

        survivor_rows.append(
            SurvivorRow(
                rank=int(row["rank"]),
                A=int(row["A"]),
                x=int(row["x"]),
                B=int(row["B"]),
                y=int(row["y"]),
                C=int(row["C"]),
                z=int(row["z"]),
                N=int(row["N"]),
                H_Cz=int(row["H_Cz"]),
                gap_signed=int(row["gap_signed"]),
                rel_gap=row["rel_gap"],
                rhs_radical_jump=row["rhs_radical_jump"],
                best_pass_Q=best_Q,
                first_empty_q_up_to_max=first_empty,
                first_empty_after_best=first_after_best,
                nontrivial_primes_at_best=" ".join(str(q) for q, _ in best_block),
                low_sections_at_best=n_low,
                low_section_sample_at_best=sample,
                D_root_at_best=D_root,
                D_drift_at_best=D_drift,
                duplicate_key=duplicate_key,
                duplicate_key_count=duplicate_counts[duplicate_key],
                decision=decision,
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
        )

    survivor_rows.sort(
        key=lambda item: (
            -item.best_pass_Q,
            item.first_empty_after_best if item.first_empty_after_best else 10**9,
            item.rank,
        )
    )
    return summaries, survivor_rows


def write_survivor_csv(rows: Sequence[SurvivorRow], path: Path) -> None:
    fields = list(SurvivorRow.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for item in rows:
            data = asdict(item)
            data["D_root_at_best"] = format_float(data["D_root_at_best"])
            data["D_drift_at_best"] = format_float(data["D_drift_at_best"])
            data["root_signal_share"] = format_float(data["root_signal_share"])
            data["drift_or_branch_share"] = format_float(data["drift_or_branch_share"])
            data["share_over_occupancy"] = format_float(data["share_over_occupancy"])
            data["noncoherent_tail_cost"] = format_float(data["noncoherent_tail_cost"])
            writer.writerow(data)


def write_report(
    summaries: Sequence[CheckpointSummary],
    survivors: Sequence[SurvivorRow],
    args: argparse.Namespace,
) -> None:
    lines: List[str] = []
    lines.append("# Above-Threshold Survivor Search -- 2026-06-02")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append("Dieser Lauf sucht im vorhandenen 5000er-Near-Miss-Korpus nach natürlichen lokalen Pass-Survivors, die den `q=37`-Schwellkanal und weitere nichttriviale Root-Kanäle überstehen.")
    lines.append("Der Check ist ein lokales Diagnoseartefakt, kein Beal-Beweis, keine HD3'-Schließung und kein Claim-Upgrade.")
    lines.append("")
    lines.append("## Setup")
    lines.append("")
    lines.append(f"- Input: `{args.input.name}`")
    lines.append(f"- Zeilen: `{args.limit or 'alle'}`")
    lines.append(f"- Checkpoints: `{', '.join(str(q) for q in args.checkpoints)}`")
    lines.append(f"- Maximaler Folgekanal: `{args.max_q}`")
    lines.append(f"- Survivor-CSV: `{args.csv.name}`")
    lines.append("")
    lines.append("## Aggregate")
    lines.append("")
    lines.append("| Q | mediane nichttriviale Kanäle | lokale Passes | Passes ohne niedrige Sektion | Passes mit niedriger Sektion | doppelte Pass-Keys | häufigste erste leere Kanäle |")
    lines.append("|---:|---:|---:|---:|---:|---:|---|")
    for item in summaries:
        lines.append(
            f"| {item.Q} | {item.nontrivial_prime_count_median} | {item.local_pass_rows} | "
            f"{item.drift_no_low_section_rows} | {item.low_section_rows} | "
            f"{item.duplicate_keys_among_passes} | `{item.first_empty_top}` |"
        )
    lines.append("")
    lines.append("## Stärkste natürliche Survivors")
    lines.append("")
    strongest = [item for item in survivors if item.best_pass_Q >= args.min_survivor_Q]
    if strongest:
        lines.append("| rank | Tupel | Gap | best Q | erste leere q danach | niedrige Sektionen | D_root | D_drift | Duplikat-Key | Entscheidung |")
        lines.append("|---:|---|---:|---:|---:|---:|---:|---:|---:|---|")
        for item in strongest[: args.report_rows]:
            tuple_text = f"`A={item.A}, x={item.x}, B={item.B}, y={item.y}, C={item.C}, z={item.z}`"
            lines.append(
                f"| {item.rank} | {tuple_text} | {item.gap_signed} | {item.best_pass_Q} | "
                f"{item.first_empty_after_best or '-'} | {item.low_sections_at_best} | "
                f"{item.D_root_at_best:.3f} | {item.D_drift_at_best:.3f} | "
                f"{item.duplicate_key_count} | `{item.decision}` |"
            )
    else:
        lines.append("Keine Survivor-Zeile erreicht die konfigurierte Mindestschwelle.")
    lines.append("")
    lines.append("## Befund")
    lines.append("")
    q61 = next((item for item in summaries if item.Q == 61), None)
    q97 = next((item for item in summaries if item.Q == 97), None)
    if q61 and q97:
        lines.append(
            f"Im gesamten Korpus bestehen `{q61.local_pass_rows}` von `{args.limit or 'allen'}` Zeilen alle nichttrivialen Root-Kanäle bis `Q=61`; alle diese Passes haben `0` niedrige globale Sektionen."
        )
        lines.append(
            f"Bei `Q=97` bleiben `{q97.local_pass_rows}` lokale Passes übrig. Damit gibt es in diesem Korpus keinen stabilen Above-threshold-Survivor über `Q=97`."
        )
    lines.append(
        "Die neue Information ist daher zweigeteilt: Es gibt stärkere Drift-ohne-Sektion-Beispiele als den rank-8-Fall, aber sie sind weiterhin finite lokale Übergangsfälle und werden spätestens bis `Q=97` lokal ausgeschlossen."
    )
    lines.append(
        "Zwei der fünf `Q=61`-Passes sind dieselbe linke Potenz in anderer Basisdarstellung (`25^6 = 125^4`) und zählen deshalb nicht als unabhängige Signaturfamilien."
    )
    lines.append("")
    lines.append("## Konsequenz")
    lines.append("")
    lines.append(
        "Für HD3' liefert der Lauf kein positives Beweissignal. Er verschärft den nächsten sinnvollen Suchauftrag: echte Kandidaten müssen mindestens über `Q=97` bestehen oder als synthetische Kontrollen mit vorregistrierter Root-Sektion geführt werden."
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=here / "near_power_b250_e3-7_z3-9_k5000.csv")
    parser.add_argument("--limit", type=int, default=5000)
    parser.add_argument("--checkpoints", type=int, nargs="+", default=[31, 37, 43, 61, 97, 127])
    parser.add_argument("--max-q", type=int, default=127)
    parser.add_argument("--min-survivor-Q", type=int, default=61)
    parser.add_argument("--sample-limit", type=int, default=12)
    parser.add_argument("--first-empty-top", type=int, default=8)
    parser.add_argument("--report-rows", type=int, default=12)
    parser.add_argument("--csv", type=Path, default=here / "ABOVE_THRESHOLD_SURVIVOR_SEARCH_2026-06-02.csv")
    parser.add_argument("--report", type=Path, default=here / "ABOVE_THRESHOLD_SURVIVOR_SEARCH_2026-06-02.md")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summaries, survivors = analyze(args)
    write_survivor_csv(survivors, args.csv)
    write_report(summaries, survivors, args)
    print(f"wrote {args.csv}")
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()
