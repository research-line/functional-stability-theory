#!/usr/bin/env python3
"""Matched open-drive control for the FST-II homochirality claim boundary.

The model is deliberately minimal and dimensionless.  It asks one narrow
question: can symmetric open driving plus autocatalysis amplify enantiomeric
excess without Frank-type mutual antagonism?  The answer is checked both by an
exact invariant and by a dependency-free RK4 sweep.  This is a model guardrail,
not empirical evidence for prebiotic homochirality.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


RUN_DATE = "2026-07-11"
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "_results"
OUT_BASE = f"HOMOCHIRALITY_OPEN_DRIVE_GATE_{RUN_DATE}"


@dataclass(frozen=True)
class Parameters:
    feed_concentration: float
    dilution_rate: float = 1.0
    autocatalysis_rate: float = 2.0
    antagonism_rate: float = 0.0


@dataclass(frozen=True)
class SimulationRow:
    case_id: str
    control_class: str
    feed_concentration: float
    antagonism_rate: float
    initial_ee: float
    final_ee: float
    absolute_ee_gain: float | None
    relative_ee_drift: float | None
    racemate_A_star: float
    racemate_L_star: float
    chiral_eigenvalue: float
    homochiral_invasion_eigenvalue: float
    final_A: float
    final_L: float
    final_D: float
    minimum_concentration: float
    decision: str


def rhs(state: tuple[float, float, float], p: Parameters) -> tuple[float, float, float]:
    """Dimensionless symmetric chemostat with optional mutual antagonism."""

    achiral, left, right = state
    q = p.dilution_rate
    k = p.autocatalysis_rate
    h = p.antagonism_rate
    production = k * achiral
    mutual_loss = h * left * right
    return (
        q * (p.feed_concentration - achiral) - production * (left + right),
        (production - q) * left - mutual_loss,
        (production - q) * right - mutual_loss,
    )


def add_scaled(
    state: tuple[float, float, float],
    slope: tuple[float, float, float],
    scale: float,
) -> tuple[float, float, float]:
    return tuple(value + scale * delta for value, delta in zip(state, slope, strict=True))


def rk4_step(
    state: tuple[float, float, float],
    dt: float,
    p: Parameters,
) -> tuple[float, float, float]:
    k1 = rhs(state, p)
    k2 = rhs(add_scaled(state, k1, 0.5 * dt), p)
    k3 = rhs(add_scaled(state, k2, 0.5 * dt), p)
    k4 = rhs(add_scaled(state, k3, dt), p)
    return tuple(
        value + (dt / 6.0) * (d1 + 2.0 * d2 + 2.0 * d3 + d4)
        for value, d1, d2, d3, d4 in zip(state, k1, k2, k3, k4, strict=True)
    )


def integrate(
    initial_state: tuple[float, float, float],
    p: Parameters,
    *,
    dt: float = 0.02,
    t_final: float = 600.0,
) -> tuple[tuple[float, float, float], float]:
    state = initial_state
    minimum = min(state)
    steps = int(round(t_final / dt))
    if not math.isclose(steps * dt, t_final, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("t_final must be an integer multiple of dt")

    for _ in range(steps):
        state = rk4_step(state, dt, p)
        minimum = min(minimum, *state)
        if not all(math.isfinite(value) for value in state):
            raise ArithmeticError("non-finite concentration encountered")
        if minimum < -1e-10:
            raise ArithmeticError(f"negative concentration encountered: {minimum}")

    return state, minimum


def enantiomeric_excess(left: float, right: float) -> float:
    total = left + right
    if total <= 0.0:
        return 0.0
    return (left - right) / total


def racemate_fixed_point(p: Parameters) -> tuple[float, float]:
    """Return (A*, x*) for the positive racemic branch L*=D*=x*."""

    q = p.dilution_rate
    k = p.autocatalysis_rate
    h = p.antagonism_rate
    threshold = q / k
    if p.feed_concentration <= threshold:
        raise ValueError("feed concentration is below the positive chiral branch")

    if h == 0.0:
        achiral = threshold
        x = 0.5 * (p.feed_concentration - achiral)
        return achiral, x

    # 2 k^2 A^2 + q(h-2k)A - h q A_in = 0; choose the positive root.
    b = q * (h - 2.0 * k)
    discriminant = b * b + 8.0 * k * k * h * q * p.feed_concentration
    achiral = (-b + math.sqrt(discriminant)) / (4.0 * k * k)
    x = (k * achiral - q) / h
    return achiral, x


def build_row(feed: float, antagonism: float, initial_ee: float) -> SimulationRow:
    p = Parameters(feed_concentration=feed, antagonism_rate=antagonism)
    threshold = p.dilution_rate / p.autocatalysis_rate
    initial_total = max(0.02, 0.2 * (feed - threshold))
    initial_left = 0.5 * initial_total * (1.0 + initial_ee)
    initial_right = 0.5 * initial_total * (1.0 - initial_ee)
    initial_state = (feed, initial_left, initial_right)

    final_state, minimum = integrate(initial_state, p)
    final_achiral, final_left, final_right = final_state
    final_ee = enantiomeric_excess(final_left, final_right)
    racemate_achiral, racemate_x = racemate_fixed_point(p)
    chiral_eigenvalue = p.autocatalysis_rate * racemate_achiral - p.dilution_rate
    homochiral_invasion = -antagonism * (feed - threshold)

    if initial_ee == 0.0:
        gain = None
        drift = abs(final_ee)
        decision = "exact_symmetric_trajectory"
    else:
        gain = abs(final_ee) / abs(initial_ee)
        drift = abs(final_ee - initial_ee) / abs(initial_ee)
        if antagonism == 0.0:
            decision = "negative_control_pass" if drift <= 1e-7 else "negative_control_fail"
        else:
            sign_preserved = math.copysign(1.0, final_ee) == math.copysign(1.0, initial_ee)
            decision = (
                "frank_positive_control_pass"
                if sign_preserved and abs(final_ee) >= 0.999
                else "frank_positive_control_fail"
            )

    label = "no_antagonism" if antagonism == 0.0 else "frank_antagonism"
    sign = "zero" if initial_ee == 0.0 else ("plus" if initial_ee > 0.0 else "minus")
    case_id = f"{label}_Ain_{feed:.1f}_ee_{sign}"
    return SimulationRow(
        case_id=case_id,
        control_class=label,
        feed_concentration=feed,
        antagonism_rate=antagonism,
        initial_ee=initial_ee,
        final_ee=final_ee,
        absolute_ee_gain=gain,
        relative_ee_drift=drift,
        racemate_A_star=racemate_achiral,
        racemate_L_star=racemate_x,
        chiral_eigenvalue=chiral_eigenvalue,
        homochiral_invasion_eigenvalue=homochiral_invasion,
        final_A=final_achiral,
        final_L=final_left,
        final_D=final_right,
        minimum_concentration=minimum,
        decision=decision,
    )


def build_rows() -> list[SimulationRow]:
    rows: list[SimulationRow] = []
    for feed in (0.6, 1.0, 1.4):
        for initial_ee in (-1e-4, 1e-4):
            rows.append(build_row(feed, 0.0, initial_ee))
            rows.append(build_row(feed, 2.0, initial_ee))
    rows.append(build_row(1.0, 2.0, 0.0))
    return rows


def finite_values(rows: Iterable[SimulationRow], field: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = getattr(row, field)
        if value is not None:
            values.append(value)
    return values


def validate(rows: list[SimulationRow]) -> dict[str, float | int | bool]:
    negative = [row for row in rows if row.control_class == "no_antagonism"]
    positive = [
        row
        for row in rows
        if row.control_class == "frank_antagonism" and row.initial_ee != 0.0
    ]
    exact = [row for row in rows if row.initial_ee == 0.0]

    max_negative_drift = max(finite_values(negative, "relative_ee_drift"))
    min_positive_final_ee = min(abs(row.final_ee) for row in positive)
    min_positive_gain = min(finite_values(positive, "absolute_ee_gain"))
    min_positive_chiral_eigenvalue = min(row.chiral_eigenvalue for row in positive)
    max_concentration_negativity = min(row.minimum_concentration for row in rows)

    summary: dict[str, float | int | bool] = {
        "rows": len(rows),
        "negative_control_rows": len(negative),
        "positive_control_rows": len(positive),
        "exact_symmetric_rows": len(exact),
        "max_no_antagonism_relative_ee_drift": max_negative_drift,
        "min_frank_final_abs_ee": min_positive_final_ee,
        "min_frank_absolute_ee_gain": min_positive_gain,
        "min_frank_chiral_eigenvalue": min_positive_chiral_eigenvalue,
        "minimum_concentration": max_concentration_negativity,
        "negative_control_passed": all(row.decision == "negative_control_pass" for row in negative),
        "frank_positive_control_passed": all(
            row.decision == "frank_positive_control_pass" for row in positive
        ),
        "exact_racemate_remains_symmetric": all(abs(row.final_ee) <= 1e-12 for row in exact),
        "claim_upgrade_allowed": False,
    }

    assert summary["negative_control_passed"]
    assert summary["frank_positive_control_passed"]
    assert summary["exact_racemate_remains_symmetric"]
    assert max_negative_drift <= 1e-7
    assert min_positive_final_ee >= 0.999
    assert min_positive_chiral_eigenvalue > 0.0
    assert max_concentration_negativity >= -1e-10
    return summary


def write_csv(rows: list[SimulationRow], path: Path) -> None:
    records = [asdict(row) for row in rows]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def format_number(value: float | None) -> str:
    if value is None:
        return "--"
    return f"{value:.6g}"


def write_markdown(
    rows: list[SimulationRow],
    summary: dict[str, float | int | bool],
    path: Path,
) -> None:
    lines = [
        f"# Homochiralitäts-Open-Drive-Gate — {RUN_DATE}",
        "",
        "**Status:** Minimalmodell-Guardrail bestanden; kein empirischer Claim-Upgrade.",
        "",
        "## Modell und exakte Wasserlinie",
        "",
        "Für achirales Substrat `A`, Enantiomere `L,D`, Feed/Dilution `q`, "
        "Autokatalyse `k` und gegenseitige Antagonismusrate `h` wird geprüft:",
        "",
        "```text",
        "A_dot = q(A_in-A) - k A(L+D)",
        "L_dot = (kA-q)L - hLD",
        "D_dot = (kA-q)D - hLD",
        "```",
        "",
        "Für `h=0` gilt bei `L,D>0` exakt `d/dt log(L/D)=0`; damit bleibt auch "
        "der Enantiomerenüberschuss `ee=(L-D)/(L+D)` konstant. Offener Feed kann "
        "die gesamte chirale Masse verändern, verstärkt aber ohne antagonistische "
        "Kopplung keine anfängliche chirale Asymmetrie.",
        "",
        "Auf der positiven racemischen Fixpunktverzweigung `L*=D*=x*>0` gilt für "
        "`h>0` dagegen `kA*-q=h x*`. Die antisymmetrische Mode `delta=L-D` erfüllt "
        "`delta_dot=(kA-q)delta`, also besitzt der racemische Fixpunkt die chirale "
        "Eigenrate `lambda_chi=h x*>0`. Das ist die Frank-Positivkontrolle. Ein exakt "
        "racemischer deterministischer Start bleibt dennoch exakt symmetrisch; die "
        "Zweigwahl benötigt eine Störung oder Rauschquelle.",
        "",
        "## Reproduzierbarer Sweep",
        "",
        "Parameter: `q=1`, `k=2`, `h in {0,2}`, `A_in in {0.6,1.0,1.4}`, "
        "`ee_0=±1e-4`; RK4 mit `dt=0.02`, `t_final=600`.",
        "",
        "| case | h | ee_0 | lambda_chi | ee_final | gain | decision |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row.case_id} | {row.antagonism_rate:.1f} | "
            f"{row.initial_ee:.6g} | {row.chiral_eigenvalue:.6g} | "
            f"{row.final_ee:.6g} | {format_number(row.absolute_ee_gain)} | "
            f"{row.decision} |"
        )

    lines.extend(
        [
            "",
            "## Zusammenfassung",
            "",
            f"- No-antagonism-Kontrollen: {summary['negative_control_rows']}; "
            f"maximale relative `ee`-Drift `{summary['max_no_antagonism_relative_ee_drift']:.3e}`.",
            f"- Frank-Positivkontrollen: {summary['positive_control_rows']}; "
            f"minimales finales `|ee|` `{summary['min_frank_final_abs_ee']:.6f}`.",
            f"- Kleinste positive chirale Eigenrate `{summary['min_frank_chiral_eigenvalue']:.6f}`.",
            f"- Exakt symmetrische Trajektorien: {summary['exact_symmetric_rows']}; "
            "die Racemat-Symmetrie blieb numerisch erhalten.",
            "- `claim_upgrade_allowed=false`.",
            "",
            "## Interpretation",
            "",
            "Das zuvor fehlende matched negative control ist im Minimalmodell geschlossen: "
            "Offenheit und symmetrische Autokatalyse allein erzeugen keine deterministische "
            "chirale Verstärkung. Die Frank-artige antagonistische Kopplung öffnet dagegen "
            "eine instabile chirale Mode und liefert im selben Rechenrahmen eine positive "
            "Kontrolle. Das schärft die P2-Wasserlinie, ersetzt aber weder Soai-/Viedma-Daten "
            "noch ein realistisches präbiotisches Reaktionsnetz und belegt keine universelle "
            "Homochiralitätsentstehung.",
            "",
            "**Skript:** `scripts/homochirality_open_drive_gate.py`",
            f"**Maschinenlesbare Daten:** `_results/{OUT_BASE}.json`, `_results/{OUT_BASE}.csv`",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    summary = validate(rows)
    payload = {
        "run_date": RUN_DATE,
        "project": "FST-II_Chemical",
        "model_scope": "dimensionless symmetric open chemostat guardrail",
        "exact_no_antagonism_invariant": "d/dt log(L/D)=0",
        "parameters": {
            "dilution_rate": 1.0,
            "autocatalysis_rate": 2.0,
            "feed_concentrations": [0.6, 1.0, 1.4],
            "antagonism_rates": [0.0, 2.0],
            "initial_ee": [-1e-4, 0.0, 1e-4],
            "dt": 0.02,
            "t_final": 600.0,
        },
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f"{OUT_BASE}.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    write_csv(rows, RESULTS / f"{OUT_BASE}.csv")
    write_markdown(rows, summary, RESULTS / f"{OUT_BASE}.md")

    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
