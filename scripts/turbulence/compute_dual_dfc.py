"""
compute_dual_dfc.py
===================

Toy/ledger evaluator for the dual Downhill Flux Condition DFC1^vee.

This is not a DNS or Navier-Stokes proof.  It turns the paper's planned
DFC1^vee test into a reproducible bookkeeping check:

    sum_j a_j Pi_j >= Pi0 * sum_j a_j
                      - C_proj * (|phi_j1| + |phi_j2|)

with a_j = phi_j - phi_{j+1}.  The examples below are deliberately small:
positive controls, same-mean negative controls, and a corridor control where
pointwise flux looks acceptable on average but fails in the Abel-dual direction.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path


PI0 = 1.0
C_PROJ_BUDGET = 0.15
SHELLS = list(range(4, 14))
OUT_DIR = Path(__file__).with_name("_results")


@dataclass
class ScenarioResult:
    scenario: str
    description: str
    shells: str
    dfc2_status: str
    pointwise_dfc1_status: str
    dual_dfc1_status: str
    weighted_flux: float
    target_flux: float
    projection_residual: float
    boundary_allowance: float
    residual_over_allowance: float | None
    weighted_flux_ratio: float
    bad_corridor_capacity: float
    worst_backscatter_shell: int | None
    same_mean_reference_risk: str
    verdict: str


def phi_profile(kind: str, n: int) -> list[float]:
    if kind == "k41":
        return [0.0 for _ in range(n)]
    if kind == "smooth":
        return [0.9 * math.exp(-0.18 * i) + 0.08 for i in range(n)]
    if kind == "front_loaded":
        return [1.2 * math.exp(-0.33 * i) + 0.03 for i in range(n)]
    if kind == "slush":
        base = [0.72 * math.exp(-0.16 * i) + 0.04 for i in range(n)]
        base[3] += 0.08
        base[4] -= 0.04
        return base
    raise ValueError(f"unknown phi profile: {kind}")


def flux_profile(kind: str, n: int) -> list[float]:
    if kind == "k41":
        return [PI0 for _ in range(n - 1)]
    if kind == "smooth_forward":
        return [1.12 + 0.04 * math.sin(i) for i in range(n - 1)]
    if kind == "alternating_tolerated":
        return [1.65, 1.50, -0.22, 1.38, -0.18, 1.25, 1.20, 1.15, 1.10]
    if kind == "bad_corridor":
        return [-0.55, -0.35, 1.90, 1.85, 1.75, 1.70, 1.60, 1.55, 1.50]
    if kind == "same_mean_shuffle":
        return [1.73, -0.58, 1.62, -0.44, 1.50, 1.42, 1.34, 1.26, 1.18]
    raise ValueError(f"unknown flux profile: {kind}")


def monotone_weights(phi: list[float]) -> tuple[list[float], bool]:
    diffs = [phi[i] - phi[i + 1] for i in range(len(phi) - 1)]
    return diffs, all(x >= -1e-12 for x in diffs)


def evaluate(name: str, description: str, phi_kind: str, flux_kind: str) -> ScenarioResult:
    phi = phi_profile(phi_kind, len(SHELLS))
    flux = flux_profile(flux_kind, len(SHELLS))
    weights, monotone = monotone_weights(phi)

    weight_mass = sum(weights)
    edge = abs(phi[0]) + abs(phi[-1])
    weighted_flux = sum(a * p for a, p in zip(weights, flux))
    target = PI0 * weight_mass
    boundary_allowance = C_PROJ_BUDGET * edge
    projection_residual = max(0.0, target - weighted_flux)

    dual_pass = projection_residual <= boundary_allowance + 1e-12
    pointwise_pass = all(p >= PI0 for p in flux)

    deficits = [max(0.0, PI0 - p) * max(0.0, a) for a, p in zip(weights, flux)]
    bad_capacity = sum(deficits) / target if target > 0 else 0.0
    worst_shell = None
    if deficits and max(deficits) > 0:
        worst_shell = SHELLS[deficits.index(max(deficits))]

    if weight_mass == 0:
        ratio = 1.0
    else:
        ratio = weighted_flux / target

    residual_ratio = None
    if boundary_allowance > 0:
        residual_ratio = projection_residual / boundary_allowance

    if name in {"bad_corridor", "same_mean_shuffle"}:
        same_mean_risk = "high"
    elif name == "alternating_tolerated":
        same_mean_risk = "medium"
    else:
        same_mean_risk = "low"

    if weight_mass == 0:
        verdict = "vacuous reference control"
    elif dual_pass and not pointwise_pass:
        verdict = "dual pass despite pointwise backscatter"
    elif dual_pass:
        verdict = "pass"
    else:
        verdict = "fail: Abel-weighted backscatter corridor"

    return ScenarioResult(
        scenario=name,
        description=description,
        shells=f"{SHELLS[0]}-{SHELLS[-1]}",
        dfc2_status="pass" if monotone else "fail",
        pointwise_dfc1_status="pass" if pointwise_pass else "fail",
        dual_dfc1_status="pass" if dual_pass else "fail",
        weighted_flux=weighted_flux,
        target_flux=target,
        projection_residual=projection_residual,
        boundary_allowance=boundary_allowance,
        residual_over_allowance=residual_ratio,
        weighted_flux_ratio=ratio,
        bad_corridor_capacity=bad_capacity,
        worst_backscatter_shell=worst_shell,
        same_mean_reference_risk=same_mean_risk,
        verdict=verdict,
    )


def format_float(value: float | None) -> str:
    if value is None:
        return "n/a"
    return f"{value:.6g}"


def build_markdown(results: list[ScenarioResult]) -> str:
    now = datetime.now().astimezone().isoformat(timespec="seconds")
    lines = [
        "# Dual-DFC1-Ledger 2026-05-27",
        "",
        f"Abschlusszeit: `{now}`",
        "",
        "## Zweck",
        "",
        "Dieser Kurzlauf operationalisiert die geplante Prüfung von `DFC1^vee`",
        "für Turbulenz / Paper B. Er ist ein Ledger- und Negativkontrolltest,",
        "kein DNS-Nachweis und kein neuer Satz über Navier-Stokes.",
        "",
        "Getestet wird die Abel-duale Bedingung",
        "",
        "```text",
        "sum_j a_j Pi_j >= Pi0 * sum_j a_j - C_proj * (|phi_j1| + |phi_j2|)",
        "a_j = phi_j - phi_{j+1},  phi_j = log(E_j/E_j*)",
        "```",
        "",
        f"Parameter: `Pi0={PI0}`, `C_proj_budget={C_PROJ_BUDGET}`, ",
        f"Fenster: Shells `{SHELLS[0]}-{SHELLS[-1]}`.",
        "",
        "## Ergebnis-Tabelle",
        "",
        "| Szenario | DFC2 | punktweise DFC1 | DFC1^vee | weighted/target | residual/allowance | bad corridor | Urteil |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in results:
        lines.append(
            "| "
            + " | ".join(
                [
                    r.scenario,
                    r.dfc2_status,
                    r.pointwise_dfc1_status,
                    r.dual_dfc1_status,
                    format_float(r.weighted_flux_ratio),
                    format_float(r.residual_over_allowance),
                    format_float(r.bad_corridor_capacity),
                    r.verdict,
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Befund",
            "",
            "- `smooth_forward` ist die Positivkontrolle: punktweise und duale DFC1 bestehen.",
            "- `alternating_tolerated` zeigt den gewünschten Unterschied: einzelne negative Shell-Flüsse zerstören punktweise DFC1, aber nicht automatisch die Abel-duale Paarung.",
            "- `bad_corridor` und `same_mean_shuffle` sind Negativkontrollen: gleicher oder hoher mittlerer Fluss reicht nicht, wenn Backscatter genau in den stark gewichteten Abel-Shells sitzt.",
            "- Damit ist der nächste empirische Test nicht ein Mittelwert-Flux, sondern ein festes Fenster-/Waterline-Ledger mit `projection_residual`, `bad_corridor_capacity` und `same_mean_reference_risk`.",
            "",
            "## Konsequenz für den Beweisstand",
            "",
            "`DFC1^vee` bleibt offen und wird durch diesen Lauf nicht bewiesen. Der Fortschritt ist operativ: Die offene DR/Eyink -> LP/Wavelet-Projektionsbrücke ist jetzt in ein prüfbares Ledger zerlegt. Ein echter Nachweis braucht weiterhin Daten oder Analyse für den Projektionsfehler und eine quantitative DFC3-Schätzung.",
            "",
            "## Maschinenlesbarer Begleitdatensatz",
            "",
            "- `DUAL_DFC_LEDGER_2026-05-27.json`",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    scenarios = [
        (
            "k41_reference",
            "reference profile; all Abel weights vanish",
            "k41",
            "k41",
        ),
        (
            "smooth_forward",
            "positive control with forward flux on all weighted shells",
            "smooth",
            "smooth_forward",
        ),
        (
            "alternating_tolerated",
            "negative pointwise shell fluxes away from dominant Abel mass",
            "smooth",
            "alternating_tolerated",
        ),
        (
            "bad_corridor",
            "backscatter concentrated in high-weight frontier shells",
            "front_loaded",
            "bad_corridor",
        ),
        (
            "same_mean_shuffle",
            "similar mean-flux profile with badly placed backscatter",
            "slush",
            "same_mean_shuffle",
        ),
    ]

    results = [evaluate(*scenario) for scenario in scenarios]

    OUT_DIR.mkdir(exist_ok=True)
    json_path = OUT_DIR / "DUAL_DFC_LEDGER_2026-05-27.json"
    md_path = OUT_DIR / "DUAL_DFC_LEDGER_2026-05-27.md"

    json_path.write_text(
        json.dumps([asdict(r) for r in results], indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    md_path.write_text(build_markdown(results), encoding="utf-8")

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    for r in results:
        print(
            f"{r.scenario}: dual={r.dual_dfc1_status}, "
            f"pointwise={r.pointwise_dfc1_status}, "
            f"ratio={r.weighted_flux_ratio:.4f}, verdict={r.verdict}"
        )


if __name__ == "__main__":
    main()
