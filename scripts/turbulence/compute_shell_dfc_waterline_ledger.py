"""
compute_shell_dfc_waterline_ledger.py
=====================================

Fast Sabra-shell smoke test for the DFC1^vee waterline ledger.

This script is not a DNS run and not a Navier-Stokes proof.  It extends the
existing toy DFC1^vee ledger with one short dynamical shell-model smoke and
matched flux-placement controls.  The purpose is guardrail bookkeeping:
positive mean flux is not enough unless the fixed window, monotone weights and
Abel-weighted placement all pass.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sabra_flux import sabra_energy_flux  # noqa: E402


# Result identity: model version + date of the corrected-dynamics run.  The
# 2026-06-05 snapshots were computed with the faulty Sabra dynamics (before
# Issue #1) and are superseded; see _results/SUPERSEDED.md.  Do not reuse the
# old tag, otherwise a new run would overwrite them under the same identity.
MODEL_TAG = "sabra-v2"
DATE_TAG = f"2026-10-04_{MODEL_TAG}"
OUT_DIR = Path(__file__).with_name("_results")

N_SHELLS = 18
LAMBDA = 2.0
K0 = 2**-4
NU = 1e-6
DT = 2e-4
T_THERM = 3.0
T_DATA = 8.0
SAMPLE_EVERY = 10
F_SHELL = 2
F_AMP = 2e-2
SEED = 42
C_PROJ_BUDGET = 0.15


@dataclass
class SimulationSummary:
    n_shells: int
    snapshots: int
    dt: float
    t_therm: float
    t_data: float
    nu: float
    forcing_shell: int
    forcing_amplitude: float
    mean_total_energy: float
    mean_dissipation: float
    runtime_seconds: float


@dataclass
class LedgerRow:
    scenario: str
    source_kind: str
    shell_window: str
    active_flux_constraints: str
    fixed_cascade_waterline: str
    dfc2_status: str
    pointwise_dfc1_status: str
    dual_dfc1_status: str
    weighted_flux_ratio: float | None
    projection_residual: float | None
    boundary_allowance: float | None
    residual_over_allowance: float | None
    bad_corridor_capacity: float | None
    worst_backscatter_shell: int | None
    same_mean_reference_risk: str
    verdict: str
    note: str


def sabra_nonlinear(u: np.ndarray, k: np.ndarray) -> np.ndarray:
    """Fix (2026-10-03, Issue #1): conjugation corrected in the local term
    (conj on u_{n-1} instead of u_{n+1}) and the backward-term coefficient
    corrected from -0.25 to +0.5, matching L'vov et al. 1998 Eq.(21) with
    a=1, b=-0.5, c=-0.5 (a+b+c=0). See compute_goy_shell_dfc.py::sabra_nonlinear
    for the full derivation and the counterexample from the issue."""
    out = np.zeros_like(u)
    out[: N_SHELLS - 2] += k[: N_SHELLS - 2] * np.conj(u[1 : N_SHELLS - 1]) * u[2:N_SHELLS]
    out[1 : N_SHELLS - 1] += -0.5 * k[: N_SHELLS - 2] * u[2:N_SHELLS] * np.conj(u[: N_SHELLS - 2])
    out[2:N_SHELLS] += 0.5 * k[: N_SHELLS - 2] * u[1 : N_SHELLS - 1] * u[: N_SHELLS - 2]
    return 1j * out


def imex_step(u: np.ndarray, k: np.ndarray) -> np.ndarray:
    forcing = np.zeros(N_SHELLS, dtype=complex)
    forcing[F_SHELL] = F_AMP
    return (u + DT * sabra_nonlinear(u, k) + DT * forcing) / (1.0 + DT * NU * k * k)


def run_sabra_smoke() -> tuple[np.ndarray, np.ndarray, SimulationSummary]:
    started = time.perf_counter()
    k = K0 * LAMBDA ** np.arange(N_SHELLS, dtype=float)
    rng = np.random.default_rng(SEED)
    u = np.array(
        [
            1e-4 * k[n] ** (-1.0 / 3.0) * np.exp(1j * 2.0 * math.pi * rng.random())
            for n in range(N_SHELLS)
        ],
        dtype=complex,
    )

    for _ in range(int(T_THERM / DT)):
        u = imex_step(u, k)

    energies: list[np.ndarray] = []
    fluxes: list[np.ndarray] = []
    dissipations: list[float] = []
    for step in range(int(T_DATA / DT)):
        u = imex_step(u, k)
        if step % SAMPLE_EVERY != 0:
            continue

        energy = np.abs(u) ** 2
        # Energiefluss: energiekonsistente Pi_n (Issue #1 Folgefix, 2026-10-03;
        # siehe sabra_flux.sabra_energy_flux). Die alte Formel
        # Im(k_n*u_n*conj(u_{n+1})*u_{n+2}) gehoerte zur VOR der Dynamik-
        # Korrektur verwendeten (fehlerhaften) Nichtlinearitaet.
        # Randkorrektur (unabhaengiger Review, 2026-10-03): sabra_energy_flux
        # liefert bereits alle N_SHELLS-1 inneren Grenzfluesse -- kein
        # Auffuellen mit 0 an der (falschen) Stelle N_SHELLS-2 mehr noetig.
        flux = sabra_energy_flux(u, k)
        energies.append(energy)
        fluxes.append(flux)
        dissipations.append(float(2.0 * NU * np.sum(k * k * energy)))

    energy_samples = np.asarray(energies)
    flux_samples = np.asarray(fluxes)
    runtime = time.perf_counter() - started
    summary = SimulationSummary(
        n_shells=N_SHELLS,
        snapshots=int(len(energy_samples)),
        dt=DT,
        t_therm=T_THERM,
        t_data=T_DATA,
        nu=NU,
        forcing_shell=F_SHELL,
        forcing_amplitude=F_AMP,
        mean_total_energy=float(np.mean(np.sum(energy_samples, axis=1))),
        mean_dissipation=float(np.mean(dissipations)),
        runtime_seconds=float(runtime),
    )
    return energy_samples, flux_samples, summary


def window_phi(energy_mean: np.ndarray, start: int, end: int) -> np.ndarray:
    k = K0 * LAMBDA ** np.arange(N_SHELLS, dtype=float)
    idx = np.arange(start, end + 1)
    e_ref = k[idx] ** (-2.0 / 3.0)
    e_ref *= float(np.sum(energy_mean[idx]) / np.sum(e_ref))
    return np.log(np.maximum(energy_mean[idx], 1e-300) / np.maximum(e_ref, 1e-300))


def blocked_row(name: str, start: int, end: int, note: str) -> LedgerRow:
    return LedgerRow(
        scenario=name,
        source_kind="sabra_smoke",
        shell_window=f"{start}-{end}",
        active_flux_constraints="fixed_window,no_posthoc_front",
        fixed_cascade_waterline="blocked_nonmonotone_phi",
        dfc2_status="fail",
        pointwise_dfc1_status="not_evaluated",
        dual_dfc1_status="blocked",
        weighted_flux_ratio=None,
        projection_residual=None,
        boundary_allowance=None,
        residual_over_allowance=None,
        bad_corridor_capacity=None,
        worst_backscatter_shell=None,
        same_mean_reference_risk="high",
        verdict="blocked: nonmonotone Free-Energy weights",
        note=note,
    )


def evaluate_window(
    name: str,
    source_kind: str,
    start: int,
    end: int,
    phi: np.ndarray,
    flux_mean: np.ndarray,
    flux_transform: str,
    note: str,
) -> LedgerRow:
    weights = phi[:-1] - phi[1:]
    if np.any(weights < -1e-12):
        return blocked_row(name, start, end, note)

    raw_flux = np.asarray(flux_mean[start:end], dtype=float)
    positive = raw_flux[raw_flux > 0]
    if len(positive) == 0:
        normalized_flux = raw_flux
        pi0 = 1.0
    else:
        pi0 = float(np.median(positive))
        normalized_flux = raw_flux / pi0 if pi0 else raw_flux

    if flux_transform == "actual":
        test_flux = normalized_flux
        same_mean_risk = "actual_order"
    elif flux_transform == "high_weight_sorted":
        test_flux = np.asarray(sorted(normalized_flux), dtype=float)
        same_mean_risk = "artificial_high_weight_pass_control"
    elif flux_transform == "low_weight_sorted":
        test_flux = np.asarray(sorted(normalized_flux, reverse=True), dtype=float)
        same_mean_risk = "matched_low_weight_fail_control"
    else:
        raise ValueError(f"unknown flux transform: {flux_transform}")

    target = float(np.sum(weights))
    weighted_flux = float(np.dot(weights, test_flux))
    edge = float(abs(phi[0]) + abs(phi[-1]))
    boundary_allowance = C_PROJ_BUDGET * edge
    projection_residual = max(0.0, target - weighted_flux)
    residual_ratio = projection_residual / boundary_allowance if boundary_allowance else None
    weighted_ratio = weighted_flux / target if target else None

    deficits = np.maximum(0.0, 1.0 - test_flux) * weights
    bad_capacity = float(np.sum(deficits) / target) if target else None
    worst_shell = None
    if len(deficits) and float(np.max(deficits)) > 0:
        worst_shell = start + int(np.argmax(deficits))

    dual_pass = projection_residual <= boundary_allowance + 1e-12
    pointwise_pass = bool(np.all(test_flux >= 1.0 - 1e-12))
    if source_kind == "sabra_matched_control" and dual_pass:
        verdict = "pass control: same flux multiset only after artificial placement"
    elif dual_pass:
        verdict = "pass"
    else:
        verdict = "fail: Abel-weighted waterline deficit"

    return LedgerRow(
        scenario=name,
        source_kind=source_kind,
        shell_window=f"{start}-{end}",
        active_flux_constraints="fixed_window,median_positive_flux_normalization,no_posthoc_front",
        fixed_cascade_waterline="tail_waterline_8_14",
        dfc2_status="pass",
        pointwise_dfc1_status="pass" if pointwise_pass else "fail",
        dual_dfc1_status="pass" if dual_pass else "fail",
        weighted_flux_ratio=weighted_ratio,
        projection_residual=projection_residual,
        boundary_allowance=boundary_allowance,
        residual_over_allowance=residual_ratio,
        bad_corridor_capacity=bad_capacity,
        worst_backscatter_shell=worst_shell,
        same_mean_reference_risk=same_mean_risk,
        verdict=verdict,
        note=note,
    )


def fmt(value: float | None) -> str:
    if value is None:
        return "n/a"
    return f"{value:.6g}"


def build_findings(rows: list[LedgerRow]) -> list[str]:
    """Findings text derived from the computed rows (no hard-coded results)."""
    by_name = {row.scenario: row for row in rows}
    transition = by_name["sabra_transition_window_4_12"]
    actual = by_name["sabra_tail_window_actual_order_8_14"]
    high = by_name["sabra_tail_window_sorted_high_weight_control"]
    low = by_name["sabra_tail_window_sorted_low_weight_control"]
    lines: list[str] = []

    if transition.dual_dfc1_status == "blocked":
        lines.extend(
            [
                "- Das feste Übergangsfenster `4-12` wird blockiert, weil die",
                "  Free-Energy-Gewichte nicht monoton sind. Das verhindert einen",
                "  nachträglichen Kaskadenfront-Claim.",
            ]
        )
    else:
        lines.extend(
            [
                "- Das feste Übergangsfenster `4-12` wird nicht blockiert; Urteil:",
                f"  `{transition.verdict}` (`weighted/target={fmt(transition.weighted_flux_ratio)}`,",
                f"  `residual/allowance={fmt(transition.residual_over_allowance)}`).",
            ]
        )

    figures = (
        f"`weighted/target={fmt(actual.weighted_flux_ratio)}`, "
        f"`residual/allowance={fmt(actual.residual_over_allowance)}`"
    )
    if actual.dual_dfc1_status == "blocked":
        lines.append("- Im Tail-Waterline-Fenster `8-14` sind die Gewichte nicht monoton (blockiert).")
    elif actual.dual_dfc1_status == "pass":
        lines.extend(
            [
                "- Im Tail-Waterline-Fenster `8-14` sind die Gewichte monoton und die",
                f"  tatsächliche Flux-Platzierung besteht `DFC1^vee` ({figures}).",
            ]
        )
    else:
        lines.extend(
            [
                "- Im Tail-Waterline-Fenster `8-14` sind die Gewichte monoton, aber die",
                f"  tatsächliche Flux-Platzierung besteht `DFC1^vee` nicht ({figures}).",
            ]
        )

    lines.append(
        f"- Sortierte High-Weight-Kontrolle: `{high.dual_dfc1_status}` "
        f"(`weighted/target={fmt(high.weighted_flux_ratio)}`, "
        f"`residual/allowance={fmt(high.residual_over_allowance)}`); "
        f"sortierte Low-Weight-Kontrolle: `{low.dual_dfc1_status}` "
        f"(`weighted/target={fmt(low.weighted_flux_ratio)}`, "
        f"`residual/allowance={fmt(low.residual_over_allowance)}`)."
    )
    if (
        actual.weighted_flux_ratio is not None
        and high.weighted_flux_ratio is not None
        and low.weighted_flux_ratio is not None
    ):
        if low.weighted_flux_ratio < actual.weighted_flux_ratio < high.weighted_flux_ratio:
            lines.append(
                "  Dasselbe Flux-Multiset liefert je nach Platzierung deutlich verschiedene"
                " Ergebnisse: Nicht der Mittelwert entscheidet, sondern die feste"
                " Abel-gewichtete Platzierung des Flux."
            )
        else:
            lines.append(
                "  Die erwartete Ordnung low < actual < high der gewichteten Flux-Verhältnisse"
                " wird in diesem Lauf nicht beobachtet; die Platzierungsaussage ist hier"
                " nicht belegt."
            )
    return lines


def build_markdown(rows: list[LedgerRow], summary: SimulationSummary) -> str:
    now = datetime.now().astimezone().isoformat(timespec="seconds")
    lines = [
        f"# Shell-DFC-Waterline-Ledger {DATE_TAG}",
        "",
        "> Modellstand `sabra-v2`: korrigierte Sabra-Dynamik und energiekonsistenter Fluss",
        "> (Issue #1). Überholt die Snapshots vom 2026-06-05 (fehlerhafte Dynamik).",
        "",
        f"Abschlusszeit: `{now}`",
        "",
        "## Zweck",
        "",
        "Dieser Lauf erweitert den bisherigen Toy-`DFC1^vee`-Ledger um einen",
        "kurzen dynamischen Sabra-Shell-Smoke und Matched Controls. Er ist kein",
        "DNS-Nachweis, kein Satz über Navier--Stokes und kein Upload-Gate.",
        "",
        "## Externer Paperstand",
        "",
        "- Benavides--Bustamante 2026 (`arXiv:2507.03397v2`) stützt Shell-Modelle",
        "  als Flux-Testbett, betont aber Phasendynamik als eigene Flussdeterminante.",
        "- Tuteri--Chibbaro--Alexakis 2026 (`arXiv:2603.11892`) zeigt, dass",
        "  Kaskadenrichtung und Shell-Geometrie bei 2D/dual-cascade-Modellen nicht",
        "  blind aus klassischen Shell-Modellen übertragen werden dürfen.",
        "- Zinchenko--Schumacher 2026 (`arXiv:2508.03401`) erweitert Duchon--Robert",
        "  auf kompressible Flüsse; das ist Kontext, aber kein direkter Eingang für",
        "  die inkompressible Paper-B-Bridge.",
        "- JHTDB Forced Isotropic Turbulence bleibt das nächste sinnvolle DNS-Ziel:",
        "  `1024^3`, `R_lambda=418`, `epsilon=0.103` und gespeicherte Velocity-/Pressure-Frames.",
        "",
        "## Sabra-Smoke-Parameter",
        "",
        "| Parameter | Wert |",
        "|---|---:|",
        f"| Shells | {summary.n_shells} |",
        f"| Snapshots | {summary.snapshots} |",
        f"| `dt` | {summary.dt:.6g} |",
        f"| Thermalisierung | {summary.t_therm:.6g} |",
        f"| Datenzeit | {summary.t_data:.6g} |",
        f"| `nu` | {summary.nu:.6g} |",
        f"| Forcing-Schale | {summary.forcing_shell} |",
        f"| Forcing-Amplitude | {summary.forcing_amplitude:.6g} |",
        f"| mittlere Energie | {summary.mean_total_energy:.6g} |",
        f"| mittlere Dissipation | {summary.mean_dissipation:.6g} |",
        f"| Laufzeit s | {summary.runtime_seconds:.3f} |",
        "",
        "## Ergebnis-Tabelle",
        "",
        "| Szenario | Fenster | DFC2 | DFC1^vee | weighted/target | residual/allowance | bad corridor | Urteil |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row.scenario,
                    row.shell_window,
                    row.dfc2_status,
                    row.dual_dfc1_status,
                    fmt(row.weighted_flux_ratio),
                    fmt(row.residual_over_allowance),
                    fmt(row.bad_corridor_capacity),
                    row.verdict,
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Befund",
            "",
            *build_findings(rows),
            "",
            "## Konsequenz für den Beweisstand",
            "",
            "`DFC1^vee` bleibt offen. Der neue Lauf stärkt die Guardrails: Shell- oder",
            "DNS-Daten zählen erst als Evidenz, wenn Fenster, Waterline, Phase-/Flux-",
            "Provenienz, Matched Controls und Projektionsrest vorab festgelegt sind.",
            "Ein v1.7-Upload sollte daraus keinen stärkeren Beweisclaim ableiten.",
            "",
            "## Maschinenlesbare Begleitdaten",
            "",
            f"- `SHELL_DFC_WATERLINE_LEDGER_{DATE_TAG}.json`",
            f"- `SHELL_DFC_WATERLINE_LEDGER_{DATE_TAG}.csv`",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    energy_samples, flux_samples, summary = run_sabra_smoke()
    energy_mean = np.mean(energy_samples, axis=0)
    flux_mean = np.mean(flux_samples, axis=0)

    transition_phi = window_phi(energy_mean, 4, 12)
    tail_phi = window_phi(energy_mean, 8, 14)

    rows = [
        blocked_row(
            "sabra_transition_window_4_12",
            4,
            12,
            "fixed transition window fails before flux evaluation because phi weights are nonmonotone",
        )
        if np.any(transition_phi[:-1] - transition_phi[1:] < -1e-12)
        else evaluate_window(
            "sabra_transition_window_4_12",
            "sabra_smoke",
            4,
            12,
            transition_phi,
            flux_mean,
            "actual",
            "fixed transition window",
        ),
        evaluate_window(
            "sabra_tail_window_actual_order_8_14",
            "sabra_smoke",
            8,
            14,
            tail_phi,
            flux_mean,
            "actual",
            "actual flux order in the fixed monotone tail waterline",
        ),
        evaluate_window(
            "sabra_tail_window_sorted_high_weight_control",
            "sabra_matched_control",
            8,
            14,
            tail_phi,
            flux_mean,
            "high_weight_sorted",
            "same flux multiset sorted so larger flux lands on larger Abel weights",
        ),
        evaluate_window(
            "sabra_tail_window_sorted_low_weight_control",
            "sabra_matched_control",
            8,
            14,
            tail_phi,
            flux_mean,
            "low_weight_sorted",
            "same flux multiset sorted so larger flux lands on smaller Abel weights",
        ),
    ]

    OUT_DIR.mkdir(exist_ok=True)
    json_path = OUT_DIR / f"SHELL_DFC_WATERLINE_LEDGER_{DATE_TAG}.json"
    csv_path = OUT_DIR / f"SHELL_DFC_WATERLINE_LEDGER_{DATE_TAG}.csv"
    md_path = OUT_DIR / f"SHELL_DFC_WATERLINE_LEDGER_{DATE_TAG}.md"

    payload = {
        "summary": asdict(summary),
        "rows": [asdict(row) for row in rows],
    }
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    md_path.write_text(build_markdown(rows, summary), encoding="utf-8")

    print(f"Wrote {json_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")
    for row in rows:
        print(
            f"{row.scenario}: dual={row.dual_dfc1_status}, "
            f"ratio={fmt(row.weighted_flux_ratio)}, verdict={row.verdict}"
        )


if __name__ == "__main__":
    main()
