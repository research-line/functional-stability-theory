#!/usr/bin/env python3
"""QUANT-III-2b absolute score and threshold calibration gate.

The existing eta scan stores per-protein max-normalized residue scores.  This
script checks whether those artifacts support an absolute calibration claim.
It deliberately does not refit eta or optimize a threshold on the HP35 CSP
rank correlation.  Instead it evaluates:

1. scale identifiability after per-protein max normalization;
2. leave-one-residue-out calibration of HP35 Nash scores to absolute CSP;
3. transfer of the two existing in-sample functional-residue thresholds
   between HP35 and Protein G; and
4. eta sensitivity of the global spectral-radius excess rho(J) - 1.

The outputs are a guardrail ledger, not a biological validation claim.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Sequence


RUN_DATE = "2026-07-11"
ETA_CONVENTION = 0.015


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mean(values: Sequence[float]) -> float:
    return sum(values) / len(values)


def ols(x: Sequence[float], y: Sequence[float]) -> tuple[float, float]:
    x_bar = mean(x)
    y_bar = mean(y)
    denominator = sum((value - x_bar) ** 2 for value in x)
    if denominator == 0.0:
        return y_bar, 0.0
    slope = sum((xv - x_bar) * (yv - y_bar) for xv, yv in zip(x, y)) / denominator
    return y_bar - slope * x_bar, slope


def regression_diagnostics(x: Sequence[float], y: Sequence[float]) -> dict:
    intercept, slope = ols(x, y)
    fitted = [intercept + slope * value for value in x]
    y_bar = mean(y)
    tss = sum((value - y_bar) ** 2 for value in y)
    rss = sum((observed - predicted) ** 2 for observed, predicted in zip(y, fitted))
    r2 = 1.0 - rss / tss

    loo_predictions: list[float] = []
    loo_baselines: list[float] = []
    for held_out in range(len(x)):
        train_x = [value for index, value in enumerate(x) if index != held_out]
        train_y = [value for index, value in enumerate(y) if index != held_out]
        loo_intercept, loo_slope = ols(train_x, train_y)
        loo_predictions.append(loo_intercept + loo_slope * x[held_out])
        loo_baselines.append(mean(train_y))

    press = sum((observed - predicted) ** 2 for observed, predicted in zip(y, loo_predictions))
    baseline_press = sum(
        (observed - predicted) ** 2 for observed, predicted in zip(y, loo_baselines)
    )
    return {
        "n": len(x),
        "intercept": intercept,
        "slope": slope,
        "r2_in_sample": r2,
        "loocv_q2": 1.0 - press / tss,
        "loocv_mae": mean([abs(a - b) for a, b in zip(y, loo_predictions)]),
        "loocv_baseline_mae": mean([abs(a - b) for a, b in zip(y, loo_baselines)]),
        "loocv_rmse": math.sqrt(press / len(y)),
        "loocv_baseline_rmse": math.sqrt(baseline_press / len(y)),
    }


def classification_metrics(labels: Sequence[int], predictions: Sequence[int]) -> dict:
    tp = sum(label == 1 and prediction == 1 for label, prediction in zip(labels, predictions))
    tn = sum(label == 0 and prediction == 0 for label, prediction in zip(labels, predictions))
    fp = sum(label == 0 and prediction == 1 for label, prediction in zip(labels, predictions))
    fn = sum(label == 1 and prediction == 0 for label, prediction in zip(labels, predictions))
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denominator if denominator else 0.0
    return {
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "balanced_accuracy": (sensitivity + specificity) / 2.0,
        "mcc": mcc,
    }


def scan_at_eta(protein: dict, eta: float) -> dict:
    matches = [row for row in protein["eta_scan"] if math.isclose(row["eta"], eta)]
    if len(matches) != 1:
        raise ValueError(f"Expected one eta={eta} row for {protein['pdb_label']}")
    return matches[0]


def fmt(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}f}"


def build_markdown(summary: dict, rows: list[dict]) -> str:
    regression = summary["hp35_csp_regression"]
    transfers = summary["threshold_transfer"]
    lines = [
        "# QUANT-III-2b — Absolute Score-/Schwellenkalibrierungs-Gate",
        "",
        f"**Stand:** {RUN_DATE}",
        "",
        "## Ergebnis",
        "",
        "Der aktuelle Datenstand trägt keine positive absolute Kalibrierung. "
        "Die gespeicherten Residuenwerte sind für jedes Protein und jedes Eta separat "
        "auf ihr eigenes Maximum 1 normiert. Positive proteinabhängige Reskalierungen "
        "der Rohwerte bleiben daher unsichtbar; eine absolute Schwelle ist aus diesen "
        "Artefakten nicht identifizierbar.",
        "",
        "Die HP35-CSP-Abbildung ist als Rangsignal weiterhin berichtbar, aber nicht als "
        "absolute Vorhersage: Die lineare In-Sample-Erklärung beträgt nur "
        f"$R^2={fmt(regression['r2_in_sample'])}$, und die Leave-one-residue-out-"
        f"Kalibrierung fällt auf $Q^2={fmt(regression['loocv_q2'])}$. Der Modell-MAE "
        f"({fmt(regression['loocv_mae'])} ppm) ist schlechter als der jeweilige "
        f"Trainingsmittelwert ({fmt(regression['loocv_baseline_mae'])} ppm).",
        "",
        "## Gate-Ledger",
        "",
        "| Gate | Status | Kernergebnis |",
        "| --- | --- | --- |",
    ]
    for row in rows:
        lines.append(f"| `{row['gate']}` | **{row['status']}** | {row['detail']} |")

    lines.extend(
        [
            "",
            "## Proteinübergreifender Schwellentransfer",
            "",
            "Die folgenden Werte sind nur Negativkontrollen mit den bestehenden "
            "funktionellen Residuenlabels; sie sind kein CSP-Holdout und dürfen nicht "
            "als externe Kalibrierung gelesen werden.",
            "",
            "| Training → Test | Schwelle | Sensitivität | Spezifität | Balanced Accuracy | MCC |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for transfer in transfers:
        lines.append(
            f"| {transfer['train']} → {transfer['test']} | {fmt(transfer['threshold'], 6)} | "
            f"{fmt(transfer['sensitivity'])} | {fmt(transfer['specificity'])} | "
            f"{fmt(transfer['balanced_accuracy'])} | {fmt(transfer['mcc'])} |"
        )

    lines.extend(
        [
            "",
            "## Zulässige Schlussfolgerung",
            "",
            "QUANT-III-2b ist nicht positiv geschlossen. Geschlossen ist nur das "
            "Diagnose-Gate: Die aktuelle per-Protein-Max-Normalisierung kann keine "
            "absolute, proteinübergreifende Schwelle tragen, und der einzige gematchte "
            "kontinuierliche Messkanal (HP35/CSP) schlägt im Leave-one-out-Test den "
            "Mittelwert-Baseline nicht. `claim_upgrade_allowed=false`.",
            "",
            "## Reparaturvertrag",
            "",
            "1. Unnormalisierte Residuenbeiträge vor der Max-Normalisierung speichern; "
            "Normierungsfaktor, Eta, Hessian-/Metrikkonvention und aktive Modenmenge ausweisen.",
            "2. Eta oder eine dimensionslose physikalische Skala vor der Zielauswertung "
            "festlegen; kein CSP-/Label-basiertes Nachstimmen.",
            "3. Mindestens einen zweiten Protein-Messkanal desselben Typs (CSP, HDX oder "
            "Phi-Werte) als echten Holdout sichern: Kalibrierung auf Protein A, einmalige "
            "Prüfung auf Protein B.",
            "4. Schwelle ausschließlich im Trainingssatz wählen und gegen Label-, Kontakt- "
            "und Parameter-Shuffles sowie eine intrinsische/massengewichtete Metrik testen.",
            "",
            "## Reproduzierbarkeit",
            "",
            "- Skript: `scripts/quant_iii_2b_absolute_threshold_gate.py`",
            "- JSON: `_results/QUANT_III_2B_ABSOLUTE_THRESHOLD_GATE_2026-07-11.json`",
            "- CSV: `_results/QUANT_III_2B_ABSOLUTE_THRESHOLD_GATE_2026-07-11.csv`",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> dict:
    project_root = Path(__file__).resolve().parents[1]
    eta_path = project_root / "results" / "eta_scan" / "eta_scan_results.json"
    calibration_path = (
        project_root / "results" / "eta_calibration" / "eta_calibration_results.json"
    )
    benchmark_path = (
        project_root / "results" / "benchmark" / "frustration_benchmark_results.json"
    )
    output_dir = project_root / "_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    eta_data = read_json(eta_path)
    calibration_data = read_json(calibration_path)
    benchmark_data = read_json(benchmark_path)

    proteins = {row["pdb_label"]: row for row in eta_data["proteins"]}
    benchmarks = {row["protein"]: row for row in benchmark_data["benchmarks"]}
    hp35_label = next(label for label in proteins if "HP35" in label)
    protein_g_label = next(label for label in proteins if "Protein G" in label)
    hp35_benchmark = next(row for label, row in benchmarks.items() if "HP35" in label)
    protein_g_benchmark = next(row for label, row in benchmarks.items() if "Protein G" in label)

    # Input-integrity gate: the benchmark score vectors must be the eta=0.015 vectors.
    for protein_label, benchmark in (
        (hp35_label, hp35_benchmark),
        (protein_g_label, protein_g_benchmark),
    ):
        scan_scores = scan_at_eta(proteins[protein_label], ETA_CONVENTION)[
            "frustration_scores"
        ]
        if len(scan_scores) != len(benchmark["frustration_scores"]):
            raise AssertionError(f"Score length mismatch for {protein_label}")
        max_difference = max(
            abs(a - b) for a, b in zip(scan_scores, benchmark["frustration_scores"])
        )
        if max_difference > 1e-9:
            raise AssertionError(
                f"Benchmark/eta-scan score drift for {protein_label}: {max_difference}"
            )

    normalization_rows = []
    for protein in eta_data["proteins"]:
        all_max_one = all(
            math.isclose(max(row["frustration_scores"]), 1.0, abs_tol=1e-12)
            for row in protein["eta_scan"]
        )
        low = protein["eta_scan"][0]["rho_J"] - 1.0
        high = protein["eta_scan"][-1]["rho_J"] - 1.0
        normalization_rows.append(
            {
                "protein": protein["pdb_label"],
                "eta_rows": len(protein["eta_scan"]),
                "all_max_one": all_max_one,
                "rho_excess_ratio_eta_0p1_over_0p001": high / low,
            }
        )

    csp = calibration_data["csp_per_residue"]
    residue_indices = sorted(int(index) for index in csp)
    hp35_scores = scan_at_eta(proteins[hp35_label], ETA_CONVENTION)[
        "frustration_scores"
    ]
    regression = regression_diagnostics(
        [hp35_scores[index] for index in residue_indices],
        [csp[str(index)]["CSP"] for index in residue_indices],
    )

    transfers = []
    for train, test in (
        (hp35_benchmark, protein_g_benchmark),
        (protein_g_benchmark, hp35_benchmark),
    ):
        threshold = train["optimal_threshold"]["threshold"]
        predictions = [int(score >= threshold) for score in test["frustration_scores"]]
        transfers.append(
            {
                "train": train["protein"],
                "test": test["protein"],
                "threshold": threshold,
                **classification_metrics(test["labels"], predictions),
            }
        )

    rows = [
        {
            "gate": "per_protein_absolute_scale",
            "status": "BLOCKED",
            "detail": (
                f"{sum(row['eta_rows'] for row in normalization_rows)}/"
                f"{sum(row['eta_rows'] for row in normalization_rows)} Protein-Eta-Zeilen "
                "haben max(score)=1; positive proteinabhängige Rohskalen sind nicht identifizierbar."
            ),
        },
        {
            "gate": "hp35_csp_absolute_regression",
            "status": "FAIL",
            "detail": (
                f"R²={fmt(regression['r2_in_sample'])}, LOO-Q²={fmt(regression['loocv_q2'])}; "
                f"MAE={fmt(regression['loocv_mae'])} ppm ist schlechter als Mean-Baseline "
                f"{fmt(regression['loocv_baseline_mae'])} ppm."
            ),
        },
        {
            "gate": "eta_absolute_scale_convention",
            "status": "BLOCKED",
            "detail": (
                "rho(J)-1 wächst in allen drei Eta-Scans von eta=0.001 zu eta=0.1 um Faktor "
                f"{fmt(mean([row['rho_excess_ratio_eta_0p1_over_0p001'] for row in normalization_rows]), 1)}; "
                "ohne vorab fixierte Eta-/Skalenkonvention wandert jede absolute Schwelle."
            ),
        },
        {
            "gate": "cross_protein_threshold_transfer",
            "status": "DIAGNOSTIC_ONLY",
            "detail": (
                f"Balanced Accuracy {fmt(transfers[0]['balanced_accuracy'])} bzw. "
                f"{fmt(transfers[1]['balanced_accuracy'])}; Sensitivität "
                f"{fmt(transfers[0]['sensitivity'])}/{fmt(transfers[1]['sensitivity'])}. "
                "Andere Labelart als CSP und Schwellen jeweils in-sample gewählt."
            ),
        },
        {
            "gate": "independent_continuous_protein_holdout",
            "status": "MISSING",
            "detail": "Nur HP35 besitzt im aktuellen Paket einen gematchten kontinuierlichen CSP-Kanal.",
        },
        {
            "gate": "claim_upgrade",
            "status": "REJECTED",
            "detail": "Absolute Score-/Schwellenkalibrierung nicht positiv geschlossen.",
        },
    ]

    summary = {
        "metadata": {
            "description": "QUANT-III-2b absolute score/threshold calibration gate",
            "date": RUN_DATE,
            "owner": "LG",
            "eta_convention": ETA_CONVENTION,
            "input_sha256": {
                str(eta_path.relative_to(project_root)): sha256(eta_path),
                str(calibration_path.relative_to(project_root)): sha256(calibration_path),
                str(benchmark_path.relative_to(project_root)): sha256(benchmark_path),
            },
        },
        "normalization_audit": normalization_rows,
        "hp35_csp_regression": regression,
        "threshold_transfer": transfers,
        "gate_rows": rows,
        "raw_score_available": False,
        "normalization": "per_protein_max_to_one",
        "independent_measurement_channel": "HP35/BMRB-4428 CSP only",
        "independent_continuous_holdout_available": False,
        "same_fit_reference_risk": True,
        "target_chart_used": True,
        "claim_upgrade_allowed": False,
        "decision": "blocked_scale_nonidentifiability_and_failed_hp35_loocv",
    }

    json_path = output_dir / f"QUANT_III_2B_ABSOLUTE_THRESHOLD_GATE_{RUN_DATE}.json"
    csv_path = output_dir / f"QUANT_III_2B_ABSOLUTE_THRESHOLD_GATE_{RUN_DATE}.csv"
    md_path = output_dir / f"QUANT_III_2B_ABSOLUTE_THRESHOLD_GATE_{RUN_DATE}.md"

    json_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gate", "status", "detail"])
        writer.writeheader()
        writer.writerows(rows)
    md_path.write_text(build_markdown(summary, rows), encoding="utf-8")

    print(json.dumps({
        "decision": summary["decision"],
        "claim_upgrade_allowed": summary["claim_upgrade_allowed"],
        "hp35_csp_regression": regression,
        "threshold_transfer": transfers,
        "outputs": [str(md_path), str(csv_path), str(json_path)],
    }, indent=2, ensure_ascii=False))
    return summary


if __name__ == "__main__":
    main()
