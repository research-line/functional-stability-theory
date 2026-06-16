"""
compute_os_capacity_ledger.py
=============================

Small follow-up to compute_birkhoff_rg.py.

Purpose:
- keep the old Kingman/Birkhoff diagnostic,
- add a separate OS-danger capacity proxy,
- include a negative control where mean contraction is negative but the
  dangerous capacity is not summable.

This is a ledger/prototype, not a Yang-Mills proof.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
RESULT_DIR = ROOT / "_results"
RESULT_DIR.mkdir(exist_ok=True)
DATA_DIR = ROOT / "_data"

V2_FIELDS = [
    "scenario",
    "k",
    "source_kind",
    "tail_model",
    "tau_B",
    "epsilon_safe",
    "eta_os_danger",
    "cap_os_path",
    "local_visible_defect",
    "nonlocal_tail_defect",
    "target_gap",
    "rp_cone_status",
    "circularity_status",
    "cofactor_certificate_status",
    "source_note",
]


def compute_transfer_matrix(beta: float, n_bins: int = 16) -> np.ndarray:
    a_vals = np.linspace(-0.99, 0.99, n_bins)
    da = a_vals[1] - a_vals[0]
    transfer = np.zeros((n_bins, n_bins), dtype=float)

    for i, ai in enumerate(a_vals):
        for j, aj in enumerate(a_vals):
            weight = math.exp(beta * ai * aj)
            haar_j = math.sqrt(max(0.0, 1.0 - aj * aj))
            transfer[i, j] = weight * haar_j * da

    rows = transfer.sum(axis=1)
    rows[rows == 0.0] = 1.0
    return transfer / rows[:, None]


def birkhoff_contraction(transfer: np.ndarray) -> tuple[float, float]:
    positive = np.maximum(transfer, 1e-300)
    log_t = np.log(positive)

    # Delta = max_{i,j,k,l} |log T_ik + log T_jl - log T_il - log T_jk|.
    delta = 0.0
    n = positive.shape[0]
    for i in range(n):
        for j in range(n):
            block = (
                log_t[i, :, None]
                + log_t[j, None, :]
                - log_t[i, None, :]
                - log_t[j, :, None]
            )
            delta = max(delta, float(np.max(np.abs(block))))

    return math.tanh(delta / 4.0), delta


def spectral_ratio_and_gap(transfer: np.ndarray) -> tuple[float, float]:
    eigvals = np.sort(np.abs(np.linalg.eigvals(transfer)))[::-1]
    if len(eigvals) < 2 or eigvals[0] <= 0:
        return 0.0, 1.0
    ratio = float(eigvals[1] / eigvals[0])
    return ratio, float(1.0 - ratio)


def rg_rows(beta: float, n_levels: int = 8, n_bins: int = 16) -> list[dict]:
    transfer = compute_transfer_matrix(beta, n_bins=n_bins)
    rows = []
    for k in range(n_levels):
        tau_b, diameter = birkhoff_contraction(transfer)
        tau_spectral, gap = spectral_ratio_and_gap(transfer)
        os_tail_proxy = max(0.0, tau_b - tau_spectral)

        safe_degradation = tau_spectral * tau_spectral / ((k + 1) ** 2)
        os_capacity = os_tail_proxy * os_tail_proxy / ((k + 1) ** 2)
        epsilon = min(0.49, safe_degradation + os_capacity)

        rows.append(
            {
                "k": k,
                "tau_B": tau_b,
                "tau_spectral": tau_spectral,
                "gap": gap,
                "diameter": diameter,
                "safe_degradation": safe_degradation,
                "os_tail_proxy": os_tail_proxy,
                "os_capacity": os_capacity,
                "epsilon_total": epsilon,
            }
        )

        transfer = transfer @ transfer
        rowsums = transfer.sum(axis=1)
        rowsums[rowsums == 0.0] = 1.0
        transfer = transfer / rowsums[:, None]

    return rows


def summarize(rows: list[dict]) -> dict:
    taus = np.array([row["tau_B"] for row in rows], dtype=float)
    eps = np.array([row["epsilon_total"] for row in rows], dtype=float)
    os_caps = np.array([row["os_capacity"] for row in rows], dtype=float)
    safe = np.array([row["safe_degradation"] for row in rows], dtype=float)
    bad_flags = [row["os_tail_proxy"] > 0.05 or row["tau_B"] > 0.995 for row in rows]

    longest_bad_run = 0
    current = 0
    for flag in bad_flags:
        current = current + 1 if flag else 0
        longest_bad_run = max(longest_bad_run, current)

    product_lower = 1.0
    for value in eps:
        product_lower *= max(0.0, 1.0 - float(value))

    return {
        "levels": len(rows),
        "mean_log_tau_B": float(np.mean(np.log(np.maximum(taus, 1e-300)))),
        "max_tau_B": float(np.max(taus)),
        "sum_safe_degradation": float(np.sum(safe)),
        "sum_os_capacity": float(np.sum(os_caps)),
        "gt2_partial_product_lower": float(product_lower),
        "longest_bad_run_proxy": int(longest_bad_run),
    }


def synthetic_controls(k_max: int = 5000) -> dict:
    ks = np.arange(1, k_max + 1, dtype=float)
    tau = np.exp(-0.08 + 0.025 * np.sin(np.log(ks + 1.0)))
    mean_log_tau = float(np.mean(np.log(tau)))

    positive_capacity = 0.25 / (ks * ks)
    negative_capacity = 0.25 / ks

    def checkpoint(capacity: np.ndarray, k: int) -> dict:
        partial = float(np.sum(capacity[:k]))
        return {
            "K": k,
            "sum_capacity": partial,
            "exp_minus_sum_capacity": float(math.exp(-partial)),
        }

    checkpoints = [50, 200, 1000, 5000]
    return {
        "mean_log_tau": mean_log_tau,
        "positive_summable": [checkpoint(positive_capacity, k) for k in checkpoints],
        "negative_harmonic": [checkpoint(negative_capacity, k) for k in checkpoints],
    }


def markdown_report(payload: dict) -> str:
    lines = [
        "# OS-Capacity-Ledger 2026-05-28",
        "",
        "Status: numerischer Ideencheck, kein Beweisclaim.",
        "",
        "Die Auswertung trennt drei Signale:",
        "",
        "- `mean_log_tau_B`: alte Kingman-/Birkhoff-Diagnostik.",
        "- `sum_safe_degradation`: lokal sichtbare, OS-sichere Degradation.",
        "- `sum_os_capacity`: Proxy für OS-gefährliche Korridorkapazität.",
        "",
        "## SU(2)-Toy-RG",
        "",
        "| beta | mean log tau_B | max tau_B | sum safe | sum OS-cap | GT2 product | bad-run |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for run in payload["runs"]:
        s = run["summary"]
        lines.append(
            f"| {run['beta']:.1f} | {s['mean_log_tau_B']:.6f} | "
            f"{s['max_tau_B']:.6f} | {s['sum_safe_degradation']:.6f} | "
            f"{s['sum_os_capacity']:.6f} | {s['gt2_partial_product_lower']:.6f} | "
            f"{s['longest_bad_run_proxy']} |"
        )

    lines.extend(
        [
            "",
            "## Negativkontrolle",
            "",
            "Die synthetische Kontrolle hält `mean(log tau_k) < 0` fest, variiert aber nur die",
            "OS-gefährliche Korridorkapazität. Bei `1/k^2` bleibt das Budget endlich; bei",
            "`1/k` wächst es harmonisch weiter. Genau dieser Fall ist der False-positive-Kanal",
            "der alten Kingman-Kurzform.",
            "",
            f"`mean_log_tau = {payload['synthetic_controls']['mean_log_tau']:.6f}`",
            "",
            "| Fall | K | Summe Kapazität | exp(-Summe) |",
            "|---|---:|---:|---:|",
        ]
    )

    for key, label in [
        ("positive_summable", "summierbar 1/k^2"),
        ("negative_harmonic", "nicht summierbar 1/k"),
    ]:
        for row in payload["synthetic_controls"][key]:
            lines.append(
                f"| {label} | {row['K']} | {row['sum_capacity']:.6f} | "
                f"{row['exp_minus_sum_capacity']:.6f} |"
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Der Test bestätigt nicht den Yang-Mills-Continuum-Transfer. Er macht aber den",
            "nächsten prüfbaren Schritt präziser: Ein Nachfolgeskript sollte echte",
            "`local-visible`- und `nonlocal-tail`-Defekte aus RG-Blockdaten statt aus diesem",
            "Proxy ableiten. Erst dann kann `Cap_K^OS` als mathematische Hypothese sinnvoll",
            "gegen U(1)-/Strong-Coupling-Positivkontrollen und synthetische Negativkontrollen",
            "getestet werden.",
            "",
        ]
    )
    return "\n".join(lines)


def generated_v2_control_rows(k_max: int = 500) -> list[dict]:
    rows: list[dict] = []

    def add_row(
        scenario: str,
        k: int,
        source_kind: str,
        tail_model: str,
        tau_b: float,
        epsilon_safe: float,
        eta_os_danger: float,
        cap_os_path: float,
        local_visible_defect: float,
        nonlocal_tail_defect: float,
        target_gap: float,
        rp_cone_status: str,
        circularity_status: str,
        cofactor_certificate_status: str,
        source_note: str,
    ) -> None:
        rows.append(
            {
                "scenario": scenario,
                "k": k,
                "source_kind": source_kind,
                "tail_model": tail_model,
                "tau_B": tau_b,
                "epsilon_safe": epsilon_safe,
                "eta_os_danger": eta_os_danger,
                "cap_os_path": cap_os_path,
                "local_visible_defect": local_visible_defect,
                "nonlocal_tail_defect": nonlocal_tail_defect,
                "target_gap": target_gap,
                "rp_cone_status": rp_cone_status,
                "circularity_status": circularity_status,
                "cofactor_certificate_status": cofactor_certificate_status,
                "source_note": source_note,
            }
        )

    for k in range(1, k_max + 1):
        add_row(
            scenario="strong_coupling_positive_control",
            k=k,
            source_kind="finite_lattice_strong_coupling",
            tail_model="summable",
            tau_b=math.exp(-0.22 - 0.02 / (k + 1.0)),
            epsilon_safe=0.020 / ((k + 1.0) ** 2),
            eta_os_danger=0.006 / ((k + 1.0) ** 2),
            cap_os_path=0.006 / ((k + 1.0) ** 2),
            local_visible_defect=0.018 / ((k + 1.0) ** 2),
            nonlocal_tail_defect=0.004 / ((k + 1.0) ** 2),
            target_gap=0.18,
            rp_cone_status="pass",
            circularity_status="pre_registered",
            cofactor_certificate_status="not_required_control",
            source_note="finite strong-coupling control; no continuum claim",
        )
        add_row(
            scenario="summable_os_capacity_control",
            k=k,
            source_kind="synthetic_rg_input_csv",
            tail_model="summable",
            tau_b=math.exp(-0.09 + 0.018 * math.sin(math.log(k + 2.0))),
            epsilon_safe=0.018 / ((k + 1.0) ** 2),
            eta_os_danger=0.030 / ((k + 1.0) ** 2),
            cap_os_path=0.030 / ((k + 1.0) ** 2),
            local_visible_defect=0.012 / ((k + 1.0) ** 2),
            nonlocal_tail_defect=0.025 / ((k + 1.0) ** 2),
            target_gap=0.10,
            rp_cone_status="pass",
            circularity_status="pre_registered",
            cofactor_certificate_status="toy_repair_present",
            source_note="summable OS-danger control, still synthetic",
        )
        add_row(
            scenario="kingman_false_positive_harmonic",
            k=k,
            source_kind="synthetic_negative_control",
            tail_model="harmonic",
            tau_b=math.exp(-0.08 + 0.025 * math.sin(math.log(k + 1.0))),
            epsilon_safe=0.012 / ((k + 1.0) ** 2),
            eta_os_danger=0.080 / k,
            cap_os_path=0.080 / k,
            local_visible_defect=0.010 / ((k + 1.0) ** 2),
            nonlocal_tail_defect=0.080 / k,
            target_gap=0.10,
            rp_cone_status="pass",
            circularity_status="pre_registered",
            cofactor_certificate_status="missing_repair",
            source_note="negative mean contraction with non-summable OS corridor",
        )
        add_row(
            scenario="rp_cone_fail_control",
            k=k,
            source_kind="synthetic_negative_control",
            tail_model="summable",
            tau_b=math.exp(-0.12),
            epsilon_safe=0.010 / ((k + 1.0) ** 2),
            eta_os_danger=0.010 / ((k + 1.0) ** 2),
            cap_os_path=0.010 / ((k + 1.0) ** 2),
            local_visible_defect=0.008 / ((k + 1.0) ** 2),
            nonlocal_tail_defect=0.006 / ((k + 1.0) ** 2),
            target_gap=0.11,
            rp_cone_status="fail",
            circularity_status="pre_registered",
            cofactor_certificate_status="irrelevant_until_rp_pass",
            source_note="summable tail but OS/RP cone compatibility deliberately fails",
        )
        add_row(
            scenario="circular_source_fail_control",
            k=k,
            source_kind="target_recycled_advice",
            tail_model="summable",
            tau_b=math.exp(-0.16),
            epsilon_safe=0.010 / ((k + 1.0) ** 2),
            eta_os_danger=0.004 / ((k + 1.0) ** 2),
            cap_os_path=0.004 / ((k + 1.0) ** 2),
            local_visible_defect=0.006 / ((k + 1.0) ** 2),
            nonlocal_tail_defect=0.003 / ((k + 1.0) ** 2),
            target_gap=0.16,
            rp_cone_status="pass",
            circularity_status="post_hoc_from_target_gap",
            cofactor_certificate_status="target_recycled",
            source_note="control for circular use of the desired continuum gap",
        )

    return rows


def write_input_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=V2_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in V2_FIELDS})


def read_input_csv(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = [field for field in V2_FIELDS if field not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Missing input columns: {', '.join(missing)}")
        for row in reader:
            rows.append(
                {
                    **{field: row[field] for field in V2_FIELDS},
                    "k": int(row["k"]),
                    "tau_B": float(row["tau_B"]),
                    "epsilon_safe": float(row["epsilon_safe"]),
                    "eta_os_danger": float(row["eta_os_danger"]),
                    "cap_os_path": float(row["cap_os_path"]),
                    "local_visible_defect": float(row["local_visible_defect"]),
                    "nonlocal_tail_defect": float(row["nonlocal_tail_defect"]),
                    "target_gap": float(row["target_gap"]),
                }
            )
    return rows


def summarize_v2_group(scenario: str, rows: list[dict]) -> dict:
    rows = sorted(rows, key=lambda item: item["k"])
    tau = np.array([row["tau_B"] for row in rows], dtype=float)
    eps = np.array([row["epsilon_safe"] for row in rows], dtype=float)
    eta = np.array([row["eta_os_danger"] for row in rows], dtype=float)
    cap = np.array([row["cap_os_path"] for row in rows], dtype=float)
    local_visible = np.array([row["local_visible_defect"] for row in rows], dtype=float)
    nonlocal_tail = np.array([row["nonlocal_tail_defect"] for row in rows], dtype=float)

    product_lower = 1.0
    for e_value, eta_value in zip(eps, eta):
        product_lower *= max(0.0, 1.0 - float(e_value)) * max(0.0, 1.0 - float(eta_value))

    rp_status = sorted({row["rp_cone_status"] for row in rows})
    circularity_status = sorted({row["circularity_status"] for row in rows})
    cofactor_status = sorted({row["cofactor_certificate_status"] for row in rows})
    tail_models = sorted({row["tail_model"] for row in rows})
    source_kinds = sorted({row["source_kind"] for row in rows})

    bad_flags = [
        row["eta_os_danger"] > 0.02
        or row["rp_cone_status"] != "pass"
        or row["circularity_status"] != "pre_registered"
        for row in rows
    ]
    longest_bad_run = 0
    current = 0
    for flag in bad_flags:
        current = current + 1 if flag else 0
        longest_bad_run = max(longest_bad_run, current)

    good_scale_density = float(
        sum(
            1
            for row in rows
            if row["rp_cone_status"] == "pass"
            and row["circularity_status"] == "pre_registered"
            and row["epsilon_safe"] + row["eta_os_danger"] < 0.02
        )
        / max(1, len(rows))
    )

    if any(status != "pre_registered" for status in circularity_status):
        decision = "rejected_circular_source"
    elif any(status != "pass" for status in rp_status):
        decision = "blocked_rp_or_os_cone"
    elif "harmonic" in tail_models:
        decision = "rejected_non_summable_os_capacity"
    elif any(status in {"missing_repair", "target_recycled"} for status in cofactor_status):
        decision = "blocked_missing_repair_certificate"
    else:
        decision = "control_pass_summable_no_claim"

    return {
        "scenario": scenario,
        "levels": len(rows),
        "source_kinds": source_kinds,
        "tail_models": tail_models,
        "mean_log_tau_B": float(np.mean(np.log(np.maximum(tau, 1e-300)))),
        "max_tau_B": float(np.max(tau)),
        "sum_epsilon_safe": float(np.sum(eps)),
        "sum_eta_os_danger": float(np.sum(eta)),
        "sum_cap_os_path": float(np.sum(cap)),
        "sum_local_visible_defect": float(np.sum(local_visible)),
        "sum_nonlocal_tail_defect": float(np.sum(nonlocal_tail)),
        "gt2_partial_product_lower": float(product_lower),
        "longest_bad_run": int(longest_bad_run),
        "good_scale_density": good_scale_density,
        "rp_cone_status": rp_status,
        "circularity_status": circularity_status,
        "cofactor_certificate_status": cofactor_status,
        "decision": decision,
    }


def v2_payload(rows: list[dict], source_path: Path | None, date_tag: str) -> dict:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[row["scenario"]].append(row)

    summaries = [
        summarize_v2_group(scenario, group_rows)
        for scenario, group_rows in sorted(grouped.items())
    ]

    return {
        "date": date_tag,
        "status": "ledger/control run only; no Yang-Mills mass-gap claim",
        "input_path": str(source_path) if source_path else None,
        "rows": rows,
        "summaries": summaries,
        "external_short_check": [
            {
                "source": "arXiv:2505.16585",
                "role": "area-law/truncated-model context; not a continuum proof",
                "url": "https://arxiv.org/abs/2505.16585",
            },
            {
                "source": "arXiv:2509.04688",
                "role": "area-law mass-gap-condition context; not a Clay closure",
                "url": "https://arxiv.org/abs/2509.04688",
            },
            {
                "source": "arXiv:2506.00284",
                "role": "withdrawn proof claim; rejected as project input",
                "url": "https://arxiv.org/abs/2506.00284",
            },
        ],
    }


def write_v2_summary_csv(summaries: list[dict], path: Path) -> None:
    fields = [
        "scenario",
        "levels",
        "mean_log_tau_B",
        "max_tau_B",
        "sum_epsilon_safe",
        "sum_eta_os_danger",
        "sum_cap_os_path",
        "gt2_partial_product_lower",
        "longest_bad_run",
        "good_scale_density",
        "rp_cone_status",
        "circularity_status",
        "cofactor_certificate_status",
        "decision",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for summary in summaries:
            writer.writerow(
                {
                    **{field: summary[field] for field in fields if field in summary},
                    "rp_cone_status": "|".join(summary["rp_cone_status"]),
                    "circularity_status": "|".join(summary["circularity_status"]),
                    "cofactor_certificate_status": "|".join(
                        summary["cofactor_certificate_status"]
                    ),
                }
            )


def markdown_v2_report(payload: dict) -> str:
    lines = [
        f"# OS-Capacity-Ledger v2 {payload['date']}",
        "",
        "Status: FORSCHER-Kontrolllauf, kein Yang-Mills-Beweis und kein Claim-Upgrade.",
        "",
        "Ziel des v2-Ledgers ist nicht ein weiterer Kingman-Mittelwert, sondern die",
        "sichtbare Trennung von Quelle, OS-sicherer Degradation, OS-gefährlicher",
        "Korridorkapazität, Zirkularität und Repair-/Cofactor-Zertifikat.",
        "",
        "## Externer Kurzcheck",
        "",
        "- `arXiv:2505.16585` und `arXiv:2509.04688` bleiben nützlich als",
        "  Area-Law-/Mass-Gap-Condition-Kontext für Lattice-Yang-Mills, schließen",
        "  aber den Kontinuums-/OS-Transfer dieses Projekts nicht.",
        "- `arXiv:2506.00284` behauptete einen konstruktiven SU(3)-Beweis, ist aber",
        "  von arXiv Admin zurückgezogen; es wird nicht als Projektnachweis genutzt.",
        "",
        "## Ledger-Entscheidungen",
        "",
        "| Szenario | mean log tau_B | Summe eps_safe | Summe eta_OS | GT2-Produkt | good-scale | RP | Zirkularität | Cofactor | Entscheidung |",
        "|---|---:|---:|---:|---:|---:|---|---|---|---|",
    ]

    for summary in payload["summaries"]:
        lines.append(
            f"| {summary['scenario']} | {summary['mean_log_tau_B']:.6f} | "
            f"{summary['sum_epsilon_safe']:.6f} | "
            f"{summary['sum_eta_os_danger']:.6f} | "
            f"{summary['gt2_partial_product_lower']:.6f} | "
            f"{summary['good_scale_density']:.3f} | "
            f"{'/'.join(summary['rp_cone_status'])} | "
            f"{'/'.join(summary['circularity_status'])} | "
            f"{'/'.join(summary['cofactor_certificate_status'])} | "
            f"{summary['decision']} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Der starke Kontrollfall und der summierbare synthetische Kontrollfall zeigen",
            "nur, dass die Ledger-Logik Positivkontrollen nicht automatisch verwirft. Der",
            "harmonische False-Positive wird trotz negativem `mean_log_tau_B` verworfen,",
            "weil die OS-gefährliche Kapazität nicht summierbar ist. Der RP-Kegel-Fail und",
            "der zirkuläre Quellen-Fail werden ebenfalls verworfen, obwohl ihre numerischen",
            "Summen harmlos aussehen.",
            "",
            "Damit ist der nächste Beweisschritt präziser, aber nicht geschlossen:",
            "Ein echter Yang-Mills-Nachzug braucht nicht-zirkuläre RG-/Gauge-Blockdaten,",
            "einen OS/RP-kompatiblen Observablenkegel und ein unabhängiges Repair- oder",
            "Cofactor-Zertifikat. Ohne diese Felder bleibt das Ledger ein Guardrail gegen",
            "falsche positive Kingman-/Mittelwertsignale.",
            "",
            "## Artefakte",
            "",
            f"- Input: `{payload['input_path']}`",
            f"- Zusammenfassung: `_results/OS_CAPACITY_LEDGER_V2_{payload['date']}.csv`",
            f"- JSON: `_results/OS_CAPACITY_LEDGER_V2_{payload['date']}.json`",
            "",
        ]
    )
    return "\n".join(lines)


def main_legacy() -> None:
    payload = {"runs": [], "synthetic_controls": synthetic_controls()}

    for beta in [2.0, 4.0, 8.0]:
        rows = rg_rows(beta)
        payload["runs"].append({"beta": beta, "rows": rows, "summary": summarize(rows)})

    json_path = RESULT_DIR / "OS_CAPACITY_LEDGER_2026-05-28.json"
    md_path = RESULT_DIR / "OS_CAPACITY_LEDGER_2026-05-28.md"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md_path.write_text(markdown_report(payload), encoding="utf-8")

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Yang-Mills OS-capacity ledgers.")
    parser.add_argument("--legacy", action="store_true", help="write the 2026-05-28 ledger")
    parser.add_argument("--input", type=Path, help="CSV input with v2 ledger fields")
    parser.add_argument("--date-tag", default="2026-06-04")
    parser.add_argument("--rows", type=int, default=500, help="rows per generated control")
    args = parser.parse_args()

    if args.legacy:
        main_legacy()
        return

    if args.input:
        input_path = args.input
        rows = read_input_csv(input_path)
    else:
        rows = generated_v2_control_rows(args.rows)
        input_path = DATA_DIR / f"OS_CAPACITY_LEDGER_V2_CONTROL_INPUT_{args.date_tag}.csv"
        write_input_csv(rows, input_path)

    payload = v2_payload(rows, input_path, args.date_tag)

    json_path = RESULT_DIR / f"OS_CAPACITY_LEDGER_V2_{args.date_tag}.json"
    csv_path = RESULT_DIR / f"OS_CAPACITY_LEDGER_V2_{args.date_tag}.csv"
    md_path = RESULT_DIR / f"OS_CAPACITY_LEDGER_V2_{args.date_tag}.md"

    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_v2_summary_csv(payload["summaries"], csv_path)
    md_path.write_text(markdown_v2_report(payload), encoding="utf-8")

    print(f"Wrote {input_path}")
    print(f"Wrote {json_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
