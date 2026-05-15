"""
prime_weil_bridge_defect.py
===========================

Prototype for the Prime-Weil bridge defect discussed in the Zeta Zoo v2
paper and its public result summaries.

The diagnostic mode computes, on the existing Zookeeper prolate/Galerkin
test space,

    B_lambda^Q(f_n) = Q_Z,lambda(f_n) - Q_arch,lambda(f_n) - Q_prime,lambda(f_n)

where Q_Z is still the zero-sum diagnostic term. This is useful for scale
calibration, but it is not a proof object.

The operator-fourier mode is a separate null-free CCM control. It does not use
zero data and computes the finite Fourier-basis matrix

    QW = W02 - WR - Wp

from the existing Zookeeper CCM builders. It is not in the same basis as the
diagnostic prolate/Galerkin run, so compare signs and scales, not individual
mode entries.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parents[3]
ZOOKEEPER_SCRIPTS = RESEARCH / ".LAB" / ".ZETA-ZOO" / "CORE" / "zookeeper" / "_scripts"
RESULTS_DIR = ROOT / "_results"

sys.path.insert(0, str(ZOOKEEPER_SCRIPTS))

from lambda_scaling_weil_bilinear import (  # noqa: E402
    RIEMANN_ZEROS,
    arch_form,
    build_basis,
    prime_form,
    rho_form_diagonal,
)


DEFAULT_LAMBDAS = [float(np.sqrt(14.0)), float(2.0 * np.sqrt(14.0)), float(4.0 * np.sqrt(14.0))]


def chi_trivial(_: int) -> complex:
    return 1.0 + 0.0j


def summarize(values: np.ndarray) -> dict:
    return {
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
        "l2": float(np.linalg.norm(values)),
        "positive_fraction": float(np.mean(values >= 0.0)),
        "negative_fraction": float(np.mean(values < 0.0)),
    }


def mp_matrix_to_numpy(matrix) -> np.ndarray:
    arr = np.zeros((matrix.rows, matrix.cols), dtype=float)
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            arr[i, j] = float(matrix[i, j])
    return arr


def run_diagnostic(lambda_values: list[float], n_galerkin: int, rho_cutoff: float) -> dict:
    runs = []
    for lam in lambda_values:
        t0 = time.time()
        basis = build_basis(lam, n_galerkin, rho_cutoff)
        arch_mat = arch_form(basis, parity=+1, q=1)
        prime_mat = prime_form(basis, chi=chi_trivial, q=1, lam=lam)
        qz_diag = rho_form_diagonal(basis, RIEMANN_ZEROS, rho_cutoff)
        arch_diag = np.real(np.diag(arch_mat))
        prime_diag = np.real(np.diag(prime_mat))
        bridge = qz_diag - arch_diag - prime_diag
        no_zero_balance = -arch_diag - prime_diag

        runs.append(
            {
                "lambda": float(lam),
                "log_lambda": float(np.log(lam)),
                "n_galerkin": int(n_galerkin),
                "rho_cutoff": float(rho_cutoff),
                "t_pw": float(basis["T_PW"]),
                "t_wide": float(basis["T_wide"]),
                "n_grid": int(basis["n_grid"]),
                "elapsed_seconds": float(time.time() - t0),
                "qz_diag": summarize(qz_diag),
                "arch_diag": summarize(arch_diag),
                "prime_diag": summarize(prime_diag),
                "bridge": summarize(bridge),
                "no_zero_balance": summarize(no_zero_balance),
                "first_modes": [
                    {
                        "n": int(i),
                        "qz": float(qz_diag[i]),
                        "arch": float(arch_diag[i]),
                        "prime": float(prime_diag[i]),
                        "bridge": float(bridge[i]),
                        "no_zero_balance": float(no_zero_balance[i]),
                    }
                    for i in range(min(8, len(bridge)))
                ],
            }
        )

    return {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "mode": "diagnostic",
        "warning": (
            "Q_Z is computed from an explicit zero list. This calibrates the "
            "Prime-Weil defect scale but is not a non-tautological OP5 bridge."
        ),
        "config": {
            "lambda_values": [float(x) for x in lambda_values],
            "n_galerkin": int(n_galerkin),
            "rho_cutoff": float(rho_cutoff),
            "zeros_used": [float(z) for z in RIEMANN_ZEROS if abs(z) <= rho_cutoff],
        },
        "runs": runs,
    }


def run_operator_fourier(lambda_values: list[float], n_fourier: int, dps: int) -> dict:
    import mpmath as mp
    from dirichlet_ccm_fourier_mp import (
        build_W02_matrix_mp,
        build_WR_matrix_mp,
        build_Wprime_matrix_mp,
        chi_trivial as ccm_chi_trivial,
        project_to_parity,
    )

    mp.mp.dps = dps
    runs = []
    for lam in lambda_values:
        t0 = time.time()
        L_mp = 2 * mp.log(mp.mpf(lam))
        size = 2 * n_fourier + 1
        dim = n_fourier + 1

        W02 = build_W02_matrix_mp(n_fourier, L_mp)
        WR = build_WR_matrix_mp(n_fourier, L_mp, parity="even", verbose=False)
        Wp = build_Wprime_matrix_mp(n_fourier, L_mp, ccm_chi_trivial, 1, verbose=False)

        dummy = mp.matrix(size, size)
        _, U_even = project_to_parity(dummy, n_fourier, parity="even")
        W02e = U_even.T * W02 * U_even
        WRe = U_even.T * WR * U_even
        Wpe = U_even.T * Wp * U_even

        W02_np = mp_matrix_to_numpy(W02e)
        WR_np = mp_matrix_to_numpy(WRe)
        Wp_np = mp_matrix_to_numpy(Wpe)
        H_np = -WR_np - Wp_np
        QW_np = W02_np + H_np
        QW_np = 0.5 * (QW_np + QW_np.T)

        w = np.linalg.eigvalsh(QW_np)
        qwe_diag = np.diag(QW_np)
        h_diag = np.diag(H_np)
        w02_diag = np.diag(W02_np)
        wr_diag = np.diag(WR_np)
        wp_diag = np.diag(Wp_np)

        runs.append(
            {
                "lambda": float(lam),
                "log_lambda": float(np.log(lam)),
                "L_ccm": float(L_mp),
                "n_fourier": int(n_fourier),
                "dimension_even": int(dim),
                "dps": int(dps),
                "elapsed_seconds": float(time.time() - t0),
                "w02_diag": summarize(w02_diag),
                "h_diag": summarize(h_diag),
                "qwe_diag": summarize(qwe_diag),
                "eigenvalues": summarize(w),
                "first_modes": [
                    {
                        "n": int(i),
                        "W02": float(w02_diag[i]),
                        "WR": float(wr_diag[i]),
                        "Wp": float(wp_diag[i]),
                        "H": float(h_diag[i]),
                        "QW": float(qwe_diag[i]),
                    }
                    for i in range(min(8, dim))
                ],
                "lowest_eigenvalues": [float(x) for x in w[: min(8, len(w))]],
            }
        )

    return {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "mode": "operator-fourier",
        "warning": (
            "Null-free Fourier-CCM control. This uses QW = W02 - WR - Wp "
            "from existing CCM matrices, not the prolate diagnostic basis."
        ),
        "config": {
            "lambda_values": [float(x) for x in lambda_values],
            "n_fourier": int(n_fourier),
            "dps": int(dps),
            "uses_zero_list": False,
        },
        "runs": runs,
    }


def write_outputs(result: dict) -> tuple[Path, Path]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    stem = "PRIME_WEIL_BRIDGE_DEFECT"
    if result.get("mode") == "operator-fourier":
        stem += "_OPERATOR_FOURIER"
    json_path = RESULTS_DIR / f"{stem}.json"
    md_path = RESULTS_DIR / f"{stem}.md"

    json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")

    if result.get("mode") == "operator-fourier":
        lines = [
            f"# {stem}",
            "",
            f"Created: `{result['created_at']}`",
            "",
            "**Status:** null-free Fourier-CCM operator control. Basis differs from diagnostic prolate/Galerkin mode.",
            "",
            "| lambda | L_ccm | min diag QW | mean diag QW | min eig QW | pos eig% | seconds |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for run in result["runs"]:
            qwe = run["qwe_diag"]
            eig = run["eigenvalues"]
            lines.append(
                "| "
                f"{run['lambda']:.6g} | {run['L_ccm']:.4f} | "
                f"{qwe['min']:.6g} | {qwe['mean']:.6g} | "
                f"{eig['min']:.6g} | {100.0 * eig['positive_fraction']:.1f}% | "
                f"{run['elapsed_seconds']:.2f} |"
            )

        lines.extend(
            [
                "",
                "## Interpretation",
                "",
                "- `QW = W02 - WR - Wp` is built without any zero list.",
                "- This is a Fourier-CCM control, not the same basis as `diagnostic`.",
                "- The useful comparison is qualitative: sign structure, scale, and lambda trend.",
                "",
                "## First Modes",
                "",
            ]
        )
        for run in result["runs"]:
            lines.extend(
                [
                    f"### lambda = {run['lambda']:.6g}",
                    "",
                    "| n | W02 | WR | Wp | H=-WR-Wp | QW |",
                    "|---:|---:|---:|---:|---:|---:|",
                ]
            )
            for row in run["first_modes"]:
                lines.append(
                    "| "
                    f"{row['n']} | {row['W02']:.6g} | {row['WR']:.6g} | "
                    f"{row['Wp']:.6g} | {row['H']:.6g} | {row['QW']:.6g} |"
                )
            lines.append("")

        md_path.write_text("\n".join(lines), encoding="utf-8")
        return json_path, md_path

    lines = [
        "# PRIME_WEIL_BRIDGE_DEFECT",
        "",
        f"Created: `{result['created_at']}`",
        "",
        "**Status:** diagnostic prototype. `Q_Z` still uses an explicit zero list; this is not a proof object.",
        "",
        "| lambda | log(lambda) | min B | mean B | L2 B | pos% | mean no-zero | seconds |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for run in result["runs"]:
        bridge = run["bridge"]
        nz = run["no_zero_balance"]
        lines.append(
            "| "
            f"{run['lambda']:.6g} | {run['log_lambda']:.4f} | "
            f"{bridge['min']:.6g} | {bridge['mean']:.6g} | {bridge['l2']:.6g} | "
            f"{100.0 * bridge['positive_fraction']:.1f}% | "
            f"{nz['mean']:.6g} | {run['elapsed_seconds']:.2f} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- `bridge = Q_Z - Q_arch - Q_prime` is the current quadratic Prime-Weil defect on the prolate/Galerkin basis.",
            "- `no_zero_balance = -Q_arch - Q_prime` shows how much of the balance is supplied before the explicit zero term is inserted.",
            "- The next non-tautological target is an `operator` mode where `Q_Z` is replaced by a CCM/prolate/Weil matrix, not by zeros.",
            "",
            "## First Modes",
            "",
        ]
    )
    for run in result["runs"]:
        lines.extend(
            [
                f"### lambda = {run['lambda']:.6g}",
                "",
                "| n | Q_Z | Q_arch | Q_prime | B | no-zero |",
                "|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for row in run["first_modes"]:
            lines.append(
                "| "
                f"{row['n']} | {row['qz']:.6g} | {row['arch']:.6g} | "
                f"{row['prime']:.6g} | {row['bridge']:.6g} | {row['no_zero_balance']:.6g} |"
            )
        lines.append("")

    md_path.write_text("\n".join(lines), encoding="utf-8")
    return json_path, md_path


def parse_lambdas(raw: str) -> list[float]:
    if raw.strip().lower() == "default":
        return DEFAULT_LAMBDAS
    return [float(part.strip()) for part in raw.split(",") if part.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description="Prime-Weil bridge defect prototype.")
    parser.add_argument("--mode", choices=["diagnostic", "operator", "operator-fourier"], default="diagnostic")
    parser.add_argument("--lambdas", default="default", help='Comma list or "default".')
    parser.add_argument("--n-galerkin", type=int, default=20)
    parser.add_argument("--n-fourier", type=int, default=10)
    parser.add_argument("--dps", type=int, default=int(os.environ.get("DPS", 50)))
    parser.add_argument("--rho-cutoff", type=float, default=40.0)
    args = parser.parse_args()

    lambdas = parse_lambdas(args.lambdas)
    if args.mode in {"operator", "operator-fourier"}:
        result = run_operator_fourier(lambdas, args.n_fourier, args.dps)
    else:
        result = run_diagnostic(lambdas, args.n_galerkin, args.rho_cutoff)

    json_path, md_path = write_outputs(result)
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
