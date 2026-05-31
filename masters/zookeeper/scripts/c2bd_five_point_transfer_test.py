"""
c2bd_five_point_transfer_test.py

Five-point transfer test for the derivative-freezing route.

It reuses the established numerical pipeline from the FST `_scripts`
folder and reports the targeted quantities for the current active
question:

1. Is there a stable cubic boundary-transfer factor
      K_lambda = - d/dL(T2_bd) / (lambda^3 * partial_L m_H) ?
2. Does the real derivative locking
      d/dL(T2_bd) / d/dL(R_bulk) ~ -1
   persist on five lambda values?
3. How small is the induced residual
      E_lambda = d/dL(R_bulk) + Kbar * lambda^3 * partial_L m_H
   with Kbar the mean transfer factor?
"""

from __future__ import annotations

import math
import statistics
import sys
from pathlib import Path

import numpy as np


SCRIPT_ROOT = Path(
    r"C:\Users\User\OneDrive\.TOPICS\.RESEARCH\.PRIO-1\DRAFT__META_RH_TREE"
    r"\02_FST_MATHEMATICS\fst_spectrum_duality\_scripts"
)
sys.path.insert(0, str(SCRIPT_ROOT))

from c2bd_derivative_matching import run  # noqa: E402


LAMBDAS = [3.0, 4.0, 5.0, 7.0, 9.0]
N_BY_LAMBDA = {
    3.0: 15,
    4.0: 15,
    5.0: 15,
    7.0: 15,
    9.0: 15,
}
DLAM = 0.01


def cv(values: list[float]) -> float:
    mean = statistics.fmean(values)
    return statistics.pstdev(values) / abs(mean)


def fit_power(xs: list[float], ys: list[float]) -> tuple[float, float]:
    logx = np.log(np.array(xs, dtype=float))
    logy = np.log(np.array(ys, dtype=float))
    p, logc = np.polyfit(logx, logy, 1)
    return float(p), float(math.exp(logc))


def main() -> None:
    print("=" * 96)
    print("  C2bd FIVE-POINT TRANSFER TEST")
    print("=" * 96)
    print(f"  lambdas = {LAMBDAS}")
    print(f"  dlam    = {DLAM}")
    print(f"  N-map   = {N_BY_LAMBDA}")

    results = []
    for lam in LAMBDAS:
        print(f"\nRunning lambda={lam:.1f} ...", flush=True)
        result = run(lam, N_BY_LAMBDA[lam], dlam=DLAM)
        results.append(result)

    k_values = []
    minus_one_devs = []
    freezing_ratios = []
    value_ratios = []
    scaled_residuals = []

    for r in results:
        lam = r["lam"]
        dL_T2bd = r["dL_T2bd"]
        dL_Rbulk = r["dL_Rbulk"]
        dL_mH = r["dL_mH"]
        dL_mua = r["dL_mua"]
        ct_rhs = r["Ct_prime_Ct2"]
        t2_bd = r["T2_bd"]
        r_bulk = r["R_bulk"]

        k_lambda = -dL_T2bd / (lam**3 * dL_mH)
        lock_ratio = dL_T2bd / dL_Rbulk
        freeze_ratio = dL_mua / ct_rhs
        value_ratio = t2_bd / r_bulk if abs(r_bulk) > 1e-30 else float("nan")

        k_values.append(k_lambda)
        minus_one_devs.append(lock_ratio + 1.0)
        freezing_ratios.append(freeze_ratio)
        value_ratios.append(value_ratio)

    k_bar = statistics.fmean(k_values)

    for r in results:
        lam = r["lam"]
        residual = r["dL_Rbulk"] + k_bar * lam**3 * r["dL_mH"]
        scaled_residuals.append(residual)

    raw_ratio = [
        abs(r["dL_T2bd"] / r["dL_mH"])
        for r in results
    ]
    p_fit, c_fit = fit_power(LAMBDAS, raw_ratio)

    print("\n" + "=" * 96)
    print("  FIVE-POINT SUMMARY")
    print("=" * 96)
    header = (
        f"{'lam':>5} {'K_lambda':>14} {'dTbd/dR':>14} {'dev(-1)':>12} "
        f"{'d(mu/a)/rhs':>14} {'T2bd/Rbulk':>14} {'E_lambda':>14}"
    )
    print(header)
    print("-" * len(header))
    for r, k_lambda, dev, freeze_ratio, value_ratio, residual in zip(
        results, k_values, minus_one_devs, freezing_ratios, value_ratios, scaled_residuals
    ):
        lock_ratio = r["dL_T2bd"] / r["dL_Rbulk"]
        print(
            f"{r['lam']:5.1f} {k_lambda:14.6e} {lock_ratio:14.6f} {dev:12.4e} "
            f"{freeze_ratio:14.6e} {value_ratio:14.6f} {residual:14.6e}"
        )

    print("\n" + "=" * 96)
    print("  AGGREGATES")
    print("=" * 96)
    print(f"mean(K_lambda)                 = {k_bar:.8e}")
    print(f"CV(K_lambda)                   = {cv(k_values):.4%}")
    print(f"max |dTbd/dR + 1|              = {max(abs(x) for x in minus_one_devs):.6e}")
    print(f"mean d(mu/a)/rhs               = {statistics.fmean(freezing_ratios):+.8e}")
    print(f"max |E_lambda|                 = {max(abs(x) for x in scaled_residuals):.6e}")
    print(f"log-log fit |dTbd/dmH| ~ C lam^p: p = {p_fit:.6f}, C = {c_fit:.8e}")

    print("\n" + "=" * 96)
    print("  POWER-CANDIDATE CHECKS")
    print("=" * 96)
    print(f"{'power':>6} {'mean':>14} {'CV':>12}")
    for p in range(0, 6):
        vals = [
            -r["dL_T2bd"] / (r["lam"] ** p * r["dL_mH"])
            for r in results
        ]
        print(f"{p:6d} {statistics.fmean(vals):14.6e} {cv(vals):12.4%}")


if __name__ == "__main__":
    main()
