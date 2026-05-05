"""
optimal_lambda_benchmark.py — Optimale Lambda-Wahl pro L-Funktion.

Berechnet F(gamma_k) fuer alle drei L-Funktionen (Riemann, chi_4, chi_3)
bei den per Kollisionsanalyse optimalen Lambda-Werten.

Optimale Lambda (max min-distance zu ersten 5 Nullstellen):
  Riemann: lambda=2.65  (min_d=0.699, k_max=7)
           lambda=10.27 (min_d=0.558, k_max=105)
  chi_4:   lambda=3.67  (min_d=0.571, k_max=13)
  chi_3:   lambda=2.68  (min_d=1.168, k_max=7)

Autor: LG (Opus 4.7 [1M], Session 19)
"""

import os
import sys
import time
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp,
    chi_trivial,
    chi_4,
    project_to_parity,
    diagonalize_mp,
)

N = int(os.environ.get("N", 30))
DPS = int(os.environ.get("DPS", 100))


def chi_3(n):
    m = n % 3
    if m == 0:
        return 0
    elif m == 1:
        return 1
    else:
        return -1


RIEMANN_ZEROS = [14.134725141734693, 21.022039638771554, 25.010857580145688,
                 30.424876125859513, 32.935061587739189]
CHI4_ZEROS = [6.020948406467813, 10.243766565498344, 12.988096097420339,
              16.343297172788686, 18.291996171925028]
CHI3_ZEROS = [8.039737159558621, 13.922653602824984, 17.821443649486849,
              21.139710893970584, 24.282206057008690]


def evaluate_F(xi, N, L_mp, z):
    size = 2 * N + 1
    z_mp = mp.mpf(z)
    s = mp.mpf(0)
    for idx in range(size):
        j = idx - N
        pole = 2 * mp.pi * j / L_mp
        s += xi[idx, 0] / (z_mp - pole)
    return float(s)


def run_full(lam, chi_func, chi_name, q_mod, include_W02, parity, target_zeros, other_zeros_dict):
    t0 = time.time()
    size = 2 * N + 1
    dim = N + 1

    M, L_mp = build_QW_mp(
        N, lam, chi_func=chi_func, q_mod=q_mod,
        include_W02=include_W02, sign_WR=-1,
        parity=parity, conductor_correction=False, verbose=False,
    )

    M_even, U_even = project_to_parity(M, N, parity="even")
    for i in range(dim):
        for j in range(i + 1, dim):
            avg = (M_even[i, j] + M_even[j, i]) / 2
            M_even[i, j] = avg
            M_even[j, i] = avg

    w, Q = diagonalize_mp(M_even, verbose=False)

    delta_N = mp.matrix(size, 1)
    inv_sqrt_L = 1 / mp.sqrt(L_mp)
    for i in range(size):
        delta_N[i, 0] = inv_sqrt_L

    delta_even = mp.matrix(dim, 1)
    for r in range(dim):
        s = mp.mpf(0)
        for row in range(size):
            s += U_even[row, r] * delta_N[row, 0]
        delta_even[r, 0] = s

    overlaps = []
    for k in range(dim):
        ov = mp.mpf(0)
        for r in range(dim):
            ov += delta_even[r, 0] * Q[r, k]
        overlaps.append(ov)

    eps = mp.mpf("1e-6")
    w_min = w[0]
    xi_even = mp.matrix(dim, 1)
    for k in range(dim):
        weight = overlaps[k] / (w[k] - w_min + eps)
        for r in range(dim):
            xi_even[r, 0] += weight * Q[r, k]
    norm2 = sum(xi_even[r, 0] ** 2 for r in range(dim))
    if norm2 > 0:
        norm = mp.sqrt(norm2)
        for r in range(dim):
            xi_even[r, 0] /= norm

    xi_full = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for k in range(dim):
            s += U_even[row, k] * xi_even[k, 0]
        xi_full[row, 0] = s

    sum_xi = sum(xi_full[i, 0] for i in range(size))
    if abs(sum_xi) > 1e-50:
        scale = mp.sqrt(L_mp) / sum_xi
        for i in range(size):
            xi_full[i, 0] *= scale

    dt = time.time() - t0

    print(f"\n{'='*80}")
    print(f"  {chi_name}  lambda={lam:.2f}  L={float(L_mp):.4f}  N={N}  ({dt:.1f}s)")
    print(f"  w_0 = {float(w[0]):+.8e}   gap(0,1) = {float(w[1]-w[0]):.3e}")
    print(f"{'='*80}")

    import math
    print(f"\n  Target Nullstellen ({chi_name}):")
    for k, g in enumerate(target_zeros[:5]):
        F = evaluate_F(xi_full, N, L_mp, g)
        nearest_j = min(range(-N, N+1), key=lambda j: abs(g - 2*math.pi*j/float(L_mp)))
        d = abs(g - 2*math.pi*nearest_j/float(L_mp))
        print(f"    gamma_{k+1} = {g:.12f}   F = {F:+.6e}   pole_d = {d:.4f} (j={nearest_j})")

    for other_name, other_zeros in other_zeros_dict.items():
        print(f"\n  Sanity: {other_name} Nullstellen (sollten NICHT null sein):")
        for k, g in enumerate(other_zeros[:2]):
            F = evaluate_F(xi_full, N, L_mp, g)
            print(f"    gamma_{k+1} = {g:.10f}   F = {F:+.6e}")

    # Scan
    z_vals = np.linspace(0.1, 35.0, 3500)
    F_vals = [evaluate_F(xi_full, N, L_mp, float(z)) for z in z_vals]
    zeros_found = []
    for k in range(1, len(F_vals)):
        if F_vals[k - 1] * F_vals[k] < 0:
            z_lo, z_hi = z_vals[k - 1], z_vals[k]
            f_lo, f_hi = F_vals[k - 1], F_vals[k]
            for _ in range(60):
                z_mid = (z_lo + z_hi) / 2
                f_mid = evaluate_F(xi_full, N, L_mp, float(z_mid))
                if f_lo * f_mid < 0:
                    z_hi, f_hi = z_mid, f_mid
                else:
                    z_lo, f_lo = z_mid, f_mid
            zeros_found.append((z_lo + z_hi) / 2)
            if len(zeros_found) >= 30:
                break

    print(f"\n  Scan-Nullstellen (Match mit {chi_name}):")
    for i, z in enumerate(zeros_found):
        match = ""
        for k, g in enumerate(target_zeros[:5]):
            if abs(z - g) < 0.5:
                match = f"  <<< gamma_{k+1} err={abs(z-g):.3e}"
                break
        if match or i < 5:
            print(f"    z_{i+1:2d} = {z:.12f}{match}")


def main():
    mp.mp.dps = DPS
    print(f"Optimales Lambda Benchmark: N={N}, dps={DPS}")
    print(f"{'='*80}")

    configs = [
        (3.0,  chi_trivial, "Riemann_zeta", 1, True,  "even",
         RIEMANN_ZEROS, {"chi_4": CHI4_ZEROS, "chi_3": CHI3_ZEROS}),
        (2.65, chi_trivial, "Riemann_zeta", 1, True,  "even",
         RIEMANN_ZEROS, {"chi_4": CHI4_ZEROS}),
        (10.27, chi_trivial, "Riemann_zeta", 1, True,  "even",
         RIEMANN_ZEROS, {"chi_4": CHI4_ZEROS}),
        (3.67, chi_4, "chi_4", 4, False, "odd",
         CHI4_ZEROS, {"Riemann": RIEMANN_ZEROS, "chi_3": CHI3_ZEROS}),
        (10.0, chi_4, "chi_4", 4, False, "odd",
         CHI4_ZEROS, {"Riemann": RIEMANN_ZEROS}),
        (2.68, chi_3, "chi_3", 3, False, "odd",
         CHI3_ZEROS, {"Riemann": RIEMANN_ZEROS, "chi_4": CHI4_ZEROS}),
        (10.0, chi_3, "chi_3", 3, False, "odd",
         CHI3_ZEROS, {"Riemann": RIEMANN_ZEROS}),
    ]

    for lam, chi_func, chi_name, q_mod, w02, parity, target, others in configs:
        run_full(lam, chi_func, chi_name, q_mod, w02, parity, target, others)


if __name__ == "__main__":
    main()
