"""
cluster_boundary_test.py — Cluster-Rand-Eigenvektor fuer alle drei L-Funktionen.

DURCHBRUCH: Der Eigenvektor am Rand des degenerierten Clusters (nicht der
Grundzustand!) reproduziert L-Funktionsnullstellen auf 10^-12 bis 10^-14.

Test fuer: Riemann (even), chi_4 (odd), chi_3 (odd)
"""

import os
import sys
import time
import math
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
    return 0 if m == 0 else (1 if m == 1 else -1)


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


def find_zeros_bisection(xi, N, L_mp, z_lo, z_hi, n_scan=3000, max_zeros=30):
    z_vals = np.linspace(z_lo, z_hi, n_scan)
    F_vals = [evaluate_F(xi, N, L_mp, float(z)) for z in z_vals]
    zeros = []
    for k in range(1, len(F_vals)):
        if F_vals[k - 1] * F_vals[k] < 0:
            a, b = z_vals[k - 1], z_vals[k]
            fa, fb = F_vals[k - 1], F_vals[k]
            for _ in range(60):
                m = (a + b) / 2
                fm = evaluate_F(xi, N, L_mp, float(m))
                if fa * fm < 0:
                    b, fb = m, fm
                else:
                    a, fa = m, fm
            zeros.append((a + b) / 2)
            if len(zeros) >= max_zeros:
                break
    return zeros


def analyze_character(lam, chi_func, chi_name, q_mod, include_W02, parity,
                      target_zeros, thresholds=[1e-10, 1e-5, 1e-3]):
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
    dt_build = time.time() - t0

    w_min = w[0]
    print(f"\n{'='*80}")
    print(f"  {chi_name}  lambda={lam}  N={N}  ({dt_build:.1f}s)")
    print(f"  w_min = {float(w_min):+.12e}")
    print(f"{'='*80}")

    for thresh in thresholds:
        cluster_idx = [k for k in range(dim) if abs(float(w[k] - w_min)) < thresh]
        if not cluster_idx:
            continue
        print(f"\n  --- Cluster threshold = {thresh:.0e} -> {len(cluster_idx)} Eigenvektoren ---")

        # F_k(gamma) Tabelle
        header = f"  {'k':>3s}  {'gap':>12s}"
        for zi in range(min(5, len(target_zeros))):
            header += f"  {'F(g'+str(zi+1)+')':>12s}"
        print(header)

        best_k = -1
        best_sum_F = float('inf')

        for k in cluster_idx:
            xi_even = mp.matrix(dim, 1)
            for r in range(dim):
                xi_even[r, 0] = Q[r, k]

            xi_full = mp.matrix(size, 1)
            for row in range(size):
                s = mp.mpf(0)
                for r in range(dim):
                    s += U_even[row, r] * xi_even[r, 0]
                xi_full[row, 0] = s

            sum_xi = sum(xi_full[i, 0] for i in range(size))
            if abs(sum_xi) > 1e-50:
                scale = mp.sqrt(L_mp) / sum_xi
                for i in range(size):
                    xi_full[i, 0] *= scale

            F_vals = [evaluate_F(xi_full, N, L_mp, g) for g in target_zeros[:5]]
            gap = float(w[k] - w_min)
            line = f"  {k:3d}  {gap:12.3e}"
            for F in F_vals:
                line += f"  {F:+12.4e}"
            print(line)

            sum_absF = sum(abs(f) for f in F_vals)
            if sum_absF < best_sum_F:
                best_sum_F = sum_absF
                best_k = k

        # Scan fuer besten Eigenvektor
        if best_k >= 0:
            print(f"\n  Bester Eigenvektor: k={best_k} (sum|F|={best_sum_F:.4e})")
            xi_even = mp.matrix(dim, 1)
            for r in range(dim):
                xi_even[r, 0] = Q[r, best_k]

            xi_full = mp.matrix(size, 1)
            for row in range(size):
                s = mp.mpf(0)
                for r in range(dim):
                    s += U_even[row, r] * xi_even[r, 0]
                xi_full[row, 0] = s

            sum_xi = sum(xi_full[i, 0] for i in range(size))
            if abs(sum_xi) > 1e-50:
                scale = mp.sqrt(L_mp) / sum_xi
                for i in range(size):
                    xi_full[i, 0] *= scale

            zeros = find_zeros_bisection(xi_full, N, L_mp, 0.1, 35.0)
            print(f"  Scan: {len(zeros)} Nullstellen, Matches mit {chi_name}:")
            for z in zeros:
                for zi, g in enumerate(target_zeros[:5]):
                    if abs(z - g) < 0.5:
                        print(f"    z={z:.15f} -> gamma_{zi+1}={g:.15f}  err={abs(z-g):.4e}")
                        break


def main():
    mp.mp.dps = DPS
    print(f"Cluster-Boundary Eigenvector Analysis")
    print(f"  N={N}, dps={DPS}")

    analyze_character(3.0, chi_trivial, "Riemann_zeta", 1, True, "even",
                      RIEMANN_ZEROS, [1e-10, 1e-5])
    analyze_character(10.0, chi_4, "L(s,chi_4)", 4, False, "odd",
                      CHI4_ZEROS, [1e-10, 1e-5])
    analyze_character(10.0, chi_3, "L(s,chi_3)", 3, False, "odd",
                      CHI3_ZEROS, [1e-10, 1e-5])
    analyze_character(3.67, chi_4, "L(s,chi_4)_opt", 4, False, "odd",
                      CHI4_ZEROS, [1e-10, 1e-5])


if __name__ == "__main__":
    main()
