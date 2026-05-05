"""
multi_eigenvector_analysis.py — Multi-Eigenvector Zero Analysis.

Statt eines einzelnen IEW-Vektors: Jeder Eigenvektor im degenerierten
Grundzustands-Cluster liefert eine eigene F_k(z). Deren Nullstellen
koennen VERSCHIEDENE L-Funktionsnullstellen enkodieren.

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
LAMBDA = float(os.environ.get("LAMBDA", 3.0))
CLUSTER_THRESHOLD = float(os.environ.get("CLUSTER_THRESHOLD", 1e-10))

RIEMANN_ZEROS = [14.134725141734693, 21.022039638771554, 25.010857580145688,
                 30.424876125859513, 32.935061587739189]


def evaluate_F(xi, N, L_mp, z):
    size = 2 * N + 1
    z_mp = mp.mpf(z)
    s = mp.mpf(0)
    for idx in range(size):
        j = idx - N
        pole = 2 * mp.pi * j / L_mp
        s += xi[idx, 0] / (z_mp - pole)
    return float(s)


def find_zeros_bisection(xi, N, L_mp, z_lo, z_hi, n_scan=2000, max_zeros=20):
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


def match_zeros(found_zeros, target_zeros, tol=0.5):
    matches = []
    for z in found_zeros:
        best_k, best_d = -1, tol
        for k, g in enumerate(target_zeros):
            d = abs(z - g)
            if d < best_d:
                best_k, best_d = k, d
        if best_k >= 0:
            matches.append((z, best_k, best_d))
    return matches


def main():
    mp.mp.dps = DPS
    size = 2 * N + 1
    dim = N + 1

    print(f"Multi-Eigenvector Zero Analysis")
    print(f"  lambda={LAMBDA}, N={N}, dps={DPS}")
    print(f"  cluster_threshold={CLUSTER_THRESHOLD}")
    print("=" * 80)

    # Riemann zeta (even, W02)
    M, L_mp = build_QW_mp(
        N, LAMBDA, chi_func=chi_trivial, q_mod=1,
        include_W02=True, sign_WR=-1, parity="even",
        conductor_correction=False, verbose=True,
    )

    M_even, U_even = project_to_parity(M, N, parity="even")
    for i in range(dim):
        for j in range(i + 1, dim):
            avg = (M_even[i, j] + M_even[j, i]) / 2
            M_even[i, j] = avg
            M_even[j, i] = avg

    w, Q = diagonalize_mp(M_even, verbose=True)

    # Degenerate cluster
    w_min = w[0]
    cluster_idx = [k for k in range(dim) if abs(float(w[k] - w_min)) < CLUSTER_THRESHOLD]
    print(f"\n  Degenerate cluster: {len(cluster_idx)} Eigenvektoren (threshold={CLUSTER_THRESHOLD})")
    print(f"  w_min = {float(w_min):+.12e}")
    for k in cluster_idx:
        print(f"    k={k}: w={float(w[k]):+.15e}  gap_from_min={float(w[k]-w_min):.3e}")

    # F_k(z) fuer jeden Eigenvektor im Cluster
    print(f"\n  F_k(gamma) fuer jeden Cluster-Eigenvektor:")
    print(f"  {'k':>3s}  {'F(g1)':>12s}  {'F(g2)':>12s}  {'F(g3)':>12s}  {'F(g4)':>12s}  {'F(g5)':>12s}")

    best_per_zero = {}

    for k in cluster_idx:
        # Embed eigenvector back to full space
        xi_even = mp.matrix(dim, 1)
        for r in range(dim):
            xi_even[r, 0] = Q[r, k]

        xi_full = mp.matrix(size, 1)
        for row in range(size):
            s = mp.mpf(0)
            for r in range(dim):
                s += U_even[row, r] * xi_even[r, 0]
            xi_full[row, 0] = s

        # Normierung: sum(xi_j) = sqrt(L)
        sum_xi = sum(xi_full[i, 0] for i in range(size))
        if abs(sum_xi) > 1e-50:
            scale = mp.sqrt(L_mp) / sum_xi
            for i in range(size):
                xi_full[i, 0] *= scale

        F_vals_at_zeros = []
        for g in RIEMANN_ZEROS[:5]:
            F = evaluate_F(xi_full, N, L_mp, g)
            F_vals_at_zeros.append(F)

        print(f"  {k:3d}  {F_vals_at_zeros[0]:+12.4e}  {F_vals_at_zeros[1]:+12.4e}  "
              f"{F_vals_at_zeros[2]:+12.4e}  {F_vals_at_zeros[3]:+12.4e}  {F_vals_at_zeros[4]:+12.4e}")

        for zi, F in enumerate(F_vals_at_zeros):
            if zi not in best_per_zero or abs(F) < abs(best_per_zero[zi][1]):
                best_per_zero[zi] = (k, F)

    print(f"\n  Bester Eigenvektor pro Nullstelle:")
    for zi in sorted(best_per_zero.keys()):
        k, F = best_per_zero[zi]
        print(f"    gamma_{zi+1} = {RIEMANN_ZEROS[zi]:.10f} -> best k={k}, |F|={abs(F):.4e}")

    # Scan fuer die besten Eigenvektoren
    print(f"\n  Scan-Zeros fuer beste Eigenvektoren:")
    scanned_ks = set()
    for zi in sorted(best_per_zero.keys())[:3]:
        k, _ = best_per_zero[zi]
        if k in scanned_ks:
            continue
        scanned_ks.add(k)

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

        zeros = find_zeros_bisection(xi_full, N, L_mp, 0.1, 35.0)
        matches = match_zeros(zeros, RIEMANN_ZEROS)

        print(f"\n  Eigenvektor k={k}: {len(zeros)} Zeros, {len(matches)} Matches")
        for z, mk, md in matches:
            print(f"    z={z:.12f} -> gamma_{mk+1}={RIEMANN_ZEROS[mk]:.12f}  err={md:.3e}")

    # IEW zum Vergleich
    print(f"\n  --- IEW zum Vergleich ---")
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
    xi_iew_even = mp.matrix(dim, 1)
    for k in range(dim):
        weight = overlaps[k] / (w[k] - w_min + eps)
        for r in range(dim):
            xi_iew_even[r, 0] += weight * Q[r, k]
    norm2 = sum(xi_iew_even[r, 0] ** 2 for r in range(dim))
    if norm2 > 0:
        norm = mp.sqrt(norm2)
        for r in range(dim):
            xi_iew_even[r, 0] /= norm

    xi_iew_full = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for k in range(dim):
            s += U_even[row, k] * xi_iew_even[k, 0]
        xi_iew_full[row, 0] = s

    sum_xi = sum(xi_iew_full[i, 0] for i in range(size))
    if abs(sum_xi) > 1e-50:
        scale = mp.sqrt(L_mp) / sum_xi
        for i in range(size):
            xi_iew_full[i, 0] *= scale

    print(f"  IEW F(gamma_k):")
    for zi, g in enumerate(RIEMANN_ZEROS[:5]):
        F = evaluate_F(xi_iew_full, N, L_mp, g)
        print(f"    gamma_{zi+1} = {g:.12f}   F = {F:+.6e}")

    zeros = find_zeros_bisection(xi_iew_full, N, L_mp, 0.1, 35.0)
    matches = match_zeros(zeros, RIEMANN_ZEROS)
    print(f"  IEW Scan: {len(zeros)} Zeros, {len(matches)} Matches")
    for z, mk, md in matches:
        print(f"    z={z:.12f} -> gamma_{mk+1}={RIEMANN_ZEROS[mk]:.12f}  err={md:.3e}")


if __name__ == "__main__":
    main()
