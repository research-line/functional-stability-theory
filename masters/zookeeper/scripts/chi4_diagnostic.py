"""
chi4_diagnostic.py — Systematische Diagnose warum chi_4-Nullstellen nicht reproduziert werden.

Teste alle Kombinationen:
  A) Kernel:    even vs odd
  B) Conductor: mit vs ohne
  C) W02:       mit vs ohne
  D) Projektion: even- vs odd-Unterraum
  E) xi-Wahl:   k=0 EV vs Inverse-Energy vs Cluster

Metrik: F(gamma_chi4,1=6.02) und F(gamma_Riemann,1=14.13)
"""

import os
import sys
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_WR_matrix_mp,
    build_W02_matrix_mp,
    build_Wprime_matrix_mp,
    chi_4,
    chi_trivial,
    project_to_parity,
    diagonalize_mp,
)

LAMBDA = 3.0
N = 20
DPS = 80

G_CHI4 = 6.020948406467813
G_RIEM = 14.134725141734693


def evaluate_F(xi, N, L_mp, z):
    size = 2 * N + 1
    z_mp = mp.mpf(z)
    s = mp.mpf(0)
    for idx in range(size):
        j = idx - N
        pole = 2 * mp.pi * j / L_mp
        s += xi[idx, 0] / (z_mp - pole)
    return float(s)


def build_and_test(label, N, lam, chi_func, q_mod, include_W02,
                   parity_kernel, conductor, proj_parity):
    L_mp = 2 * mp.log(mp.mpf(lam))
    size = 2 * N + 1

    WR = build_WR_matrix_mp(N, L_mp, parity=parity_kernel, verbose=False)

    if conductor:
        delta = mp.log(mp.mpf(q_mod) / mp.pi) / L_mp
        for i in range(size):
            for j in range(size):
                WR[i, j] += delta

    Wp = build_Wprime_matrix_mp(N, L_mp, chi_func, q_mod, verbose=False)

    M = mp.matrix(size, size)
    for i in range(size):
        for j in range(size):
            M[i, j] = -WR[i, j] - Wp[i, j]

    if include_W02:
        W02 = build_W02_matrix_mp(N, L_mp)
        for i in range(size):
            for j in range(size):
                M[i, j] += W02[i, j]

    for i in range(size):
        for j in range(i + 1, size):
            avg = (M[i, j] + M[j, i]) / 2
            M[i, j] = avg
            M[j, i] = avg

    M_sub, U_sub = project_to_parity(M, N, parity=proj_parity)
    dim = M_sub.rows

    for i in range(dim):
        for j in range(i + 1, dim):
            avg = (M_sub[i, j] + M_sub[j, i]) / 2
            M_sub[i, j] = avg
            M_sub[j, i] = avg

    w, Q = diagonalize_mp(M_sub, verbose=False)

    delta_N = mp.matrix(size, 1)
    inv_sqrt_L = 1 / mp.sqrt(L_mp)
    for i in range(size):
        delta_N[i, 0] = inv_sqrt_L

    delta_sub = mp.matrix(dim, 1)
    for r in range(dim):
        s = mp.mpf(0)
        for row in range(size):
            s += U_sub[row, r] * delta_N[row, 0]
        delta_sub[r, 0] = s

    overlaps = []
    for k in range(dim):
        ov = mp.mpf(0)
        for r in range(dim):
            ov += delta_sub[r, 0] * Q[r, k]
        overlaps.append(ov)

    results = {}

    # xi = k=0 eigenvector
    xi_sub = mp.matrix(dim, 1)
    for r in range(dim):
        xi_sub[r, 0] = Q[r, 0]
    xi_full = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for k in range(dim):
            s += U_sub[row, k] * xi_sub[k, 0]
        xi_full[row, 0] = s
    sum_xi = sum(xi_full[i, 0] for i in range(size))
    if abs(sum_xi) > 1e-50:
        scale = mp.sqrt(L_mp) / sum_xi
        for i in range(size):
            xi_full[i, 0] *= scale
    results['k0_chi4'] = evaluate_F(xi_full, N, L_mp, G_CHI4)
    results['k0_riem'] = evaluate_F(xi_full, N, L_mp, G_RIEM)

    # xi = Inverse-Energy
    eps = mp.mpf("1e-6")
    w_min = w[0]
    xi_sub = mp.matrix(dim, 1)
    for k in range(dim):
        weight = overlaps[k] / (w[k] - w_min + eps)
        for r in range(dim):
            xi_sub[r, 0] += weight * Q[r, k]
    norm2 = sum(xi_sub[r, 0] ** 2 for r in range(dim))
    if norm2 > 0:
        norm = mp.sqrt(norm2)
        for r in range(dim):
            xi_sub[r, 0] /= norm
    xi_full = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for k in range(dim):
            s += U_sub[row, k] * xi_sub[k, 0]
        xi_full[row, 0] = s
    sum_xi = sum(xi_full[i, 0] for i in range(size))
    if abs(sum_xi) > 1e-50:
        scale = mp.sqrt(L_mp) / sum_xi
        for i in range(size):
            xi_full[i, 0] *= scale
    results['ie_chi4'] = evaluate_F(xi_full, N, L_mp, G_CHI4)
    results['ie_riem'] = evaluate_F(xi_full, N, L_mp, G_RIEM)

    print(f"  {label:45s}  "
          f"F(g_c4)={results['ie_chi4']:+.3e}  "
          f"F(g_R)={results['ie_riem']:+.3e}  "
          f"w0={float(w[0]):+.6e}  "
          f"|ov0|={float(abs(overlaps[0])):.3e}")


def main():
    mp.mp.dps = DPS
    print(f"chi4 Diagnostik: N={N}, lambda={LAMBDA}, dps={DPS}")
    print(f"  G_chi4={G_CHI4:.6f}, G_Riem={G_RIEM:.6f}")
    print("=" * 100)

    # Referenz: Riemann (trivial, even kernel, W02, even proj)
    print("\n--- Referenz: Riemann (trivial character) ---")
    build_and_test("triv/even-kern/W02/even-proj/conductor=no",
                   N, LAMBDA, chi_trivial, 1, True, "even", False, "even")

    # chi_4 Varianten
    print("\n--- chi_4 Varianten ---")

    configs = [
        ("chi4/odd-kern/noW02/even-proj/cond=yes",  "odd",  False, True,  "even"),
        ("chi4/odd-kern/noW02/even-proj/cond=no",   "odd",  False, False, "even"),
        ("chi4/odd-kern/noW02/odd-proj/cond=yes",   "odd",  False, True,  "odd"),
        ("chi4/odd-kern/noW02/odd-proj/cond=no",    "odd",  False, False, "odd"),
        ("chi4/even-kern/noW02/even-proj/cond=yes",  "even", False, True,  "even"),
        ("chi4/even-kern/noW02/even-proj/cond=no",  "even", False, False, "even"),
        ("chi4/even-kern/noW02/odd-proj/cond=yes",  "even", False, True,  "odd"),
        ("chi4/even-kern/noW02/odd-proj/cond=no",   "even", False, False, "odd"),
        ("chi4/odd-kern/W02/even-proj/cond=yes",    "odd",  True,  True,  "even"),
        ("chi4/odd-kern/W02/even-proj/cond=no",     "odd",  True,  False, "even"),
        ("chi4/even-kern/W02/even-proj/cond=no",    "even", True,  False, "even"),
    ]

    for label, kern, w02, cond, proj in configs:
        build_and_test(label, N, LAMBDA, chi_4, 4, w02, kern, cond, proj)


if __name__ == "__main__":
    main()
