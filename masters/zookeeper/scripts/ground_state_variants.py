"""
ground_state_variants.py
========================

Testet verschiedene Grundzustand-Definitionen auf der Fourier-Matrix QW_lambda
(in Even-Unterraum), um zu sehen, welche am besten die Riemann-Zeros als
F(z)-Nullstellen reproduziert.

Varianten:
  (a) xi = Q_even[:, k]  einzeln fuer k = 0..10 (jeder EV als Kandidat)
  (b) xi = Cluster-Projektion (wie Standard) — Baseline
  (c) xi = sum_k Q_even[:,k] * <delta_N, Q_even[:,k]> / (w_k - w_min + eps)
      (Inverse-Energy-Weighting)
  (d) xi = EW-Vektor mit max |<delta_N, EV>|

Autor: LG (Opus 4.7, Session 18)
"""

import os
import sys
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp,
    chi_trivial,
    project_to_parity,
    diagonalize_mp,
    RIEMANN_ZEROS,
)

LAMBDA = float(os.environ.get("LAMBDA", float(np.sqrt(14.0))))
N = int(os.environ.get("N", 20))
DPS = int(os.environ.get("DPS", 80))


def evaluate_F(xi_full, N, L_mp, z):
    size = 2 * N + 1
    z_mp = mp.mpf(z)
    s = mp.mpf(0)
    for idx in range(size):
        j = idx - N
        pole = 2 * mp.pi * j / L_mp
        s += xi_full[idx, 0] / (z_mp - pole)
    return float(s)


def overlap_with_delta(xi_full, delta_N, size):
    ov = mp.mpf(0)
    for i in range(size):
        ov += delta_N[i, 0] * xi_full[i, 0]
    return float(abs(ov))


def embed_even_to_full(xi_even, U_even, N, size):
    xi_full = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for k in range(N + 1):
            s += U_even[row, k] * xi_even[k, 0]
        xi_full[row, 0] = s
    return xi_full


def normalize_sum(xi_full, size, L_mp):
    sum_xi = mp.mpf(0)
    for i in range(size):
        sum_xi += xi_full[i, 0]
    if abs(sum_xi) > 1e-80:
        scale = mp.sqrt(L_mp) / sum_xi
        xi_norm = mp.matrix(size, 1)
        for i in range(size):
            xi_norm[i, 0] = xi_full[i, 0] * scale
        return xi_norm
    return xi_full


def main():
    mp.mp.dps = DPS
    print(f"Grundzustand-Varianten bei N={N}, dps={DPS}, lambda=sqrt(14)")
    print("=" * 70)

    # Baue Matrix einmal
    M, L_mp = build_QW_mp(
        N, LAMBDA, chi_func=chi_trivial, q_mod=1,
        include_W02=True, sign_WR=-1, verbose=True
    )

    # Even-Unterraum
    M_even, U_even = project_to_parity(M, N, parity="even")
    for i in range(N + 1):
        for j in range(i + 1, N + 1):
            avg = (M_even[i, j] + M_even[j, i]) / 2
            M_even[i, j] = avg
            M_even[j, i] = avg

    w_even, Q_even = diagonalize_mp(M_even, verbose=True)

    # delta_N: 1/sqrt(L) in allen Eintraegen
    size = 2 * N + 1
    delta_N = mp.matrix(size, 1)
    inv_sqrt_L = 1 / mp.sqrt(L_mp)
    for i in range(size):
        delta_N[i, 0] = inv_sqrt_L

    # delta_N im Even-Basis
    delta_N_even = mp.matrix(N + 1, 1)
    for r in range(N + 1):
        s = mp.mpf(0)
        for row in range(size):
            s += U_even[row, r] * delta_N[row, 0]
        delta_N_even[r, 0] = s

    # Overlaps <delta_N_even, Q_even[:,k]> berechnen
    overlaps = []
    for k in range(N + 1):
        ov = mp.mpf(0)
        for r in range(N + 1):
            ov += delta_N_even[r, 0] * Q_even[r, k]
        overlaps.append(ov)

    print(f"\n  Overlaps |<delta_N_even, Q_even[:,k]>| fuer k=0..15:")
    for k in range(min(16, N + 1)):
        print(f"    k={k:3d}: w = {float(w_even[k]):+.12e}   "
              f"|ov| = {float(abs(overlaps[k])):.6e}")

    print(f"\n{'=' * 70}\nVariante (a): Einzelne Eigenvektoren als xi")
    print("=" * 70)
    for k in range(min(15, N + 1)):
        # xi_even = Q_even[:, k], als column
        xi_even = mp.matrix(N + 1, 1)
        for r in range(N + 1):
            xi_even[r, 0] = Q_even[r, k]
        xi_full = embed_even_to_full(xi_even, U_even, N, size)
        ov = overlap_with_delta(xi_full, delta_N, size)
        if ov < 1e-60:
            continue
        xi_n = normalize_sum(xi_full, size, L_mp)
        F1 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[0])
        F2 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[1])
        F3 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[2])
        print(f"  k={k:3d}: w={float(w_even[k]):+.8e}  |ov|={ov:.3e}  "
              f"F(g1)={F1:+.3e}  F(g2)={F2:+.3e}  F(g3)={F3:+.3e}")

    print(f"\n{'=' * 70}\nVariante (b): Cluster-Projektion (Baseline)")
    print("=" * 70)
    cluster_tol = 1e-3
    degen_idx = [k for k in range(N + 1) if abs(float(w_even[k] - w_even[0])) < cluster_tol]
    xi_even = mp.matrix(N + 1, 1)
    for k in degen_idx:
        for r in range(N + 1):
            xi_even[r, 0] += overlaps[k] * Q_even[r, k]
    norm2 = sum(xi_even[r, 0] ** 2 for r in range(N + 1))
    if norm2 > 0:
        norm = mp.sqrt(norm2)
        for r in range(N + 1):
            xi_even[r, 0] /= norm
    xi_full = embed_even_to_full(xi_even, U_even, N, size)
    xi_n = normalize_sum(xi_full, size, L_mp)
    ov = overlap_with_delta(xi_full, delta_N, size)
    F1 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[0])
    F2 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[1])
    F3 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[2])
    print(f"  Cluster k=0..{max(degen_idx)}: |ov|={ov:.3e}  "
          f"F(g1)={F1:+.3e}  F(g2)={F2:+.3e}  F(g3)={F3:+.3e}")

    print(f"\n{'=' * 70}\nVariante (c): Inverse-Energy-Weighting")
    print("=" * 70)
    eps = mp.mpf("1e-6")
    w_min = w_even[0]
    xi_even = mp.matrix(N + 1, 1)
    for k in range(N + 1):
        weight = overlaps[k] / (w_even[k] - w_min + eps)
        for r in range(N + 1):
            xi_even[r, 0] += weight * Q_even[r, k]
    norm2 = sum(xi_even[r, 0] ** 2 for r in range(N + 1))
    if norm2 > 0:
        norm = mp.sqrt(norm2)
        for r in range(N + 1):
            xi_even[r, 0] /= norm
    xi_full = embed_even_to_full(xi_even, U_even, N, size)
    xi_n = normalize_sum(xi_full, size, L_mp)
    ov = overlap_with_delta(xi_full, delta_N, size)
    F1 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[0])
    F2 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[1])
    F3 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[2])
    print(f"  Inv-Energy: |ov|={ov:.3e}  "
          f"F(g1)={F1:+.3e}  F(g2)={F2:+.3e}  F(g3)={F3:+.3e}")

    print(f"\n{'=' * 70}\nVariante (d): Max-Overlap-EV allein")
    print("=" * 70)
    k_max = max(range(N + 1), key=lambda k: float(abs(overlaps[k])))
    xi_even = mp.matrix(N + 1, 1)
    for r in range(N + 1):
        xi_even[r, 0] = Q_even[r, k_max]
    xi_full = embed_even_to_full(xi_even, U_even, N, size)
    xi_n = normalize_sum(xi_full, size, L_mp)
    ov = overlap_with_delta(xi_full, delta_N, size)
    F1 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[0])
    F2 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[1])
    F3 = evaluate_F(xi_n, N, L_mp, RIEMANN_ZEROS[2])
    print(f"  Max-ov k={k_max}: w={float(w_even[k_max]):+.8e}  |ov|={ov:.3e}  "
          f"F(g1)={F1:+.3e}  F(g2)={F2:+.3e}  F(g3)={F3:+.3e}")


if __name__ == "__main__":
    main()
