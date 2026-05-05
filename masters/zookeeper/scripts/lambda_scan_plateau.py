"""
lambda_scan_plateau.py — Scanne lambda-Werte um Plateau-Ursache zu finden.

Hypothese: Das Plateau bei ~10^{-8} (Riemann) bzw ~10^{-5} (chi_4) kommt von
der Integral-Abschneidung bei L = 2 log(lambda). Groesseres lambda => groesseres L
=> kleinere Schwanzfehler in alpha_L, beta_L.

Gegenargument: Groesseres lambda => mehr Primen (k_max = lambda^2), was die
prime-Summe veraendert.

Test: F(gamma_1) fuer Riemann und chi_4 bei lambda = 3, 5, 7, 10, 15, 20.
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
LAMBDAS = [3.0, 5.0, 7.0, 10.0, 15.0, 20.0]

G_CHI4 = [6.020948406467813, 10.243766565498344, 12.988096097420339]
G_RIEM = [14.134725141734693, 21.022039638771554, 25.010857580145688]


def evaluate_F(xi, N, L_mp, z):
    size = 2 * N + 1
    z_mp = mp.mpf(z)
    s = mp.mpf(0)
    for idx in range(size):
        j = idx - N
        pole = 2 * mp.pi * j / L_mp
        s += xi[idx, 0] / (z_mp - pole)
    return float(s)


def run_config(lam, chi_func, chi_name, q_mod, include_W02, parity):
    t0 = time.time()
    L_mp = 2 * mp.log(mp.mpf(lam))
    size = 2 * N + 1

    M, L_mp = build_QW_mp(
        N, lam,
        chi_func=chi_func, q_mod=q_mod,
        include_W02=include_W02, sign_WR=-1,
        parity=parity, conductor_correction=False,
        verbose=False,
    )

    M_even, U_even = project_to_parity(M, N, parity="even")
    dim = N + 1
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

    # Inverse-Energy xi
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

    gammas = G_RIEM if chi_name == "trivial" else G_CHI4
    F_vals = [evaluate_F(xi_full, N, L_mp, g) for g in gammas[:3]]

    k_max = int(mp.exp(L_mp))
    print(f"  lam={lam:5.1f}  L={float(L_mp):5.3f}  k_max={k_max:5d}  "
          f"w0={float(w[0]):+.6e}  "
          f"F(g1)={F_vals[0]:+.3e}  F(g2)={F_vals[1]:+.3e}  F(g3)={F_vals[2]:+.3e}  "
          f"({dt:.1f}s)")
    return F_vals


def main():
    mp.mp.dps = DPS
    print(f"Lambda-Scan: N={N}, dps={DPS}")
    print("=" * 100)

    print(f"\n--- Riemann zeta (trivial, even, W02) ---")
    for lam in LAMBDAS:
        run_config(lam, chi_trivial, "trivial", 1, True, "even")

    print(f"\n--- Dirichlet L(s, chi_4) (odd kernel, no W02, no conductor) ---")
    for lam in LAMBDAS:
        run_config(lam, chi_4, "chi_4", 4, False, "odd")


if __name__ == "__main__":
    main()
