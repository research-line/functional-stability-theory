"""
chi3_test.py — Dirichlet chi_3 (mod 3, EVEN character) Validierung.

chi_3(1) = 1, chi_3(2) = -1, chi_3(0) = 0 (Legendre-Symbol mod 3)
chi_3(-1) = chi_3(2) = -1 => chi_3 ist ein EVEN character? Nein:
chi_3(-1) = chi_3(2) = -1 => ODD!

Korrektur: chi_3 mod 3 mit chi_3(2) = -1 hat chi_3(-1) = chi_3(2) = -1 => ODD.
Also: parity="odd", kein W02, ODD-Kernel.

Erste Nullstellen von L(s, chi_3) (LMFDB):
  gamma_1 = 8.039...
  gamma_2 = 13.922...
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
    project_to_parity,
    diagonalize_mp,
)

LAMBDA = float(os.environ.get("LAMBDA", 10.0))
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


# L(s, chi_3) Nullstellen (LMFDB: Character 3.2)
CHI3_ZEROS = [
    8.03973715955862085712590,
    13.92265360282498372950510,
    17.82144364948684860213940,
    21.13971089397058245478250,
    24.28220605700869128199300,
]

RIEMANN_ZEROS = [14.134725141734693, 21.022039638771554]
CHI4_ZEROS = [6.020948406467813, 10.243766565498344]


def evaluate_F(xi, N, L_mp, z):
    size = 2 * N + 1
    z_mp = mp.mpf(z)
    s = mp.mpf(0)
    for idx in range(size):
        j = idx - N
        pole = 2 * mp.pi * j / L_mp
        s += xi[idx, 0] / (z_mp - pole)
    return float(s)


def main():
    mp.mp.dps = DPS
    print("=" * 70)
    print(f"Dirichlet chi_3 CCM-Test")
    print(f"  lambda = {LAMBDA}, N = {N}, dps = {DPS}")
    print(f"  parity = odd  (chi_3(-1) = -1)")
    print(f"  q_mod = 3, kein W02")
    print("=" * 70)

    M, L_mp = build_QW_mp(
        N, LAMBDA,
        chi_func=chi_3,
        q_mod=3,
        include_W02=False,
        sign_WR=-1,
        parity="odd",
        conductor_correction=False,
        verbose=True,
    )

    size = 2 * N + 1
    M_even, U_even = project_to_parity(M, N, parity="even")
    dim = N + 1
    for i in range(dim):
        for j in range(i + 1, dim):
            avg = (M_even[i, j] + M_even[j, i]) / 2
            M_even[i, j] = avg
            M_even[j, i] = avg

    w, Q = diagonalize_mp(M_even, verbose=True)

    print(f"\n  Niedrigste 10 Eigenwerte:")
    for k in range(min(10, dim)):
        gap = float(w[k + 1] - w[k]) if k + 1 < dim else 0
        print(f"    k={k:3d}: w = {float(w[k]):+.12e}   gap = {gap:+.3e}")

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

    print(f"\n  F(gamma) bei chi_3-Nullstellen:")
    for k, g in enumerate(CHI3_ZEROS):
        F = evaluate_F(xi_full, N, L_mp, g)
        print(f"    gamma_chi3,{k+1} = {g:.15f}   F = {F:+.6e}")

    print(f"\n  F(gamma) bei Riemann-Nullstellen (Sanity — sollten NICHT null sein):")
    for k, g in enumerate(RIEMANN_ZEROS):
        F = evaluate_F(xi_full, N, L_mp, g)
        print(f"    gamma_R,{k+1} = {g:.10f}   F = {F:+.6e}")

    print(f"\n  F(gamma) bei chi_4-Nullstellen (Sanity — sollten NICHT null sein):")
    for k, g in enumerate(CHI4_ZEROS):
        F = evaluate_F(xi_full, N, L_mp, g)
        print(f"    gamma_chi4,{k+1} = {g:.10f}   F = {F:+.6e}")

    # Scan
    print(f"\n  Scan F(z) auf [0.1, 30] nach Nullstellen:")
    z_vals = np.linspace(0.1, 30.0, 3000)
    F_vals = [evaluate_F(xi_full, N, L_mp, float(z)) for z in z_vals]
    zeros_found = []
    for k in range(1, len(F_vals) - 1):
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
            if len(zeros_found) >= 20:
                break

    # Match mit chi_3-Nullstellen
    print(f"  Gefundene Nullstellen (Match mit chi_3):")
    for i, z in enumerate(zeros_found[:20]):
        match = ""
        for k, g in enumerate(CHI3_ZEROS):
            if abs(z - g) < 0.5:
                match = f"  *** chi3,{k+1} err={abs(z-g):.3e}"
                break
        print(f"    z_{i+1:2d} = {z:.12f}{match}")


if __name__ == "__main__":
    main()
