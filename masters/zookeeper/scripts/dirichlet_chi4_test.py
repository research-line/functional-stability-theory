"""
dirichlet_chi4_test.py
======================

Dirichlet chi_4-Test via CCM-Konstruktion (Fourier-Basis, full mpmath).

Erwartete erste Nullstelle von L(s, chi_4) auf kritischer Linie:
    gamma_chi4,1 = 6.020948406467813248841013122025116830842932188530...

Konfiguration (basierend auf Riemann-Sanity-Ergebnissen):
  - lambda = 3.0 (CCM's Test-Beispiel, beste Konvergenz)
  - N = 30 (ausreichend bei 10^-8 Plateau)
  - dps = 100 (ausreichend)
  - parity = 'odd' (chi_4(-1) = -1)
  - include_W02 = False (L-Funktion hat keinen Pol)
  - Conductor: log(4/pi)/L uniform auf WR
  - Inverse-Energy-Weighting fuer xi

Autor: LG (Opus 4.7 [1M], Session 18 — P1)
"""

import os
import sys
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp,
    chi_4,
    project_to_parity,
    diagonalize_mp,
)

LAMBDA = float(os.environ.get("LAMBDA", 3.0))
N = int(os.environ.get("N", 30))
DPS = int(os.environ.get("DPS", 100))

# Erste L(s, chi_4) Nullstellen auf kritischer Linie (LMFDB)
CHI4_ZEROS = [
    6.020948406467813248841013122025,
    10.243766565498344019708550001906,
    12.988096097420338478712946059894,
    16.343297172788686238881030410924,
    18.291996171925028354003894030405,
]


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
    print(f"Dirichlet chi_4 CCM-Test")
    print(f"  lambda = {LAMBDA}")
    print(f"  N      = {N}")
    print(f"  dps    = {DPS}")
    print(f"  parity = odd")
    print(f"  q_mod  = 4")
    print("=" * 70)

    # Baue Matrix fuer chi_4 (ODD character)
    # conductor_correction=False: Diagnostik zeigt Doppelzaehlung mit c_plus_w
    M, L_mp = build_QW_mp(
        N, LAMBDA,
        chi_func=chi_4,
        q_mod=4,
        include_W02=False,
        sign_WR=-1,
        parity="odd",
        conductor_correction=False,
        verbose=True,
    )

    size = 2 * N + 1

    # Even-Projektion
    # Hinweis: Fuer ODD character chi_4 ist der Grundzustand trotzdem
    # in der symmetrisierten Matrix zu suchen. Die Paritaet beschreibt die
    # arch-Kern-Wahl, nicht die Symmetrie des xi-Vektors im Fourier-Raum.
    print(f"\n  Projektion auf EVEN-Unterraum (Dim N+1 = {N+1})...")
    M_even, U_even = project_to_parity(M, N, parity="even")
    for i in range(N + 1):
        for j in range(i + 1, N + 1):
            avg = (M_even[i, j] + M_even[j, i]) / 2
            M_even[i, j] = avg
            M_even[j, i] = avg

    print(f"  Diagonalisiere EVEN-Matrix ({N+1}x{N+1})...")
    w_even, Q_even = diagonalize_mp(M_even, verbose=True)

    print(f"\n  Niedrigste 15 Eigenwerte:")
    for k in range(min(15, len(w_even))):
        if k + 1 < len(w_even):
            gap = float(w_even[k + 1] - w_even[k])
            print(f"    k={k:3d}: w = {float(w_even[k]):+.12e}   gap = {gap:+.3e}")
        else:
            print(f"    k={k:3d}: w = {float(w_even[k]):+.12e}")

    # delta_N im Even-Basis
    delta_N = mp.matrix(size, 1)
    inv_sqrt_L = 1 / mp.sqrt(L_mp)
    for i in range(size):
        delta_N[i, 0] = inv_sqrt_L
    delta_N_even = mp.matrix(N + 1, 1)
    for r in range(N + 1):
        s = mp.mpf(0)
        for row in range(size):
            s += U_even[row, r] * delta_N[row, 0]
        delta_N_even[r, 0] = s

    # Overlaps
    overlaps = []
    for k in range(N + 1):
        ov = mp.mpf(0)
        for r in range(N + 1):
            ov += delta_N_even[r, 0] * Q_even[r, k]
        overlaps.append(ov)

    print(f"\n  Overlaps |<delta_N, Q_even[:,k]>| fuer k=0..10:")
    for k in range(min(11, N + 1)):
        print(f"    k={k:3d}: w = {float(w_even[k]):+.12e}   |ov| = {float(abs(overlaps[k])):.6e}")

    # Inverse-Energy-Weighting als xi
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

    # Einbetten
    xi_full = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for k_inner in range(N + 1):
            s += U_even[row, k_inner] * xi_even[k_inner, 0]
        xi_full[row, 0] = s

    # Normierung: sum(xi_j) = sqrt(L)
    sum_xi = mp.mpf(0)
    for i in range(size):
        sum_xi += xi_full[i, 0]
    print(f"\n  xi (Inverse-Energy): sum = {float(sum_xi):+.6e}   "
          f"sqrt(L) = {float(mp.sqrt(L_mp)):.6f}")
    if abs(sum_xi) > 1e-50:
        scale = mp.sqrt(L_mp) / sum_xi
        for i in range(size):
            xi_full[i, 0] *= scale

    # F(gamma) direkt auswerten fuer chi_4-Nullstellen
    print(f"\n  F(gamma) bei erwarteten chi_4-Nullstellen:")
    for k, gamma_k in enumerate(CHI4_ZEROS):
        F_val = evaluate_F(xi_full, N, L_mp, gamma_k)
        print(f"    k={k+1}: gamma_chi4,{k+1} = {gamma_k:.15f}, F = {F_val:+.6e}")

    # Zum Vergleich: Riemann-gamma_1 (sollte NICHT Nullstelle sein hier)
    print(f"\n  F(gamma) bei Riemann-gamma_k (Sanity — sollten NICHT null sein):")
    riemann = [14.134725141734693, 21.022039638771554]
    for k, gamma_k in enumerate(riemann):
        F_val = evaluate_F(xi_full, N, L_mp, gamma_k)
        print(f"    gamma_Riemann,{k+1} = {gamma_k:.10f}, F = {F_val:+.6e}")

    # Scan nach Vorzeichenwechseln
    print(f"\n  Scan F(z) auf [0.1, 30] nach Nullstellen:")
    z_vals = np.linspace(0.1, 30.0, 2000)
    F_vals = [evaluate_F(xi_full, N, L_mp, float(z)) for z in z_vals]
    zeros_found = []
    for k in range(1, len(F_vals) - 1):
        if F_vals[k - 1] * F_vals[k] < 0:
            # Bisection
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
            if len(zeros_found) >= 15:
                break

    print(f"  Gefundene Nullstellen:")
    for i, z in enumerate(zeros_found[:15]):
        if i < len(CHI4_ZEROS):
            err = abs(z - CHI4_ZEROS[i])
            print(f"    z_{i+1} = {z:.12f}   erwartet = {CHI4_ZEROS[i]:.12f}   err = {err:.3e}")
        else:
            print(f"    z_{i+1} = {z:.12f}")


if __name__ == "__main__":
    main()
