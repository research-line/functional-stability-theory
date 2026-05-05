"""
cluster_mechanism_analysis.py — Mathematische Analyse des Cluster-Rand-Mechanismus.

Zentrale Frage: WARUM kodiert der Cluster-Rand-Eigenvektor die Nullstellen?

Zerlegung:
  M_even = W02_even + H_even,  wo H = -WR - W'
  W02_even = C * u_tilde @ u_tilde^T  (Rang 1 im Even-Unterraum!)
  P = I - u_tilde @ u_tilde^T  (Nullraum-Projektor von W02_even)
  PHP = Einschraenkung von H auf den Nullraum von W02

Hypothese: Die Cluster-Eigenwerte sind die Eigenwerte von PHP,
und der Rand-Eigenvektor ist maximal an H gekoppelt.
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
    chi_trivial,
    project_to_parity,
    diagonalize_mp,
)

N = int(os.environ.get("N", 30))
DPS = int(os.environ.get("DPS", 100))
LAM = float(os.environ.get("LAMBDA", 3.0))

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


def main():
    mp.mp.dps = DPS
    L_mp = 2 * mp.log(mp.mpf(LAM))
    size = 2 * N + 1
    dim = N + 1

    print(f"Cluster-Mechanismus-Analyse: lambda={LAM}, N={N}, dps={DPS}")
    print(f"L = {float(L_mp):.6f}")
    print("=" * 80)

    # Schritt 1: Baue Einzelkomponenten
    print("\n1. Baue W02, WR, W' separat...")
    W02 = build_W02_matrix_mp(N, L_mp)
    WR = build_WR_matrix_mp(N, L_mp, parity="even", verbose=False)
    Wp = build_Wprime_matrix_mp(N, L_mp, chi_trivial, 1, verbose=False)

    # Schritt 2: Projiziere auf Even-Unterraum
    print("2. Even-Projektion...")

    # Erstelle Dummy-Matrix um U_even zu bekommen
    dummy = mp.matrix(size, size)
    _, U_even = project_to_parity(dummy, N, parity="even")

    # Projiziere jede Matrix
    W02_even = U_even.T * W02 * U_even
    WR_even = U_even.T * WR * U_even
    Wp_even = U_even.T * Wp * U_even

    # Symmetrisiere
    for mat in [W02_even, WR_even, Wp_even]:
        for i in range(dim):
            for j in range(i + 1, dim):
                avg = (mat[i, j] + mat[j, i]) / 2
                mat[i, j] = avg
                mat[j, i] = avg

    # Schritt 3: Analysiere W02_even — sollte Rang 1 sein
    print("\n3. W02_even Eigenwerte (sollte Rang 1 sein):")
    w02_evals, w02_evecs = diagonalize_mp(W02_even, verbose=False)

    nonzero_evals = [(k, float(w02_evals[k])) for k in range(dim)
                     if abs(float(w02_evals[k])) > 1e-20]
    print(f"   Nicht-Null-Eigenwerte: {len(nonzero_evals)}")
    for k, ev in nonzero_evals:
        print(f"     k={k}: eigenvalue = {ev:+.6e}")

    # u_tilde = Eigenvektor zum groessten Eigenwert
    max_k = max(range(dim), key=lambda i: float(w02_evals[i]))
    u_tilde = mp.matrix(dim, 1)
    for r in range(dim):
        u_tilde[r, 0] = w02_evecs[r, max_k]
    C_tilde = float(w02_evals[max_k])
    print(f"\n   W02-Richtung: u_tilde (k={max_k}), C_tilde = {C_tilde:+.6e}")

    # Schritt 4: Baue volle Matrix und Cluster
    print("\n4. Volle Matrix M_even = W02_even - WR_even - Wp_even")
    H_even = mp.matrix(dim, dim)
    M_even = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            H_even[i, j] = -WR_even[i, j] - Wp_even[i, j]
            M_even[i, j] = W02_even[i, j] + H_even[i, j]

    w_full, Q_full = diagonalize_mp(M_even, verbose=False)
    w_min = w_full[0]
    cluster = [k for k in range(dim) if abs(float(w_full[k] - w_min)) < 1e-10]
    print(f"   Cluster: {len(cluster)} Eigenvektoren bei w_min = {float(w_min):+.12e}")

    # Schritt 5: PHP — Restriction von H auf Nullraum von W02_even
    print("\n5. PHP: H eingeschraenkt auf Nullraum von W02_even")

    # Nullraum-Basis von W02_even: alle Eigenvektoren mit EW ≈ 0
    null_idx = [k for k in range(dim) if abs(float(w02_evals[k])) < 1e-20]
    print(f"   Nullraum-Dimension: {len(null_idx)}")

    # Projector P in der Eigenbasis von W02_even
    # PHP in der Nullraum-Basis
    null_dim = len(null_idx)
    PHP = mp.matrix(null_dim, null_dim)
    for i_idx, i in enumerate(null_idx):
        for j_idx, j in enumerate(null_idx):
            # <w02_evec_i, H, w02_evec_j>
            s = mp.mpf(0)
            for a in range(dim):
                for b in range(dim):
                    s += w02_evecs[a, i] * H_even[a, b] * w02_evecs[b, j]
            PHP[i_idx, j_idx] = s

    # Symmetrisiere
    for i in range(null_dim):
        for j in range(i + 1, null_dim):
            avg = (PHP[i, j] + PHP[j, i]) / 2
            PHP[i, j] = avg
            PHP[j, i] = avg

    w_php, Q_php = diagonalize_mp(PHP, verbose=False)

    print(f"\n   PHP Eigenwerte (unterste 10):")
    for k in range(min(10, null_dim)):
        gap = float(w_php[k + 1] - w_php[k]) if k + 1 < null_dim else 0
        print(f"     k={k}: w_php = {float(w_php[k]):+.12e}   gap = {gap:.3e}")

    # Vergleiche PHP-Eigenwerte mit Cluster-Eigenwerten
    print(f"\n6. Vergleich: Cluster-EW von M_even vs. PHP-EW")
    print(f"   {'k':>3s}  {'M_even EW':>16s}  {'PHP EW':>16s}  {'Differenz':>12s}")
    for k in range(min(len(cluster), null_dim)):
        m_ew = float(w_full[cluster[k]])
        php_ew = float(w_php[k])
        diff = abs(m_ew - php_ew)
        print(f"   {k:3d}  {m_ew:+16.12e}  {php_ew:+16.12e}  {diff:12.3e}")

    # Schritt 7: Kopplung <q_k, u_tilde> und <q_k, H, q_k>
    print(f"\n7. Kopplungs-Analyse der Cluster-Eigenvektoren:")
    print(f"   {'k':>3s}  {'<q_k,u_tilde>':>15s}  {'<q_k,H,q_k>':>15s}  {'F(g1)':>12s}  {'F(g2)':>12s}")

    for k in cluster:
        # <q_k, u_tilde>
        ov_u = mp.mpf(0)
        for r in range(dim):
            ov_u += Q_full[r, k] * u_tilde[r, 0]

        # <q_k, H_even, q_k>
        Hq = mp.matrix(dim, 1)
        for r in range(dim):
            s = mp.mpf(0)
            for c in range(dim):
                s += H_even[r, c] * Q_full[c, k]
            Hq[r, 0] = s
        qHq = mp.mpf(0)
        for r in range(dim):
            qHq += Q_full[r, k] * Hq[r, 0]

        # F-Werte
        xi_full = mp.matrix(size, 1)
        for row in range(size):
            s = mp.mpf(0)
            for r in range(dim):
                s += U_even[row, r] * Q_full[r, k]
            xi_full[row, 0] = s
        F1 = evaluate_F(xi_full, N, L_mp, RIEMANN_ZEROS[0])
        F2 = evaluate_F(xi_full, N, L_mp, RIEMANN_ZEROS[1])

        print(f"   {k:3d}  {float(ov_u):+15.6e}  {float(qHq):+15.6e}  {F1:+12.4e}  {F2:+12.4e}")

    # Schritt 8: W02-Anteil im Eigenvektor
    print(f"\n8. Zerlegung: Wie viel W02-Richtung steckt in jedem Cluster-EV?")
    print(f"   W02 = C_tilde * |u_tilde><u_tilde|,  C_tilde = {C_tilde:+.6e}")
    print(f"   Fuer jeden Cluster-EV q_k: ||P q_k||^2 = 1 - |<q_k,u_tilde>|^2")
    for k in cluster:
        ov = mp.mpf(0)
        for r in range(dim):
            ov += Q_full[r, k] * u_tilde[r, 0]
        null_comp = 1 - float(ov)**2
        print(f"   k={k}: |<q_k,u_tilde>|^2 = {float(ov)**2:.6e}   "
              f"||P q_k||^2 = {null_comp:.6e}")


if __name__ == "__main__":
    main()
