"""
weil_positivity_analysis.py — Numerische Verifikation des Weil-Positivitaets-Zusammenhangs.

Zentrale Frage: WARUM kodiert der Cluster-Nullraum die Nullstellen?

Methodik:
  1. Evaluationsvektoren e(gamma_j) im Even-Unterraum berechnen
  2. Zerlegung: e(gamma) = <u_tilde, e> u_tilde  +  P e(gamma)   (Pol + Nullraum)
  3. F(gamma_j) = <e(gamma_j), q_k> fuer jeden Cluster-EV q_k
  4. Gram-Matrix G_ij = <P e(gamma_i), P e(gamma_j)> — wie orthogonal sind die Nullraum-Anteile?
  5. Spektrale Approximation: PHP ≈ Sigma_rho |Pc(rho)><Pc(rho)|
  6. Positvitaets-Bilanz: W02 vs. PHP Eigenwerte

Theorem (Weil, 1952): QW(f,f) >= 0 fuer alle zul. f  <=>  RH.
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


def build_eval_vector_even(gamma, N, L_mp):
    """Evaluationsvektor e(gamma) im Even-Unterraum (dim N+1).
    KORREKT: F(gamma) = xi_0/gamma + (1/sqrt2) sum_{n>=1} xi_n * 2gamma/(gamma^2 - p_n^2)
    wobei xi = U_even * q, also xi_0 = q[0], xi_n = q[n]/sqrt2.
    => F = q[0]/gamma + (1/sqrt2) sum q[n] * 2gamma/(gamma^2 - p_n^2)
    => e_0 = 1/gamma, e_n = (1/sqrt2) * 2gamma/(gamma^2 - p_n^2)
    """
    dim = N + 1
    g = mp.mpf(gamma)
    L = L_mp
    inv_sqrt2 = 1 / mp.sqrt(mp.mpf(2))
    e = mp.matrix(dim, 1)
    e[0, 0] = 1 / g
    for n in range(1, N + 1):
        pole = 2 * mp.pi * n / L
        e[n, 0] = inv_sqrt2 * 2 * g / (g**2 - pole**2)
    return e


def evaluate_F_fullspace(q_k, U_even, N, L_mp, gamma):
    """Referenzberechnung: F(gamma) im vollen Fourier-Raum."""
    size = 2 * N + 1
    dim = N + 1
    xi = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for r in range(dim):
            s += U_even[row, r] * q_k[r, 0]
        xi[row, 0] = s
    g = mp.mpf(gamma)
    F = mp.mpf(0)
    for idx in range(size):
        j = idx - N
        pole = 2 * mp.pi * j / L_mp
        F += xi[idx, 0] / (g - pole)
    return float(F)


def inner_product(a, b, dim):
    """<a, b> fuer mp.matrix Spaltenvektor."""
    s = mp.mpf(0)
    for i in range(dim):
        s += a[i, 0] * b[i, 0]
    return s


def mat_vec(M, v, dim):
    """M @ v fuer mp.matrix."""
    r = mp.matrix(dim, 1)
    for i in range(dim):
        s = mp.mpf(0)
        for j in range(dim):
            s += M[i, j] * v[j, 0]
        r[i, 0] = s
    return r


def main():
    mp.mp.dps = DPS
    L_mp = 2 * mp.log(mp.mpf(LAM))
    size = 2 * N + 1
    dim = N + 1

    print(f"Weil-Positivitaets-Analyse: lambda={LAM}, N={N}, dps={DPS}")
    print(f"L = {float(L_mp):.6f}")
    print("=" * 80)

    # === PHASE 1: Baue alle Matrizen ===
    print("\n[Phase 1] Matrix-Konstruktion...")
    W02 = build_W02_matrix_mp(N, L_mp)
    WR = build_WR_matrix_mp(N, L_mp, parity="even", verbose=False)
    Wp = build_Wprime_matrix_mp(N, L_mp, chi_trivial, 1, verbose=False)

    dummy = mp.matrix(size, size)
    _, U_even = project_to_parity(dummy, N, parity="even")

    W02_even = U_even.T * W02 * U_even
    WR_even = U_even.T * WR * U_even
    Wp_even = U_even.T * Wp * U_even

    for mat in [W02_even, WR_even, Wp_even]:
        for i in range(dim):
            for j in range(i + 1, dim):
                avg = (mat[i, j] + mat[j, i]) / 2
                mat[i, j] = avg
                mat[j, i] = avg

    # H und M
    H_even = mp.matrix(dim, dim)
    M_even = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            H_even[i, j] = -WR_even[i, j] - Wp_even[i, j]
            M_even[i, j] = W02_even[i, j] + H_even[i, j]

    # Diagonalisierungen
    w02_evals, w02_evecs = diagonalize_mp(W02_even, verbose=False)
    w_full, Q_full = diagonalize_mp(M_even, verbose=False)

    # u_tilde und C_tilde
    max_k = max(range(dim), key=lambda i: float(w02_evals[i]))
    u_tilde = mp.matrix(dim, 1)
    for r in range(dim):
        u_tilde[r, 0] = w02_evecs[r, max_k]
    C_tilde = float(w02_evals[max_k])

    # Cluster identifizieren
    w_min = w_full[0]
    cluster = [k for k in range(dim) if abs(float(w_full[k] - w_min)) < 1e-10]
    n_cluster = len(cluster)

    print(f"  W02: C_tilde = {C_tilde:+.6e} (Rang 1)")
    print(f"  Cluster: {n_cluster} Eigenvektoren bei w_min = {float(w_min):+.12e}")

    # === PHASE 2: Evaluationsvektoren an Riemann-Nullstellen ===
    print(f"\n[Phase 2] Evaluationsvektoren e(gamma_j) im Even-Unterraum")
    n_zeros = len(RIEMANN_ZEROS)
    eval_vecs = []
    for j, gamma in enumerate(RIEMANN_ZEROS):
        e_j = build_eval_vector_even(gamma, N, L_mp)
        eval_vecs.append(e_j)

    # === PHASE 3: Zerlegung e = <u,e> u  +  P e ===
    print(f"\n[Phase 3] Zerlegung der Evaluationsvektoren")
    print(f"  e(gamma) = <u_tilde, e> * u_tilde  +  P * e(gamma)")
    print(f"  {'gamma':>12s}  {'||e||^2':>12s}  {'|<u,e>|^2':>12s}  "
          f"{'||Pe||^2':>12s}  {'Pol-Anteil':>10s}")

    e_pole_overlaps = []
    e_null_components = []
    for j, gamma in enumerate(RIEMANN_ZEROS):
        e_j = eval_vecs[j]
        norm2 = float(inner_product(e_j, e_j, dim))
        overlap_u = float(inner_product(u_tilde, e_j, dim))
        pole_frac = overlap_u**2 / norm2

        Pe_j = mp.matrix(dim, 1)
        for r in range(dim):
            Pe_j[r, 0] = e_j[r, 0] - overlap_u * u_tilde[r, 0]
        Pe_norm2 = float(inner_product(Pe_j, Pe_j, dim))

        e_pole_overlaps.append(overlap_u)
        e_null_components.append(Pe_j)

        print(f"  {gamma:12.6f}  {norm2:12.4e}  {overlap_u**2:12.4e}  "
              f"{Pe_norm2:12.4e}  {100*pole_frac:8.2f}%")

    # === PHASE 3b: Verifikation ===
    print(f"\n[Phase 3b] Verifikation: <e_even, q> vs. Vollraum-F")
    print(f"  {'k':>3s}  {'<e,q> (even)':>15s}  {'F_full':>15s}  {'Match?':>8s}")
    for k in cluster[:2]:
        q_k = mp.matrix(dim, 1)
        for r in range(dim):
            q_k[r, 0] = Q_full[r, k]
        F_even = float(inner_product(eval_vecs[0], q_k, dim))
        F_full = evaluate_F_fullspace(q_k, U_even, N, L_mp, RIEMANN_ZEROS[0])
        match = "OK" if abs(F_even - F_full) < 1e-10 * max(abs(F_even), abs(F_full), 1e-30) else "FAIL"
        print(f"  {k:3d}  {F_even:+15.6e}  {F_full:+15.6e}  {match:>8s}")

    # === PHASE 4: F(gamma_j) Zerlegung fuer Cluster-EVs ===
    print(f"\n[Phase 4] F(gamma_j) = <e(gamma_j), q_k> Zerlegung")
    print(f"  F = <u,e>*<u,q>  +  <Pe, Pq>   (Pol-Term + Nullraum-Term)")
    print(f"\n  {'k':>3s} {'gamma':>9s}  {'F=<e,q>':>13s}  "
          f"{'Pol-Term':>13s}  {'Null-Term':>13s}  {'Cancel-Digits':>13s}")

    for k in cluster:
        q_k = mp.matrix(dim, 1)
        for r in range(dim):
            q_k[r, 0] = Q_full[r, k]
        uq = float(inner_product(u_tilde, q_k, dim))

        Pq_k = mp.matrix(dim, 1)
        for r in range(dim):
            Pq_k[r, 0] = q_k[r, 0] - uq * u_tilde[r, 0]

        for j, gamma in enumerate(RIEMANN_ZEROS[:3]):
            F_val = float(inner_product(eval_vecs[j], q_k, dim))
            pol_term = e_pole_overlaps[j] * uq
            null_term = float(inner_product(e_null_components[j], Pq_k, dim))
            if abs(F_val) > 1e-60 and abs(pol_term) > 1e-60:
                cancel = -np.log10(abs(F_val) / max(abs(pol_term), abs(null_term)))
            else:
                cancel = 0
            print(f"  {k:3d} {gamma:9.3f}  {F_val:+13.4e}  "
                  f"{pol_term:+13.4e}  {null_term:+13.4e}  {cancel:13.1f}")

    # === PHASE 5: Gram-Matrix der Nullraum-Evaluationsvektoren ===
    print(f"\n[Phase 5] Gram-Matrix G_ij = <P e(gamma_i), P e(gamma_j)>")
    G = np.zeros((n_zeros, n_zeros))
    for i in range(n_zeros):
        for j in range(n_zeros):
            G[i, j] = float(inner_product(e_null_components[i], e_null_components[j], dim))

    print(f"  Gram-Matrix (5x5):")
    header = "       " + "  ".join(f"{'g'+str(i+1):>10s}" for i in range(n_zeros))
    print(header)
    for i in range(n_zeros):
        row = f"  g{i+1}  " + "  ".join(f"{G[i,j]:+10.4e}" for j in range(n_zeros))
        print(row)

    G_evals = np.linalg.eigvalsh(G)
    print(f"\n  Gram-Eigenwerte: {['%.3e' % v for v in sorted(G_evals)]}")
    cond = G_evals[-1] / max(G_evals[0], 1e-50)
    print(f"  Konditionszahl: {cond:.2e}")

    # === PHASE 6: PHP-Spektrum und Evaluationsvektoren ===
    print(f"\n[Phase 6] PHP-Spektrum vs. Evaluationsvektor-Projektion")

    null_idx = [k for k in range(dim) if abs(float(w02_evals[k])) < 1e-20]
    null_dim = len(null_idx)

    PHP = mp.matrix(null_dim, null_dim)
    for i_idx, i in enumerate(null_idx):
        for j_idx, j in enumerate(null_idx):
            s = mp.mpf(0)
            for a in range(dim):
                for b in range(dim):
                    s += w02_evecs[a, i] * H_even[a, b] * w02_evecs[b, j]
            PHP[i_idx, j_idx] = s
    for i in range(null_dim):
        for j in range(i + 1, null_dim):
            avg = (PHP[i, j] + PHP[j, i]) / 2
            PHP[i, j] = avg
            PHP[j, i] = avg

    w_php, Q_php = diagonalize_mp(PHP, verbose=False)

    # Projiziere eval-Vektoren in den Nullraum (W02-Eigenbasis)
    print(f"  Projektion der Evaluationsvektoren in den W02-Nullraum (dim={null_dim}):")
    print(f"  {'gamma':>12s}  {'||P_null e||^2':>14s}  {'Overlap mit PHP-EV0..4':>30s}")

    eval_null = []
    for j, gamma in enumerate(RIEMANN_ZEROS):
        e_null_j = mp.matrix(null_dim, 1)
        for idx_k, k in enumerate(null_idx):
            s = mp.mpf(0)
            for r in range(dim):
                s += w02_evecs[r, k] * eval_vecs[j][r, 0]
            e_null_j[idx_k, 0] = s
        eval_null.append(e_null_j)
        norm2 = float(inner_product(e_null_j, e_null_j, null_dim))

        overlaps = []
        for kk in range(min(5, null_dim)):
            php_ev = mp.matrix(null_dim, 1)
            for r in range(null_dim):
                php_ev[r, 0] = Q_php[r, kk]
            ov = float(inner_product(e_null_j, php_ev, null_dim))
            overlaps.append(ov)

        ov_str = "  ".join(f"{v:+.3e}" for v in overlaps)
        print(f"  {gamma:12.6f}  {norm2:14.4e}  {ov_str}")

    # === PHASE 7: Weil-Positivitaets-Bilanz ===
    print(f"\n[Phase 7] Weil-Positivitaets-Bilanz")
    print(f"  QW = W02 + H = C_tilde |u><u| + H")
    print(f"  Weil-Kriterium: QW >= 0  <=>  RH")
    print(f"\n  Eigenwerte von M_even = QW (unterste 15 + oberster):")
    for k in range(min(15, dim)):
        print(f"    k={k:2d}: w = {float(w_full[k]):+16.12e}")
    print(f"    k={dim-1:2d}: w = {float(w_full[dim-1]):+16.12e}  (max)")

    n_neg = sum(1 for k in range(dim) if float(w_full[k]) < -1e-15)
    n_pos = sum(1 for k in range(dim) if float(w_full[k]) > 1e-15)
    n_zero = dim - n_neg - n_pos
    print(f"\n  Vorzeichen-Verteilung: {n_neg} negativ, {n_zero} null, {n_pos} positiv")
    print(f"  => QW ist {'NICHT positiv semidefinit' if n_neg > 0 else 'positiv semidefinit'} bei N={N}")

    # Relation: minimaler EW von QW vs. C_tilde
    print(f"\n  Positvitaets-Bedingung:")
    print(f"    C_tilde (Pol-Beitrag):     {C_tilde:+.6e}")
    print(f"    w_min (Cluster):           {float(w_min):+.6e}")
    print(f"    w_min / C_tilde:           {float(w_min)/C_tilde:+.6f}")
    print(f"    Abstand zu Positivitaet:   {-float(w_min):+.6e}")

    # === PHASE 8: Cancellation-Analyse ===
    print(f"\n[Phase 8] Cancellation-Mechanismus: Warum verschwindet F(gamma)?")
    print(f"  F(gamma) = [Pol-Term] + [Nullraum-Term]")
    print(f"  Pol-Term  = <u, e(gamma)> * <u, q_k>")
    print(f"  Null-Term = <Pe(gamma), Pq_k>")
    print(f"  Wenn Pol ≫ Null: Cancellation kommt aus <u,q_k> ~ 0")
    print(f"  Wenn Pol ~ Null: Exakte Cancellation zwischen beiden Termen\n")

    for k in cluster:
        q_k = mp.matrix(dim, 1)
        for r in range(dim):
            q_k[r, 0] = Q_full[r, k]
        uq = float(inner_product(u_tilde, q_k, dim))

        print(f"  Cluster-EV k={k}:  <u, q_k> = {uq:+.6e},  |<u,q>|^2 = {uq**2:.6e}")

        for j in range(min(3, n_zeros)):
            gamma = RIEMANN_ZEROS[j]
            F_val = float(inner_product(eval_vecs[j], q_k, dim))
            pol = e_pole_overlaps[j] * uq
            null = F_val - pol
            if abs(null) > 1e-60:
                cancel_ratio = abs(pol) / abs(null)
                cancel_digits = -np.log10(abs(F_val) / abs(pol)) if abs(pol) > 1e-60 else 0
                print(f"    gamma_{j+1}={gamma:.3f}: F={F_val:+.3e}  "
                      f"Pol={pol:+.3e}  Null={null:+.3e}  "
                      f"Cancellation: {cancel_digits:.1f} Stellen")

    # === PHASE 9: Kontrafaktische Analyse ===
    print(f"\n[Phase 9] Kontrafaktisch: Was waere bei Nicht-Nullstelle?")
    test_points = [10.0, 17.5, 23.0, 28.0, 35.0]
    print(f"  F(z) an Nicht-Nullstellen vs. Nullstellen fuer k={cluster[-1]}:")
    print(f"  {'z':>8s}  {'F(z)':>13s}  {'Typ':>12s}")

    k_best = cluster[-1]
    q_best = mp.matrix(dim, 1)
    for r in range(dim):
        q_best[r, 0] = Q_full[r, k_best]

    for z in test_points:
        e_z = build_eval_vector_even(z, N, L_mp)
        F_z = float(inner_product(e_z, q_best, dim))
        print(f"  {z:8.3f}  {F_z:+13.4e}  Nicht-Null")
    for z in RIEMANN_ZEROS:
        e_z = build_eval_vector_even(z, N, L_mp)
        F_z = float(inner_product(e_z, q_best, dim))
        print(f"  {z:8.3f}  {F_z:+13.4e}  NULLSTELLE")

    # === PHASE 10: N-Abhaengigkeit von QW-Positivitaet ===
    print(f"\n[Phase 10] Gesamt-Zusammenfassung")
    print(f"  Die Cluster-Eigenvektoren kodieren die Nullstellen durch die")
    print(f"  Weil-Explizitformel: QW = W02 + H enthaelt die arithmetische")
    print(f"  Information (Primzahlen) in H, die ueber die Explizitformel")
    print(f"  mit den Nullstellen von zeta(s) verknuepft ist.")
    print(f"\n  W02 (Pol bei s=1) erzeugt die strukturelle Degeneracy.")
    print(f"  PHP (H im Nullraum) traegt die Nullstellen-Information.")
    print(f"  F(gamma) ~ 0 weil die Cluster-EVs orthogonal zu den")
    print(f"  Evaluationsvektoren e(gamma) sind — eine Konsequenz der")
    print(f"  Eigenwertgleichung QW*q = w*q und der Weil-Formel.")


if __name__ == "__main__":
    main()
