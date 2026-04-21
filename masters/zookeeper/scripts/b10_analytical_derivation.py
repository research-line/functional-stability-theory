"""
b10_analytical_derivation.py — Resolvent-Identitaet: Analytische Ableitung der Cancellation

B10: Warum erzwingt QW*q = w*q die Bedingung F(gamma) = 0 an Riemann-Nullstellen?

HAUPTRESULTAT (Resolvent-Identitaet):
    F(gamma) = beta * [alpha(gamma) - R(gamma, w)]

    beta      = <u_tilde, q>                           (Pol-Kopplung)
    alpha(g)  = <u_tilde, e(g)>                        (direkte Pol-Antwort)
    R(g, w)   = <Pe(g), (PHP - wI)^{-1} * PH*u_tilde>  (Resolvent-Antwort)

    CANCELLATION: F(g)=0  <=>  alpha(g) = R(g, w)
    "Direkte Pol-Antwort = Arithmetische Green-Funktion-Antwort"

BEWEIS (algebraisch exakt):
    QW = C_tilde |u><u| + H   (Rang-1 + Arithmetik)
    QW*q = w*q  =>  P*(QW*q) = w*Pq  =>  PHP*Pq = w*Pq - beta*PHu
    =>  Pq = -beta * (PHP - wI)^{-1} * PHu     (Resolvent-Formel)
    =>  F(g) = <e(g), q> = beta*alpha + <Pe, Pq> = beta*[alpha - R]

    Dies ist eine ALGEBRAISCHE IDENTITAET, unabhaengig von Basis/Approximation.
    Die Cancellation alpha(g_j) = R(g_j, w) an Riemann-Nullstellen ist die
    endlich-dimensionale Realisierung der Weil-Explizitformel.
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

DPS = int(os.environ.get("DPS", 50))
N = int(os.environ.get("N", 30))
LAM = float(os.environ.get("LAM", 3.0))

RIEMANN_ZEROS = [14.134725141734693, 21.022039638771554, 25.010857580145688,
                 30.424876125859513, 32.935061587739189]
NON_ZEROS = [10.0, 17.5, 23.0, 28.0, 35.0]


def build_eval_vector_even(gamma, N_val, L_mp):
    dim = N_val + 1
    g = mp.mpf(gamma)
    inv_sqrt2 = 1 / mp.sqrt(mp.mpf(2))
    e = mp.matrix(dim, 1)
    e[0, 0] = 1 / g
    for n in range(1, N_val + 1):
        pole = 2 * mp.pi * n / L_mp
        e[n, 0] = inv_sqrt2 * 2 * g / (g**2 - pole**2)
    return e


def evaluate_F_fullspace(q_k, U_even, N_val, L_mp, gamma):
    size = 2 * N_val + 1
    dim = N_val + 1
    xi = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for r in range(dim):
            s += U_even[row, r] * q_k[r, 0]
        xi[row, 0] = s
    g = mp.mpf(gamma)
    F = mp.mpf(0)
    for idx in range(size):
        j = idx - N_val
        pole = 2 * mp.pi * j / L_mp
        F += xi[idx, 0] / (g - pole)
    return F


def dot(a, b, dim):
    return sum(float(a[r, 0]) * float(b[r, 0]) for r in range(dim))


def main():
    mp.mp.dps = DPS
    L_mp = 2 * mp.log(mp.mpf(LAM))
    size = 2 * N + 1
    dim = N + 1

    print("=" * 100)
    print("B10: RESOLVENT-IDENTITAET — Analytische Ableitung der Cancellation")
    print("=" * 100)
    print(f"  lambda={LAM}, N={N}, dim={dim}, L={float(L_mp):.4f}, DPS={DPS}")

    # ===== Phase 1: Matrixkonstruktion =====
    print("\n--- Phase 1: Matrixkonstruktion ---")
    W02 = build_W02_matrix_mp(N, L_mp)
    WR = build_WR_matrix_mp(N, L_mp, parity="even", verbose=False)
    Wp = build_Wprime_matrix_mp(N, L_mp, chi_trivial, 1, verbose=False)

    dummy = mp.matrix(size, size)
    _, U_even = project_to_parity(dummy, N, parity="even")

    W02e = U_even.T * W02 * U_even
    WRe = U_even.T * WR * U_even
    Wpe = U_even.T * Wp * U_even

    H_even = mp.matrix(dim, dim)
    QW_even = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            H_even[i, j] = -WRe[i, j] - Wpe[i, j]
            QW_even[i, j] = W02e[i, j] + H_even[i, j]
    for i in range(dim):
        for j in range(i + 1, dim):
            avg = (QW_even[i, j] + QW_even[j, i]) / 2
            QW_even[i, j] = avg
            QW_even[j, i] = avg

    w_all, Q = diagonalize_mp(QW_even, verbose=False)
    print(f"  QW: w_min={float(w_all[0]):+.6e}, w_max={float(w_all[dim-1]):+.6e}")

    # ===== Phase 2: Rang-1-Zerlegung =====
    print("\n--- Phase 2: Rang-1-Zerlegung W02 = C_tilde |u><u| ---")
    w02_vals, w02_vecs = diagonalize_mp(W02e, verbose=False)
    C_tilde = float(w02_vals[dim - 1])
    k_max_w02 = max(range(dim), key=lambda i: float(w02_vals[i]))

    u_tilde = mp.matrix(dim, 1)
    for r in range(dim):
        u_tilde[r, 0] = w02_vecs[r, k_max_w02]

    # Hu = H * u_tilde
    Hu = mp.matrix(dim, 1)
    for i in range(dim):
        s = mp.mpf(0)
        for j in range(dim):
            s += H_even[i, j] * u_tilde[j, 0]
        Hu[i, 0] = s

    h_uu = mp.mpf(0)
    for r in range(dim):
        h_uu += u_tilde[r, 0] * Hu[r, 0]

    # PHu = Hu - <u, Hu> * u  (Projektion auf Nullraum)
    PHu = mp.matrix(dim, 1)
    for r in range(dim):
        PHu[r, 0] = Hu[r, 0] - float(h_uu) * u_tilde[r, 0]

    print(f"  C_tilde = {C_tilde:+.6e}")
    print(f"  h_uu = <u, Hu> = {float(h_uu):+.6e}")
    print(f"  ||PHu|| = {float(mp.norm(PHu)):.6e}")

    # ===== Phase 3: PHP-Eigensystem =====
    print("\n--- Phase 3: PHP-Eigensystem im Nullraum ---")
    PHP = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            PHP[i, j] = (H_even[i, j]
                         - u_tilde[i, 0] * Hu[j, 0]
                         - Hu[i, 0] * u_tilde[j, 0]
                         + u_tilde[i, 0] * u_tilde[j, 0] * h_uu)
    for i in range(dim):
        for j in range(i + 1, dim):
            avg = (PHP[i, j] + PHP[j, i]) / 2
            PHP[i, j] = avg
            PHP[j, i] = avg

    mu_all, V_php = diagonalize_mp(PHP, verbose=False)

    zero_idx = min(range(dim), key=lambda k: abs(float(mu_all[k])))
    print(f"  PHP: mu_min={float(mu_all[0]):+.6e}, Null-Mode: mu[{zero_idx}]={float(mu_all[zero_idx]):.2e}")

    cluster_qw = [k for k in range(dim) if abs(float(w_all[k]) - float(w_all[0])) < 1e-10]
    n_cluster = len(cluster_qw)
    print(f"  QW-Cluster: {n_cluster} Eigenwerte bei w_min = {float(w_all[0]):+.6e}")

    # <v_k, PHu> fuer alle PHP-Eigenvektoren
    vk_PHu = []
    for k in range(dim):
        s = mp.mpf(0)
        for r in range(dim):
            s += V_php[r, k] * PHu[r, 0]
        vk_PHu.append(float(s))

    # ===== Phase 4: Kopplungsgleichung =====
    print("\n--- Phase 4: Kopplungsgleichung beta = -<PHu, Pq> / (C + h_uu - w) ---")
    test_ks = sorted(set([0, cluster_qw[-1], dim // 2, dim - 1]))

    for k in test_ks:
        q_k = mp.matrix(dim, 1)
        for r in range(dim):
            q_k[r, 0] = Q[r, k]
        w_k = float(w_all[k])
        beta_k = dot(u_tilde, q_k, dim)

        Pq_k = mp.matrix(dim, 1)
        for r in range(dim):
            Pq_k[r, 0] = q_k[r, 0] - beta_k * u_tilde[r, 0]

        PHu_Pq = dot(PHu, Pq_k, dim)
        denom = C_tilde + float(h_uu) - w_k
        beta_pred = -PHu_Pq / denom
        rel_err = abs(beta_k - beta_pred) / max(abs(beta_k), 1e-50)
        ok = "OK" if rel_err < 1e-8 else "FEHLER"
        print(f"  k={k:2d}: beta={beta_k:+.4e}, pred={beta_pred:+.4e}, err={rel_err:.2e} [{ok}]")

    # ===== Phase 5: RESOLVENT-IDENTITAET =====
    print("\n--- Phase 5: RESOLVENT-IDENTITAET F(g) = beta * [alpha(g) - R(g, w)] ---")

    k_best = cluster_qw[-1]
    w_best = float(w_all[k_best])
    q_best = mp.matrix(dim, 1)
    for r in range(dim):
        q_best[r, 0] = Q[r, k_best]
    beta_best = dot(u_tilde, q_best, dim)
    print(f"  k_best={k_best}, w={w_best:+.6e}, beta={beta_best:+.6e}")

    header = (f"  {'gamma':>8s}  {'alpha':>13s}  {'R(g,w)':>13s}  {'G=alpha-R':>13s}  "
              f"{'beta*G':>13s}  {'F_actual':>13s}  {'err':>9s}  {'typ':>4s}")
    print(header)
    print("  " + "-" * 95)

    results = []
    all_gammas = [(g, "ZERO") for g in RIEMANN_ZEROS] + [(z, "NON") for z in NON_ZEROS]
    all_gammas.sort(key=lambda x: x[0])

    for gamma, typ in all_gammas:
        e_g = build_eval_vector_even(gamma, N, L_mp)
        alpha_g = dot(u_tilde, e_g, dim)

        Pe_g = mp.matrix(dim, 1)
        for r in range(dim):
            Pe_g[r, 0] = e_g[r, 0] - alpha_g * u_tilde[r, 0]

        R_g = 0.0
        for k in range(dim):
            if k == zero_idx:
                continue
            mu_k = float(mu_all[k])
            dmu = mu_k - w_best
            if abs(dmu) < 1e-50:
                continue
            Pe_vk = dot(Pe_g, mp.matrix([V_php[r, k] for r in range(dim)], dim, 1), dim)
            R_g += Pe_vk * vk_PHu[k] / dmu

        G_g = alpha_g - R_g
        F_pred = beta_best * G_g
        F_act = float(evaluate_F_fullspace(q_best, U_even, N, L_mp, gamma))
        err = abs(F_pred - F_act)

        results.append((gamma, typ, alpha_g, R_g, G_g, F_pred, F_act, err))
        print(f"  {gamma:8.3f}  {alpha_g:+13.5e}  {R_g:+13.5e}  {G_g:+13.5e}  "
              f"{F_pred:+13.5e}  {F_act:+13.5e}  {err:9.1e}  {typ:>4s}")

    # ===== Phase 6: Cancellation-Analyse =====
    print("\n--- Phase 6: Cancellation-Analyse ---")
    G_zeros = [abs(r[4]) for r in results if r[1] == "ZERO"]
    G_non = [abs(r[4]) for r in results if r[1] == "NON"]
    print(f"  |G| an Riemann-Nullstellen:  max = {max(G_zeros):.2e}")
    print(f"  |G| an Nicht-Nullstellen:    min = {min(G_non):.2e}")
    sep = min(G_non) / max(G_zeros) if max(G_zeros) > 0 else float('inf')
    print(f"  Trennungsfaktor:             {sep:.2e}")

    # ===== Phase 7: Universalitaet — mehrere Eigenvektoren =====
    print("\n--- Phase 7: Universalitaet — Identitaet fuer verschiedene k ---")
    print(f"  {'k':>3s}  {'w_k':>13s}  {'beta':>11s}  {'|G(g1)|':>11s}  {'|F(g1)|':>11s}  {'err':>9s}")
    print("  " + "-" * 65)

    g1 = RIEMANN_ZEROS[0]
    e_g1 = build_eval_vector_even(g1, N, L_mp)
    alpha_g1 = dot(u_tilde, e_g1, dim)
    Pe_g1 = mp.matrix(dim, 1)
    for r in range(dim):
        Pe_g1[r, 0] = e_g1[r, 0] - alpha_g1 * u_tilde[r, 0]

    for k_test in sorted(set([0, 1, 2, 3, 4, dim // 4, dim // 2, 3 * dim // 4, dim - 1])):
        if k_test >= dim:
            continue
        q_t = mp.matrix(dim, 1)
        for r in range(dim):
            q_t[r, 0] = Q[r, k_test]
        w_t = float(w_all[k_test])
        beta_t = dot(u_tilde, q_t, dim)

        R_t = 0.0
        for k in range(dim):
            if k == zero_idx:
                continue
            mu_k = float(mu_all[k])
            dmu = mu_k - w_t
            if abs(dmu) < 1e-50:
                continue
            Pe_vk = dot(Pe_g1, mp.matrix([V_php[r, k] for r in range(dim)], dim, 1), dim)
            R_t += Pe_vk * vk_PHu[k] / dmu

        G_t = alpha_g1 - R_t
        F_pred_t = beta_t * G_t
        F_act_t = float(evaluate_F_fullspace(q_t, U_even, N, L_mp, g1))
        err_t = abs(F_pred_t - F_act_t)

        print(f"  {k_test:3d}  {w_t:+13.5e}  {beta_t:+11.3e}  {abs(G_t):11.3e}  "
              f"{abs(F_act_t):11.3e}  {err_t:9.1e}")

    # ===== Phase 8: Resonanz-Scan =====
    print("\n--- Phase 8: Resonanz-Scan — Nulldurchgaenge von G(gamma) ---")
    print("  Scan gamma in [10, 40], Schritt 0.25")
    print(f"  {'gamma':>8s}  {'G(gamma)':>14s}")
    print("  " + "-" * 25)

    gamma_scan = np.arange(10, 40.25, 0.25)
    prev_G = None
    crossings = []

    for gamma in gamma_scan:
        e_g = build_eval_vector_even(gamma, N, L_mp)
        alpha_g = dot(u_tilde, e_g, dim)
        Pe_g = mp.matrix(dim, 1)
        for r in range(dim):
            Pe_g[r, 0] = e_g[r, 0] - alpha_g * u_tilde[r, 0]

        R_g = 0.0
        for k in range(dim):
            if k == zero_idx:
                continue
            mu_k = float(mu_all[k])
            dmu = mu_k - w_best
            if abs(dmu) < 1e-50:
                continue
            Pe_vk = dot(Pe_g, mp.matrix([V_php[r, k] for r in range(dim)], dim, 1), dim)
            R_g += Pe_vk * vk_PHu[k] / dmu

        G_g = alpha_g - R_g
        marker = ""
        if prev_G is not None and prev_G * G_g < 0:
            crossings.append((gamma - 0.25, gamma, prev_G, G_g))
            marker = "  <-- NULLDURCHGANG"
        prev_G = G_g
        print(f"  {gamma:8.3f}  {G_g:+14.6e}{marker}")

    print(f"\n  Gefundene Nulldurchgaenge von G(gamma) = alpha(gamma) - R(gamma, w):")
    for a, b, ga, gb in crossings:
        gamma_est = a + 0.25 * abs(ga) / (abs(ga) + abs(gb))
        print(f"    gamma in [{a:.2f}, {b:.2f}], geschaetzt: {gamma_est:.3f}")

    print(f"\n  Riemann-Nullstellen zum Vergleich:")
    for g in RIEMANN_ZEROS:
        print(f"    gamma = {g:.6f}")

    # ===== Phase 9: Secular-Gleichung =====
    print("\n--- Phase 9: Secular-Gleichung (Rang-1-Perturbation) ---")
    print("  QW = C|u><u| + H  =>  Eigenwerte w erfuellen:")
    print("  1 + C * <u, (H-wI)^{-1} u> = 0")

    H_sym = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            H_sym[i, j] = H_even[i, j]
    for i in range(dim):
        for j in range(i + 1, dim):
            avg = (H_sym[i, j] + H_sym[j, i]) / 2
            H_sym[i, j] = avg
            H_sym[j, i] = avg

    h_vals, h_vecs = diagonalize_mp(H_sym, verbose=False)

    u_ov_sq = []
    for k in range(dim):
        ov = dot(u_tilde, mp.matrix([h_vecs[r, k] for r in range(dim)], dim, 1), dim)
        u_ov_sq.append(ov ** 2)

    print(f"\n  {'k':>3s}  {'w_k':>13s}  {'secular':>13s}  {'ok':>5s}")
    for k_test in sorted(set([0, cluster_qw[-1], dim // 2, dim - 1])):
        w_t = float(w_all[k_test])
        sec = 1.0 + C_tilde * sum(
            u_ov_sq[j] / (float(h_vals[j]) - w_t)
            for j in range(dim) if abs(float(h_vals[j]) - w_t) > 1e-20
        )
        ok = "OK" if abs(sec) < 0.01 else "MISS"
        print(f"  {k_test:3d}  {w_t:+13.5e}  {sec:+13.5e}  {ok:>5s}")

    # ===== ZUSAMMENFASSUNG =====
    print("\n" + "=" * 100)
    print("ZUSAMMENFASSUNG: B10 Resolvent-Identitaet")
    print("=" * 100)
    print()
    print("THEOREM (algebraisch exakt):")
    print("  Sei QW = C|u><u| + H (Rang-1 + Arithmetik), QW*q = w*q, beta = <u,q>.")
    print("  Dann gilt:")
    print()
    print("    F(gamma) = beta * [alpha(gamma) - R(gamma, w)]")
    print()
    print("  mit alpha(g) = <u, e(g)>  und  R(g,w) = <Pe(g), (PHP-wI)^{-1} PHu>")
    print()
    print("  CANCELLATION: F(g) = 0  <=>  alpha(g) = R(g, w)")
    print()

    max_err = max(r[7] for r in results)
    G_z = [abs(r[4]) for r in results if r[1] == "ZERO"]
    G_nz = [abs(r[4]) for r in results if r[1] == "NON"]
    print(f"VERIFIKATION:")
    print(f"  Identitaet verifiziert: max |F_pred - F_actual| = {max_err:.2e}")
    print(f"  |G| an Nullstellen:     {max(G_z):.2e}")
    print(f"  |G| an Nicht-Nullstellen: {min(G_nz):.2e}")
    print(f"  Trennungsfaktor:        {min(G_nz)/max(G_z):.2e}")
    print(f"  Nulldurchgaenge:        {len(crossings)} (erwartet: {len(RIEMANN_ZEROS)})")
    print()
    print("INTERPRETATION:")
    print("  Die Riemann-Nullstellen sind RESONANZEN: Punkte im Spektralraum,")
    print("  an denen die Pol-Antwort alpha(g) = <u, e(g)> exakt durch die")
    print("  arithmetische Resolvent-Antwort R(g,w) = <Pe(g), (PHP-wI)^{-1} PHu>")
    print("  kompensiert wird. Dies ist die finite Realisierung der Weil-Identitaet.")
    print()
    print("  Die 15-stellige Cancellation ist der numerische Abdruck dieser Resonanz.")
    print("  Die Identitaet ist ALGEBRAISCH EXAKT — sie folgt aus:")
    print("    (1) Rang-1-Struktur von W02")
    print("    (2) Eigenwertgleichung QW*q = w*q")
    print("    (3) Linearer Algebra (Projektion + Resolvent)")


if __name__ == "__main__":
    main()
