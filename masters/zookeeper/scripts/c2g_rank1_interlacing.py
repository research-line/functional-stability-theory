"""
c2g_rank1_interlacing.py — Rang-1-Stoerungstheorie und Cauchy-Interlacing

A = H + C_tilde |u_tilde><u_tilde| ist Rang-1-Perturbation von H.
Aequivalent: A hat Block-Zerlegung in {u_tilde, e_1,...,e_N} Basis:

    A = [ uHu + C_tilde   f^T       ]
        [ f               diag(t_j)  ]

Secular-Gleichung: uHu + C_tilde - w - sum f_j^2/(t_j - w) = 0

Cauchy-Interlacing: w_min < t_1 < w_2 < t_2 < ... < w_N < t_N < w_{N+1}

=> A-Eigenvektor bei w_min hat P_0-Komponenten xi_j = -f_j * alpha / (t_j - w_min)
=> Da w_min < t_j: sgn(xi_j) = -sgn(f_j)   [ERZWUNGEN durch Interlacing]
=> phi_j = xi_j/f_j = -alpha/(t_j - w_min) < 0   [AUTOMATISCH]

Test: Gilt das auch fuer k_full (nicht nur fuer xi_0)?
=> Ja, wenn k_full genuegend nahe an xi_0 ist (dominanter c_0-Koeffizient)

Autor: LG (Opus 4.6, Session 27)
Datum: 2026-04-19
"""

import os, sys, time
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, diagonalize_mp, chi_trivial,
)
from c2_approximation_test import (
    k_lambda_value, inner_product, norm,
)
from c2_poisson_decomposition import project_to_fourier

DPS = int(os.environ.get("DPS", 50))
CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
    {"lam": 7.0, "N": 85},
]


def run(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    print(f"\n{'='*72}")
    print(f"  RANG-1 INTERLACING:  lambda={lam}, N={N}, dim={dim}")
    print(f"{'='*72}")

    # === Build operators ===
    t0 = time.time()
    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    W02 = mp.matrix(dim, dim)
    H = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            W02[i, j] = Aq[i, j] - Ah[i, j]
            H[i, j] = Ah[i, j]

    col0 = mp.matrix(dim, 1)
    for i in range(dim):
        col0[i, 0] = W02[i, 0]
    cn = norm(col0, dim)
    ut = mp.matrix(dim, 1)
    for i in range(dim):
        ut[i, 0] = col0[i, 0] / cn
    C_tilde = float(cn * cn)

    # H u_tilde
    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))
    print(f"  Build: {time.time()-t0:.0f}s, C_tilde={C_tilde:.4f}, uHu={uHu:.6f}")

    # === k_full ===
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
    alpha = float(inner_product(ut, kn, dim))

    # f = P_0 H ut, g = P_0 k
    f_lam = mp.matrix(dim, 1)
    g_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]
        g_lam[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]

    # === Diagonalize T and A ===
    T = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            T[i, j] = (H[i, j]
                        - ut[i, 0] * Hut[j, 0]
                        - Hut[i, 0] * ut[j, 0]
                        + mp.mpf(uHu) * ut[i, 0] * ut[j, 0])

    t0 = time.time()
    ts, Es = diagonalize_mp(T, verbose=False)
    ws, Qs = diagonalize_mp(Aq, verbose=False)
    print(f"  Diag T+A: {time.time()-t0:.0f}s")

    # T-eigenvalues (sorted, skip zero-mode)
    t_vals = []
    t_evecs = []
    for j in range(dim):
        tj = float(ts[j])
        if abs(tj) < 1e-12:
            continue
        ev = mp.matrix(dim, 1)
        for i in range(dim):
            ev[i, 0] = Es[i, j]
        t_vals.append(tj)
        t_evecs.append(ev)

    # A-eigenvalues (sorted)
    w_vals = [float(ws[k]) for k in range(dim)]
    w_min = w_vals[0]

    # f_j, g_j in T-eigenbasis
    fj_arr = []
    gj_arr = []
    for idx, ev in enumerate(t_evecs):
        fj_arr.append(float(inner_product(ev, f_lam, dim)))
        gj_arr.append(float(inner_product(ev, g_lam, dim)))

    n_t = len(t_vals)
    n_w = len(w_vals)

    print(f"\n  alpha = {alpha:.6f}")
    print(f"  w_min = {w_min:.6f}")
    print(f"  t_1   = {t_vals[0]:.6f}")
    print(f"  Gap w_min < t_1: {t_vals[0] - w_min:.6e}")

    # === PART 1: INTERLACING VERIFICATION ===
    print(f"\n  === CAUCHY-INTERLACING ===")
    print(f"  Pattern: w_1 < t_1 < w_2 < t_2 < ... < w_N < t_N < w_{n_t+1}")
    print(f"  {'k':>4}  {'w_k':>12}  {'t_k':>12}  {'w_{k+1}':>12}  {'OK?':>4}")
    print(f"  {'-'*52}")

    violations = 0
    for k in range(min(n_t, n_w - 1)):
        wk = w_vals[k]
        tk = t_vals[k] if k < n_t else None
        wk1 = w_vals[k + 1] if k + 1 < n_w else None

        if tk is not None and wk1 is not None:
            ok = wk < tk < wk1
            if not ok:
                violations += 1
            marker = "OK" if ok else "FAIL"
        else:
            marker = "---"

        t_str = f"{tk:12.6f}" if tk is not None else f"{'---':>12}"
        w1_str = f"{wk1:12.6f}" if wk1 is not None else f"{'---':>12}"
        if k < 8 or k > n_t - 3:
            print(f"  {k:4d}  {wk:12.6f}  {t_str}  {w1_str}  {marker:>4}")
        elif k == 8:
            print(f"  {'...':>4}")

    print(f"\n  Violations: {violations}/{min(n_t, n_w-1)}")

    # === PART 2: SECULAR EQUATION ===
    print(f"\n  === SECULAR-GLEICHUNG ===")
    print(f"  S(w) = uHu + C_tilde - w - sum f_j^2/(t_j - w)")

    for k in range(min(5, n_w)):
        w = w_vals[k]
        S_val = uHu + C_tilde - w
        for idx in range(n_t):
            denom = t_vals[idx] - w
            if abs(denom) > 1e-30:
                S_val -= fj_arr[idx]**2 / denom
        print(f"  S(w_{k}) = S({w:.6f}) = {S_val:.4e}")

    # === PART 3: EIGENVECTOR RECONSTRUCTION ===
    print(f"\n  === RANG-1 EIGENVEKTOR-REKONSTRUKTION ===")
    print(f"  xi_j^theory = -f_j * alpha / (t_j - w_min)")
    print(f"  vs. g_j = <e_j, P_0 k_full> (beobachtet)")

    # Lowest A-eigenvector
    xi0 = mp.matrix(dim, 1)
    for i in range(dim):
        xi0[i, 0] = Qs[i, 0]
    xi0_u = float(inner_product(ut, xi0, dim))

    print(f"\n  xi_0 u-Komponente: {xi0_u:.6f} (vs alpha={alpha:.6f})")

    print(f"\n  {'j':>4}  {'t_j':>10}  {'f_j':>10}  {'g_j':>10}  {'xi_j^R1':>10}  "
          f"{'xi_j^A':>10}  {'g/xi_R1':>10}  {'sgn_match':>10}")
    print(f"  {'-'*88}")

    sign_match_count = 0
    sign_total = 0
    ratio_g_xi = []

    for idx in range(n_t):
        tj = t_vals[idx]
        fj = fj_arr[idx]
        gj = gj_arr[idx]

        denom = tj - w_min
        if abs(denom) < 1e-12:
            continue

        # Rank-1 forced form: xi_j = -f_j * xi0_u / (t_j - w_min)
        xi_j_r1 = -fj * alpha / denom

        # Actual A-eigenvector component
        xi_j_actual = float(inner_product(t_evecs[idx], xi0, dim))

        # Ratio g_j / xi_j_r1
        if abs(xi_j_r1) > 1e-30:
            g_over_xi = gj / xi_j_r1
        else:
            g_over_xi = float('nan')

        # Sign match
        if abs(gj) > 1e-25 and abs(xi_j_r1) > 1e-25:
            sign_total += 1
            sm = (gj * xi_j_r1 > 0)
            if sm:
                sign_match_count += 1
        else:
            sm = None

        sm_str = "YES" if sm else ("NO" if sm is not None else "---")

        if idx < 10 or idx > n_t - 5 or (sm is not None and not sm):
            print(f"  {idx:4d}  {tj:10.4f}  {fj:10.4e}  {gj:10.4e}  "
                  f"{xi_j_r1:10.4e}  {xi_j_actual:10.4e}  "
                  f"{g_over_xi:10.6f}  {sm_str:>10}")

        if abs(g_over_xi) < 1e10 and not np.isnan(g_over_xi):
            ratio_g_xi.append(g_over_xi)

    print(f"\n  Vorzeichen-Match g_j vs xi_j^R1: {sign_match_count}/{sign_total}")

    if ratio_g_xi:
        arr = np.array(ratio_g_xi)
        print(f"  g_j/xi_j^R1 Statistik ({len(arr)} Moden):")
        print(f"    Median = {np.median(arr):.6f}")
        print(f"    Mean   = {np.mean(arr):.6f}")
        print(f"    Std    = {np.std(arr):.6f}")
        print(f"    Min    = {np.min(arr):.6f}")
        print(f"    Max    = {np.max(arr):.6f}")

    # === PART 4: SIGN-FORCING ANALYSE ===
    print(f"\n  === VORZEICHEN-FORCING ANALYSE ===")
    print(f"  Frage: Wie robust ist sgn(g_j) = -sgn(f_j)?")
    print(f"  Rang-1 GARANTIERT dies fuer xi_0 (da w_min < t_j => t_j-w_min > 0)")
    print(f"  Fuer k_full gilt es, wenn |c_0| die Summe |c_k>=1| dominiert\n")

    # k_full in A-eigenbasis
    c_coeffs = []
    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]
        ck = float(inner_product(qk, kn, dim))
        c_coeffs.append(ck)

    c0 = abs(c_coeffs[0])
    c_rest = np.sqrt(sum(c**2 for c in c_coeffs[1:]))
    c_rest_max = max(abs(c) for c in c_coeffs[1:])

    print(f"  k_full in A-Eigenbasis:")
    print(f"    |c_0| = {c0:.6f}  (Grundzustand)")
    print(f"    ||c_rest|| = {c_rest:.6f}  (alle k>=1)")
    print(f"    max|c_k| (k>=1) = {c_rest_max:.6f}")
    print(f"    |c_0|/||c_rest|| = {c0/c_rest:.4f}")

    # Per-mode sign stability
    print(f"\n  Per-Moden Vorzeichen-Stabilitaet:")
    print(f"  {'j':>4}  {'|c0*xi0_j|':>12}  {'sum|ck*xik_j|':>14}  {'Ratio':>10}  {'sgn_safe':>10}")
    print(f"  {'-'*58}")

    sign_safe = 0
    sign_unsafe = 0
    min_ratio = 1e30

    for idx in range(n_t):
        tj = t_vals[idx]
        fj = fj_arr[idx]
        if abs(fj) < 1e-25:
            continue

        # c_0 * xi0_j contribution
        xi0_j = float(inner_product(t_evecs[idx], xi0, dim))
        dominant = abs(c_coeffs[0] * xi0_j)

        # sum_{k>=1} |c_k * xi_k_j|
        perturb = 0.0
        for k in range(1, dim):
            qk = mp.matrix(dim, 1)
            for i in range(dim):
                qk[i, 0] = Qs[i, k]
            xik_j = float(inner_product(t_evecs[idx], qk, dim))
            perturb += abs(c_coeffs[k] * xik_j)

        ratio = dominant / perturb if perturb > 1e-30 else float('inf')
        safe = ratio > 1.0

        if safe:
            sign_safe += 1
        else:
            sign_unsafe += 1

        if ratio < min_ratio:
            min_ratio = ratio

        if idx < 8 or not safe or idx > n_t - 3:
            safe_str = "SAFE" if safe else "UNSAFE"
            print(f"  {idx:4d}  {dominant:12.4e}  {perturb:14.4e}  {ratio:10.4f}  {safe_str:>10}")
        elif idx == 8:
            print(f"  {'...':>4}")

    print(f"\n  Sign-safe Moden: {sign_safe}/{sign_safe+sign_unsafe}")
    print(f"  Min Dominanz-Ratio: {min_ratio:.4f}")
    print(f"  => Ratio > 1 bedeutet: c_0*xi0_j dominiert => Vorzeichen ERZWUNGEN")

    # === PART 5: INTERLACING => SIGN THEOREM ===
    print(f"\n  === ZUSAMMENFASSUNG: INTERLACING => VORZEICHEN ===")
    print(f"  1. A = H + C_tilde|u><u| (Rang-1)")
    print(f"  2. Cauchy-Interlacing: w_min < t_1 < ... (verifiziert, {violations} Verletzungen)")
    print(f"  3. => xi_0 hat sgn(xi_j) = -sgn(f_j) ERZWUNGEN")
    print(f"  4. k_full: |c_0|={c0:.4f}, ||c_rest||={c_rest:.4f}")
    print(f"  5. => Vorzeichen-Match: {sign_match_count}/{sign_total}")
    print(f"  6. => Sign-safe Moden: {sign_safe}/{sign_safe+sign_unsafe}")
    if violations == 0 and sign_match_count == sign_total:
        print(f"  ERGEBNIS: Rang-1 Interlacing erklaert phi_j < 0 VOLLSTAENDIG")


if __name__ == "__main__":
    for cfg in CONFIGS:
        run(cfg["lam"], cfg["N"])
    print(f"\n=== FERTIG ===")
