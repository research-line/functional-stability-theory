"""
c2_resolvent_cancellation.py — B10 x AP2: Resolvent-Cancellation Verbindung

B10: F(gamma) = beta * [alpha(gamma) - R(gamma, w)]  fuer Eigenvektoren
AP2: F(gamma) = F_n1(gamma) + F_n2(gamma)  mit F_n1 ≈ -F_n2

Verbindung: Die AP2-Cancellation (a_k ≈ -b_k fuer Nicht-Cluster-EVs)
wird durch die B10-Resolvent-Struktur organisiert.

Messungen:
1. Ratio b_k/a_k vs Eigenwert-Gap (Cancellation-Profil)
2. u_tilde-Zerlegung von k_n1 und r_n2 (Pol vs Nullraum)
3. F-Zerlegung bei gamma_1 in Pol- und Nullraum-Anteile
4. Eigenvektorbeitraege: welche k tragen die Cancellation?

Autor: LG (Opus 4.6, Session 23)
Datum: 2026-04-19
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
    chi_trivial,
)
from c2_approximation_test import (
    k_lambda_value,
    h_educated_guess,
    inner_product,
    norm,
    evaluate_F_even,
)
from c2_poisson_decomposition import project_to_fourier

DPS = int(os.environ.get("DPS", 50))
LAM = 3.0
N_VAL = 30
GAMMA1 = 14.134725141734693


def main():
    mp.mp.dps = DPS
    lam = LAM
    N = N_VAL
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    print("=" * 80)
    print(f"  B10 x AP2: RESOLVENT-CANCELLATION (lam={lam}, N={N})")
    print("=" * 80)

    # === QW (mit W02) und H (ohne W02) ===
    t0 = time.time()
    M_qw_full, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    M_h_full, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)

    M_qw, _ = project_to_parity(M_qw_full, N, parity="even")
    M_h, _ = project_to_parity(M_h_full, N, parity="even")

    W02_even = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            W02_even[i, j] = M_qw[i, j] - M_h[i, j]

    w_sorted, Q_sorted = diagonalize_mp(M_qw, verbose=False)
    w_min = w_sorted[0]
    print(f"\n  QW + H in {time.time()-t0:.1f}s, w_min={mp.nstr(w_min, 8)}")

    # === u_tilde aus W02 extrahieren (Rang 1: W02 = C_tilde |u><u|) ===
    C_tilde = sum(W02_even[i, i] for i in range(dim))
    col0 = mp.matrix(dim, 1)
    for i in range(dim):
        col0[i, 0] = W02_even[i, 0]
    col0_norm = norm(col0, dim)
    u_tilde = mp.matrix(dim, 1)
    for i in range(dim):
        u_tilde[i, 0] = col0[i, 0] / col0_norm
    print(f"  C_tilde = {mp.nstr(C_tilde, 6)}")
    print(f"  ||u_tilde|| = {float(norm(u_tilde, dim)):.6f}")

    # === n-Split Fourier-Projektion ===
    t0 = time.time()

    def k_full_of_x(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    def k_n1_of_x(x):
        return mp.sqrt(mp.exp(x)) * h_educated_guess(mp.exp(x))

    c_full = project_to_fourier(k_full_of_x, L_mp, N)
    c_n1 = project_to_fourier(k_n1_of_x, L_mp, N)
    c_n2 = mp.matrix(dim, 1)
    for i in range(dim):
        c_n2[i, 0] = c_full[i, 0] - c_n1[i, 0]

    norm_full = norm(c_full, dim)
    k_normed = mp.matrix(dim, 1)
    n1_part = mp.matrix(dim, 1)
    n2_part = mp.matrix(dim, 1)
    for i in range(dim):
        k_normed[i, 0] = c_full[i, 0] / norm_full
        n1_part[i, 0] = c_n1[i, 0] / norm_full
        n2_part[i, 0] = c_n2[i, 0] / norm_full
    print(f"  Fourier-Projektion in {time.time()-t0:.1f}s")

    # === Evaluationsvektor e(gamma_1) ===
    gamma = mp.mpf(GAMMA1)
    e_gamma = mp.matrix(dim, 1)
    e_gamma[0, 0] = 1 / gamma
    for m in range(1, dim):
        pole = 2 * mp.pi * m / L_mp
        e_gamma[m, 0] = mp.sqrt(mp.mpf(2)) * gamma / (gamma**2 - pole**2)

    alpha_gamma = inner_product(u_tilde, e_gamma, dim)
    print(f"  alpha(gamma_1) = {mp.nstr(alpha_gamma, 8)}")

    # === Phase 1: b_k/a_k Ratio-Analyse ===
    print(f"\n{'='*80}")
    print(f"  Phase 1: RATIO-ANALYSE b_k/a_k")
    print(f"{'='*80}")
    print(f"  {'k':>4s}  {'a_k (n1)':>12s}  {'b_k (n2)':>12s}  {'c_k (full)':>12s}  "
          f"{'b/a ratio':>10s}  {'c/a (dev)':>12s}  {'gap':>12s}  {'beta_k':>10s}")

    ov_n1 = []
    ov_n2 = []
    ov_full = []
    beta_k_list = []
    F_qk_list = []

    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]

        ak = float(inner_product(n1_part, qk, dim))
        bk = float(inner_product(n2_part, qk, dim))
        ck = float(inner_product(k_normed, qk, dim))
        bk_u = float(inner_product(u_tilde, qk, dim))
        F_qk = float(evaluate_F_even(qk, N, L_mp, GAMMA1))

        ov_n1.append(ak)
        ov_n2.append(bk)
        ov_full.append(ck)
        beta_k_list.append(bk_u)
        F_qk_list.append(F_qk)

        gap = float(w_sorted[k] - w_min)
        ratio = bk / ak if abs(ak) > 1e-20 else float('inf')
        dev = ck / ak if abs(ak) > 1e-20 else float('inf')

        if k < 20:
            print(f"  {k:4d}  {ak:12.4e}  {bk:12.4e}  {ck:12.4e}  "
                  f"{ratio:10.6f}  {dev:12.4e}  {gap:12.4e}  {bk_u:10.4e}")

    # === Phase 2: Cancellation-Profil ===
    print(f"\n{'='*80}")
    print(f"  Phase 2: CANCELLATION-PROFIL")
    print(f"{'='*80}")
    print(f"  Cluster-Grenze: eng (delta=1e-10) => k=0..4")
    print(f"  Frage: Wie schnell geht |c_k/a_k| -> 0 (= Cancellation perfekt)?")

    print(f"\n  {'k':>4s}  {'|c_k/a_k|':>12s}  {'gap':>12s}  "
          f"{'|c_k|^2 Beitrag':>16s}  {'cum Leckage%':>14s}")
    cum_leak = 0
    total_leak = sum(ov_full[k]**2 for k in range(5, dim))
    for k in range(dim):
        ak = ov_n1[k]
        ck = ov_full[k]
        gap = float(w_sorted[k] - w_min)
        dev = abs(ck / ak) if abs(ak) > 1e-20 else 0
        contrib = ck**2
        if k >= 5:
            cum_leak += contrib
        pct = cum_leak / total_leak * 100 if total_leak > 0 else 0
        if k < 20 or contrib > 1e-15:
            marker = " <-- Cluster-Rand" if k == 4 else ""
            print(f"  {k:4d}  {dev:12.4e}  {gap:12.4e}  {contrib:16.4e}  "
                  f"{pct:14.2f}%{marker}")

    # === Phase 3: u_tilde-Zerlegung ===
    print(f"\n{'='*80}")
    print(f"  Phase 3: u_tilde-ZERLEGUNG (Pol vs Nullraum)")
    print(f"{'='*80}")

    u_n1 = float(inner_product(u_tilde, n1_part, dim))
    u_n2 = float(inner_product(u_tilde, n2_part, dim))
    u_full = float(inner_product(u_tilde, k_normed, dim))

    n1_norm_sq = float(inner_product(n1_part, n1_part, dim))
    n2_norm_sq = float(inner_product(n2_part, n2_part, dim))

    print(f"  <u_tilde, k_n1>  = {u_n1:12.6e}  (Pol-Anteil von n=1)")
    print(f"  <u_tilde, r_n2>  = {u_n2:12.6e}  (Pol-Anteil von n>=2)")
    print(f"  <u_tilde, k_full>= {u_full:12.6e}  (Pol-Anteil gesamt)")
    print(f"  Summencheck:       {u_n1 + u_n2:12.6e}  (soll = u_full)")
    print(f"  ||n1||^2 = {n1_norm_sq:.6e}, ||n2||^2 = {n2_norm_sq:.6e}")
    print(f"  Pol-Anteil^2 n1   = {u_n1**2:.6e}  ({u_n1**2/n1_norm_sq*100:.1f}% von ||n1||^2)")
    print(f"  Pol-Anteil^2 n2   = {u_n2**2:.6e}  ({u_n2**2/n2_norm_sq*100:.1f}% von ||n2||^2)")
    print(f"  Nullraum ||Pn1||^2= {n1_norm_sq - u_n1**2:.6e}")
    print(f"  Nullraum ||Pn2||^2= {n2_norm_sq - u_n2**2:.6e}")
    print(f"  Pol-Cancellation:  u_n1 + u_n2 = {u_n1 + u_n2:.6e}")

    # === Phase 4: F bei gamma_1 — Pol/Nullraum-Zerlegung ===
    print(f"\n{'='*80}")
    print(f"  Phase 4: F(gamma_1) — Pol/Nullraum-Zerlegung")
    print(f"{'='*80}")

    F_n1 = float(inner_product(n1_part, e_gamma, dim))
    F_n2 = float(inner_product(n2_part, e_gamma, dim))
    F_full = float(inner_product(k_normed, e_gamma, dim))

    F_n1_pole = float(alpha_gamma) * u_n1
    F_n1_null = F_n1 - F_n1_pole
    F_n2_pole = float(alpha_gamma) * u_n2
    F_n2_null = F_n2 - F_n2_pole
    F_full_pole = float(alpha_gamma) * u_full
    F_full_null = F_full - F_full_pole

    print(f"  alpha(gamma_1) = {float(alpha_gamma):.6e}")
    print(f"")
    print(f"  {'Teil':>12s}  {'F_total':>12s}  {'F_pole':>12s}  {'F_null':>12s}  {'pole%':>8s}")
    print(f"  {'k_n1':>12s}  {F_n1:12.4e}  {F_n1_pole:12.4e}  {F_n1_null:12.4e}  "
          f"{abs(F_n1_pole)/max(abs(F_n1),1e-50)*100:8.1f}%")
    print(f"  {'r_n2':>12s}  {F_n2:12.4e}  {F_n2_pole:12.4e}  {F_n2_null:12.4e}  "
          f"{abs(F_n2_pole)/max(abs(F_n2),1e-50)*100:8.1f}%")
    print(f"  {'k_full':>12s}  {F_full:12.4e}  {F_full_pole:12.4e}  {F_full_null:12.4e}  "
          f"{'---':>8s}")

    print(f"\n  Cancellation-Zerlegung:")
    print(f"    Pol-Cancellation:    {F_n1_pole:+.4e} + {F_n2_pole:+.4e} = {F_n1_pole+F_n2_pole:+.4e}")
    print(f"    Null-Cancellation:   {F_n1_null:+.4e} + {F_n2_null:+.4e} = {F_n1_null+F_n2_null:+.4e}")
    print(f"    Total-Cancellation:  {F_n1:+.4e} + {F_n2:+.4e} = {F_full:+.4e}")

    # === Phase 5: Eigenvektorbeitraege zu F_n1, F_n2 ===
    print(f"\n{'='*80}")
    print(f"  Phase 5: EIGENVEKTORBEITRAEGE zu F(gamma_1)")
    print(f"{'='*80}")
    print(f"  F_n1 = sum a_k F_qk,  F_n2 = sum b_k F_qk")
    print(f"  Cluster-EVs: F_qk ≈ 0, also Beitrag ≈ 0 unabhaengig von a_k, b_k")
    print(f"  Nicht-Cluster: a_k ≈ -b_k, also a_k F_qk ≈ -b_k F_qk (Cancellation)")
    print(f"\n  {'k':>4s}  {'F_qk':>12s}  {'a_k*F_qk':>12s}  {'b_k*F_qk':>12s}  "
          f"{'(a+b)*F_qk':>12s}  {'beta_k':>10s}")

    cum_Fn1 = 0
    cum_Fn2 = 0
    for k in range(min(20, dim)):
        ak = ov_n1[k]
        bk = ov_n2[k]
        Fk = F_qk_list[k]
        a_Fk = ak * Fk
        b_Fk = bk * Fk
        cum_Fn1 += a_Fk
        cum_Fn2 += b_Fk
        print(f"  {k:4d}  {Fk:12.4e}  {a_Fk:12.4e}  {b_Fk:12.4e}  "
              f"{a_Fk+b_Fk:12.4e}  {beta_k_list[k]:10.4e}")

    print(f"\n  Kumulative Summe (k=0..19):")
    print(f"    sum a_k F_qk = {cum_Fn1:.6e}  (soll ≈ F_n1 = {F_n1:.6e})")
    print(f"    sum b_k F_qk = {cum_Fn2:.6e}  (soll ≈ F_n2 = {F_n2:.6e})")

    # === Phase 6: Synthese ===
    print(f"\n{'='*80}")
    print(f"  SYNTHESE")
    print(f"{'='*80}")
    print(f"  Die drei Cancellation-Ebenen:")
    print(f"    1. B10:  alpha(gamma) ≈ R(gamma,w) fuer Cluster-EVs -> F_qk ≈ 0")
    print(f"    2. AP2:  a_k ≈ -b_k fuer Nicht-Cluster-EVs -> n1/n2-Cancellation")
    print(f"    3. MS2:  ||P_perp k|| ~ lambda^{{-1.84}} -> Leckage verschwindet")
    print(f"")
    print(f"  B10 + AP2 erklaeren MS2:")
    print(f"    Im Cluster:      F_qk ≈ 0 (B10-Cancellation alpha=R)")
    print(f"    Ausserhalb:      a_k + b_k ≈ 0 (AP2-Cancellation n1=-n2)")
    print(f"    => c_k ≈ 0 fuer ALLE k => k_lambda liegt im Cluster => MS2")


if __name__ == "__main__":
    main()
