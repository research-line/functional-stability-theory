"""
c2m_modal_comparison.py — Modale Boundary/Bulk-Vergleichsanalyse (C2m.B)

Prüft GPTs "Boundary Comparison Lemma":
  |⟨e_j, B_L/nf⟩| ≤ κ_bd · α · |f_j|   auf J_reg
  |⟨e_j, E_bulk⟩| ≤ κ_bulk · α · |f_j|  auf J_reg

Schlüsseldiagnostik:
  - Sind κ_bd, κ_bulk < 1 auf J_reg?
  - Ist κ_bd ungefähr konstant über j? (GPT: "selbe Shift-Familie")
  - Cancellation: b_j und E_bulk_j gegen-Vorzeichen?
  - Korrelation b_j mit f_j?

Autor: LG (Opus 4.6, Session 27)
Datum: 2026-04-19
"""

import os, sys, time, math
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
from c2m_boundary_defect_explicit import B_prime_L, B_R_L

DPS = int(os.environ.get("DPS", 50))
CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
]
NGRID = 200
F_THRESHOLD = 1e-20


def run(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

    print(f"\n{'='*72}")
    print(f"  MODAL COMPARISON: lambda={lam}, N={N}, L={L:.4f}")
    print(f"{'='*72}")

    # === Phase 1: Operatoren ===
    t0 = time.time()
    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    H = mp.matrix(dim, dim)
    W02 = mp.matrix(dim, dim)
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

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = float(norm(cf, dim))
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / mp.mpf(nf)
    alpha = float(inner_product(ut, kn, dim))

    ws_qw, _ = diagonalize_mp(Aq, verbose=False)
    w_min = float(ws_qw[0])

    # Hũ und PHP = P₀HP₀
    Hu = mp.matrix(dim, 1)
    for i in range(dim):
        Hu[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hu, dim))

    PHP = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            PHP[i, j] = (H[i, j]
                         - Hu[i, 0] * ut[j, 0]
                         - ut[i, 0] * Hu[j, 0]
                         + mp.mpf(uHu) * ut[i, 0] * ut[j, 0])

    ts, Vs = diagonalize_mp(PHP, verbose=False)

    # Residual
    Hk = mp.matrix(dim, 1)
    for i in range(dim):
        Hk[i, 0] = sum(H[i, j] * kn[j, 0] for j in range(dim))

    res = mp.matrix(dim, 1)
    for i in range(dim):
        res[i, 0] = Hk[i, 0] - mp.mpf(w_min) * kn[i, 0]

    res_u = float(inner_product(ut, res, dim))
    res_P0 = mp.matrix(dim, 1)
    for i in range(dim):
        res_P0[i, 0] = res[i, 0] - mp.mpf(res_u) * ut[i, 0]

    print(f"  Phase 1 in {time.time()-t0:.0f}s")
    print(f"  w_min={w_min:.6f}, alpha={alpha:.6f}, ||K_L||={nf:.6f}")
    print(f"  |<ut, res>| = {abs(res_u):.6e}")
    print(f"  ||P0 res||  = {float(norm(res_P0, dim)):.6e}")

    # === Phase 2: Modale Projektionen (f_j, res_j) ===
    t_vals = np.array([float(ts[j]) for j in range(dim)])
    f_arr = np.zeros(dim)
    res_arr = np.zeros(dim)

    for j in range(dim):
        vj = mp.matrix(dim, 1)
        for i in range(dim):
            vj[i, 0] = Vs[i, j]
        f_arr[j] = float(inner_product(vj, Hu, dim))
        res_arr[j] = float(inner_product(vj, res_P0, dim))

    j_zero = int(np.argmin(np.abs(t_vals)))

    # Cluster-Erkennung via Eigenvalue-Nähe zu w_min
    # PHP-Cluster = QW-Cluster projiziert auf Ran(P₀)
    t_offsets = np.abs(t_vals - w_min)
    t_offsets[j_zero] = 1e10
    t_sorted_off = np.sort(t_offsets)
    for gap_idx in range(len(t_sorted_off) - 1):
        if t_sorted_off[gap_idx + 1] / max(t_sorted_off[gap_idx], 1e-30) > 10:
            cluster_gap = (t_sorted_off[gap_idx] + t_sorted_off[gap_idx + 1]) / 2
            break
    else:
        cluster_gap = 1e-3

    j_cluster = [j for j in range(dim) if j != j_zero and t_offsets[j] < cluster_gap]
    j_reg = [j for j in range(dim) if j != j_zero and t_offsets[j] >= cluster_gap]
    j_sing = j_cluster

    print(f"  PHP zero-EV: index {j_zero} (t={t_vals[j_zero]:.2e})")
    print(f"  Cluster-Gap: {cluster_gap:.4e}, Cluster-Moden: {len(j_cluster)}")
    print(f"  J_reg: {len(j_reg)}, J_sing(Cluster): {len(j_sing)}")

    # === Phase 3: B_L Grid ===
    print(f"\n  --- B_L Grid ({NGRID} Punkte) ---")
    t2 = time.time()
    x_grid = np.linspace(0.01 * L, 0.99 * L, NGRID)
    B_values = np.zeros(NGRID)

    test_xs = [0.05, 0.30, 0.50, 0.70, 0.95]
    bp_max = 0.0
    for frac in test_xs:
        bp = float(B_prime_L(mp.mpf(frac * L), L_mp))
        bp_max = max(bp_max, abs(bp))
    skip_prime = bp_max < 1e-15
    if skip_prime:
        print(f"  B'_L vernachlässigbar (max={bp_max:.1e})")

    for i, x in enumerate(x_grid):
        x_mp = mp.mpf(x)
        bp = B_prime_L(x_mp, L_mp) if not skip_prime else mp.mpf(0)
        br1, br2 = B_R_L(x_mp, L_mp)
        B_values[i] = float(bp + br1 + br2)
        if (i + 1) % 50 == 0:
            print(f"    {i+1}/{NGRID} ({time.time()-t2:.0f}s)")

    print(f"  Grid in {time.time()-t2:.0f}s")

    # Fourier-Projektion von B_L, normiert durch ||K_L||
    dx = x_grid[1] - x_grid[0]
    inv_sqrt_L = 1.0 / np.sqrt(L)
    sqrt_2_over_L = np.sqrt(2.0 / L)

    b_fourier = np.zeros(dim)
    b_fourier[0] = np.sum(B_values * inv_sqrt_L) * dx
    for n in range(1, dim):
        basis_n = sqrt_2_over_L * np.cos(2.0 * np.pi * n * x_grid / L)
        b_fourier[n] = np.sum(B_values * basis_n) * dx
    b_fourier /= nf

    # P₀-Projektion
    b_mp = mp.matrix(dim, 1)
    for i in range(dim):
        b_mp[i, 0] = mp.mpf(b_fourier[i])
    b_u = float(inner_product(ut, b_mp, dim))
    for i in range(dim):
        b_mp[i, 0] -= mp.mpf(b_u) * ut[i, 0]

    # b_j = ⟨v_j, P₀(B_L/nf)⟩
    b_arr = np.zeros(dim)
    for j in range(dim):
        vj = mp.matrix(dim, 1)
        for i in range(dim):
            vj[i, 0] = Vs[i, j]
        b_arr[j] = float(inner_product(vj, b_mp, dim))

    ebulk_arr = res_arr - b_arr

    # === Phase 4: κ-Analyse ===
    print(f"\n  {'='*72}")
    print(f"  MODALE κ-ANALYSE  (Boundary Comparison Lemma)")
    print(f"  {'='*72}")

    kappa_tot = np.zeros(len(j_reg))
    kappa_bd = np.zeros(len(j_reg))
    kappa_blk = np.zeros(len(j_reg))

    for idx, j in enumerate(j_reg):
        denom = abs(alpha * f_arr[j])
        kappa_tot[idx] = abs(res_arr[j]) / denom
        kappa_bd[idx] = abs(b_arr[j]) / denom
        kappa_blk[idx] = abs(ebulk_arr[j]) / denom

    order = np.argsort(kappa_tot)[::-1]

    print(f"\n  Top-20 Moden nach κ_total (absteigend):")
    print(f"  {'j':>4}  {'t_j':>10}  {'|f_j|':>10}  {'|res_j|':>10}  {'|b_j|':>10}  {'|E_blk|':>10}  {'κ_tot':>8}  {'κ_bd':>8}  {'κ_blk':>8}")
    print(f"  {'-'*92}")

    for rank in range(min(20, len(order))):
        sidx = order[rank]
        j = j_reg[sidx]
        print(f"  {j:4d}  {t_vals[j]:10.4f}  {abs(f_arr[j]):10.3e}  "
              f"{abs(res_arr[j]):10.3e}  {abs(b_arr[j]):10.3e}  "
              f"{abs(ebulk_arr[j]):10.3e}  {kappa_tot[sidx]:8.2e}  "
              f"{kappa_bd[sidx]:8.2e}  {kappa_blk[sidx]:8.2e}")

    # Zusammenfassung
    print(f"\n  --- κ Summary ---")
    for name, arr in [("κ_total", kappa_tot), ("κ_bd", kappa_bd), ("κ_bulk", kappa_blk)]:
        i_max = np.argmax(arr)
        print(f"  {name:>10}:  max={arr[i_max]:.4e} (j={j_reg[i_max]}),  "
              f"median={np.median(arr):.4e},  mean={np.mean(arr):.4e}")

    all_lt_1 = np.all(kappa_tot < 1.0)
    print(f"\n  κ_total < 1 auf ALLEN J_reg-Moden?  {'JA' if all_lt_1 else 'NEIN'}")
    print(f"  Marge:  1 - max(κ_total) = {1.0 - np.max(kappa_tot):.4e}")

    # Vorzeichen-Analyse
    same_sign = sum(1 for j in j_reg if b_arr[j] * ebulk_arr[j] > 0)
    opp_sign = len(j_reg) - same_sign
    print(f"\n  --- Vorzeichen b_j vs E_bulk_j ---")
    print(f"  Gleich-VZ:  {same_sign}/{len(j_reg)}")
    print(f"  Gegen-VZ:   {opp_sign}/{len(j_reg)}  (= partielle Cancellation)")

    # Korrelationen
    f_reg = np.array([f_arr[j] for j in j_reg])
    b_reg = np.array([b_arr[j] for j in j_reg])
    eb_reg = np.array([ebulk_arr[j] for j in j_reg])
    r_reg = np.array([res_arr[j] for j in j_reg])

    if np.std(f_reg) > 0:
        corr_bf = np.corrcoef(f_reg, b_reg)[0, 1]
        corr_ef = np.corrcoef(f_reg, eb_reg)[0, 1]
        corr_rf = np.corrcoef(f_reg, r_reg)[0, 1]
        corr_be = np.corrcoef(b_reg, eb_reg)[0, 1]
        print(f"\n  --- Korrelationen ---")
        print(f"  corr(b_j, f_j):       {corr_bf:+.4f}  (GPT: selbe Shift-Familie → ≈1?)")
        print(f"  corr(E_bulk_j, f_j):  {corr_ef:+.4f}")
        print(f"  corr(res_j, f_j):     {corr_rf:+.4f}")
        print(f"  corr(b_j, E_bulk_j):  {corr_be:+.4f}")

    # κ_bd Verteilung (GPT: sollte ~konstant sein)
    print(f"\n  --- κ_bd Verteilung ---")
    print(f"  min:   {np.min(kappa_bd):.4e}")
    print(f"  max:   {np.max(kappa_bd):.4e}")
    print(f"  mean:  {np.mean(kappa_bd):.4e}")
    print(f"  std:   {np.std(kappa_bd):.4e}")
    cv = np.std(kappa_bd) / np.mean(kappa_bd) if np.mean(kappa_bd) > 0 else 0
    print(f"  CV:    {cv:.2f}  (→ 0 = perfekt konstant)")

    # Ratio b_j / f_j auf J_reg
    ratio_bf = b_reg / f_reg
    print(f"\n  --- Ratio b_j / f_j ---")
    print(f"  mean:  {np.mean(ratio_bf):.6e}")
    print(f"  std:   {np.std(ratio_bf):.6e}")
    print(f"  min:   {np.min(ratio_bf):.6e}")
    print(f"  max:   {np.max(ratio_bf):.6e}")
    print(f"  Ist b_j ∝ f_j?  CV = {np.std(ratio_bf)/abs(np.mean(ratio_bf)):.2f}" if abs(np.mean(ratio_bf)) > 0 else "")

    # === C2n.4: Globaler Norm-Bound ===
    norm_R = float(norm(res_P0, dim))

    # Optimales τ_λ: scanne alle |f_j| als Schwellenwert, minimiere κ_global
    f_abs_all = np.array([abs(f_arr[j]) for j in range(dim) if j != j_zero])
    f_unique = np.sort(np.unique(f_abs_all))

    print(f"\n  {'='*72}")
    print(f"  C2n.4 — GLOBAL NORM BOUND (Residual Norm Bound)")
    print(f"  {'='*72}")
    print(f"  ||R_λ|| = ||P₀ res||    = {norm_R:.6e}")
    print(f"  α_cl                      = {abs(alpha):.6f}")

    # Scan: für jeden möglichen τ den resultierenden κ und modeweisen max-κ berechnen
    print(f"\n  --- τ_λ-Scan (alle Schwellenwerte) ---")
    print(f"  {'τ_λ':>12}  {'|J_reg|':>7}  {'|J_sing|':>7}  {'κ_global':>12}  {'max κ_mode':>12}  {'ALLE<1':>6}")

    best_kappa = 1e30
    best_tau = None
    best_nreg = 0
    best_nsing = 0
    best_max_mode = 0

    for fi in range(len(f_unique)):
        tau_test = f_unique[fi]
        if tau_test < 1e-20:
            continue
        jreg_test = [j for j in range(dim) if j != j_zero and abs(f_arr[j]) >= tau_test]
        if len(jreg_test) == 0:
            continue
        tau_eff = np.min([abs(f_arr[j]) for j in jreg_test])
        kg = norm_R / (abs(alpha) * tau_eff)
        kap_modes = [abs(res_arr[j]) / abs(alpha * f_arr[j]) for j in jreg_test if abs(f_arr[j]) > 0]
        max_km = max(kap_modes) if kap_modes else 0
        all_ok = all(k < 1 for k in kap_modes)
        nsing = dim - 1 - len(jreg_test)

        if kg < 1 and all_ok and len(jreg_test) > best_nreg:
            best_kappa = kg
            best_tau = tau_eff
            best_nreg = len(jreg_test)
            best_nsing = nsing
            best_max_mode = max_km

        if tau_eff > 1e-8:
            print(f"  {tau_eff:12.4e}  {len(jreg_test):7d}  {nsing:7d}  {kg:12.4e}  {max_km:12.4e}  {'JA' if all_ok else 'NEIN':>6}")

    print(f"\n  --- Optimales τ_λ ---")
    if best_tau is not None:
        print(f"  τ_λ*                      = {best_tau:.6e}")
        print(f"  J_reg: {best_nreg},  J_sing: {best_nsing}")
        print(f"  κ_global*                 = {best_kappa:.6e}")
        print(f"  max κ_mode                = {best_max_mode:.6e}")
        print(f"  κ < 1?  {'JA' if best_kappa < 1 else 'NEIN'}  (Marge: {1-best_kappa:.4e})")
    else:
        print(f"  Kein τ_λ mit κ < 1 gefunden!")

    tau_min = best_tau if best_tau else np.min(f_abs_all[f_abs_all > 1e-20])
    kappa_global = best_kappa if best_tau else norm_R / (abs(alpha) * tau_min)
    j_reg_tau = [j for j in range(dim) if j != j_zero and abs(f_arr[j]) >= tau_min]
    j_sing_tau = [j for j in range(dim) if j != j_zero and abs(f_arr[j]) < tau_min]

    return {
        "lam": lam,
        "n_reg": len(j_reg_tau),
        "n_sing": len(j_sing_tau),
        "max_kappa_total": float(best_max_mode),
        "max_kappa_bd": float(np.max(kappa_bd)),
        "max_kappa_bulk": float(np.max(kappa_blk)),
        "median_kappa_bd": float(np.median(kappa_bd)),
        "corr_bf": float(corr_bf) if np.std(f_reg) > 0 else 0.0,
        "cv_kappa_bd": float(cv),
        "kappa_global": float(kappa_global),
        "tau_lambda": float(tau_min),
    }


if __name__ == "__main__":
    results = []
    for cfg in CONFIGS:
        r = run(cfg["lam"], cfg["N"])
        results.append(r)

    if len(results) >= 2:
        print(f"\n{'='*72}")
        print(f"  λ-VERGLEICH")
        print(f"{'='*72}")
        print(f"\n  {'Metrik':>20}  ", end="")
        for r in results:
            print(f"{'lam='+str(r['lam']):>14}  ", end="")
        print(f"{'Exponent':>10}")

        for name, key in [
            ("max κ_total", "max_kappa_total"),
            ("κ_global (C2n.4)", "kappa_global"),
            ("τ_λ", "tau_lambda"),
            ("max κ_bd", "max_kappa_bd"),
            ("max κ_bulk", "max_kappa_bulk"),
            ("corr(b,f)", "corr_bf"),
            ("CV(κ_bd)", "cv_kappa_bd"),
        ]:
            vals = [r[key] for r in results]
            print(f"  {name:>20}  ", end="")
            for v in vals:
                print(f"{v:14.4e}  ", end="")
            if len(vals) >= 2 and abs(vals[0]) > 0 and abs(vals[-1]) > 0:
                lams = [r["lam"] for r in results]
                try:
                    exp = math.log(abs(vals[-1]) / abs(vals[0])) / math.log(lams[-1] / lams[0])
                    print(f"  {exp:8.2f}")
                except (ValueError, ZeroDivisionError):
                    print()
            else:
                print()

    print(f"\n=== FERTIG ===")
