"""
c2ae_information_theory_test.py — Informationstheoretische Diagnostik

Berechnet:
1. Off-Cluster-Massenverteilung p_a = b_a / Σ b_n
2. Shannon-Entropie H_off = -Σ p_a log p_a
3. Effektives Spektrum N_eff = e^{H_off}
4. Lokale Hazard-Rate: C_{m+1} / (u_m - s_m)
5. Fisher-like Sensitivität: (m_A^off)'(s_j)
6. Shoaling × Fisher Korrelation

Motivation: GPT-Vorschlag (Session 27b++). Wenn N_eff klein bleibt
obwohl n_off wächst → formale Massenkonzentration = Zwischenlemma für C2ad.

Autor: LG (Opus 4.6, Session 27b++)
Datum: 2026-04-20
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from c2q_first_shell_dominance import build_operators

DPS = int(os.environ.get("DPS", 50))
CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
    {"lam": 7.0, "N": 77},
]


def info_theory_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*120}", flush=True)
    print(f"C2ae Information Theory Diagnostics: λ={lam}, N={N}", flush=True)
    print(f"{'='*120}", flush=True)

    t0 = time.time()
    print("Building operators...", flush=True)
    ops = build_operators(lam, N)
    print(f"Build complete in {time.time()-t0:.1f}s", flush=True)

    alpha = ops['alpha']
    w_arr = ops['w_arr']
    c_arr = ops['c_arr']
    alpha_a = ops['alpha_a']
    noncl_A = ops['noncl_A']
    cl_A = ops['cl_A']
    t_vals = ops['t_vals']
    f_arr = ops['f_arr']
    j_reg = ops['j_reg']
    w_min = ops['w_min']

    A_cl = sum(alpha_a[a]**2 for a in cl_A)
    n_off = len(noncl_A)

    # --- Compressed pole set ---
    u_poles = [w_min]
    beta_poles = [A_cl]
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    for a in noncl_sorted:
        u_poles.append(w_arr[a])
        beta_poles.append(alpha_a[a]**2)
    u_poles = np.array(u_poles)
    beta_poles = np.array(beta_poles)
    M = len(u_poles) - 1

    # A_λ
    R_norm_sq = sum((c_arr[a] * (w_arr[a] - w_min))**2 for a in noncl_A)
    R_norm = np.sqrt(R_norm_sq)
    A_lam = R_norm / abs(alpha)

    # --- 1. Off-Cluster Mass Distribution ---
    b_raw = np.array([c_arr[a] * alpha_a[a] / alpha for a in noncl_sorted])
    b_abs = np.abs(b_raw)
    b_total = np.sum(b_abs)

    p_a = b_abs / b_total if b_total > 0 else np.zeros_like(b_abs)

    print(f"\n--- 1. Off-Cluster Mass Distribution ---", flush=True)
    print(f"  n_off (Off-Cluster)     = {n_off}", flush=True)
    print(f"  Σ|b_a|                  = {b_total:.6e}", flush=True)
    print(f"  Σ b_a (signed)          = {np.sum(b_raw):.6e}", flush=True)

    print(f"\n  Top-10 masses (by |b_a|):", flush=True)
    idx_sorted = np.argsort(-b_abs)
    cumul = 0
    for rank, i in enumerate(idx_sorted[:10]):
        cumul += p_a[i]
        a = noncl_sorted[i]
        print(f"    rank {rank+1:2d}: a={a:3d}, w_a-w_min={w_arr[a]-w_min:10.4e}, "
              f"|b_a|={b_abs[i]:10.4e}, p_a={p_a[i]:8.4f}, cumul={cumul:8.4f}", flush=True)

    # --- 2. Shannon Entropy and N_eff ---
    H_off = 0
    for p in p_a:
        if p > 1e-50:
            H_off -= p * np.log(p)
    N_eff = np.exp(H_off)

    print(f"\n--- 2. Shannon Entropy ---", flush=True)
    print(f"  H_off                   = {H_off:.4f} nats", flush=True)
    print(f"  H_off / log(n_off)      = {H_off / np.log(n_off):.4f}  (1.0 = uniform)", flush=True)
    print(f"  N_eff = e^H             = {N_eff:.2f}", flush=True)
    print(f"  N_eff / n_off           = {N_eff / n_off:.4f}  (fraction active)", flush=True)
    print(f"  n_off                   = {n_off}", flush=True)

    # --- 3. Cumulative distribution (CDF) ---
    cdf = np.cumsum(p_a[idx_sorted])
    pct50 = np.searchsorted(cdf, 0.5) + 1
    pct90 = np.searchsorted(cdf, 0.9) + 1
    pct99 = np.searchsorted(cdf, 0.99) + 1

    print(f"\n--- 3. Cumulative Mass Concentration ---", flush=True)
    print(f"  50% mass in top         = {pct50} modes ({pct50/n_off*100:.1f}%)", flush=True)
    print(f"  90% mass in top         = {pct90} modes ({pct90/n_off*100:.1f}%)", flush=True)
    print(f"  99% mass in top         = {pct99} modes ({pct99/n_off*100:.1f}%)", flush=True)

    # --- 4. Interlacing + Hazard Rate ---
    t_all_reg = sorted([t_vals[j] for j in j_reg])
    s_zeros = []
    t_matched = set()

    for m in range(1, M + 1):
        u_left = u_poles[m - 1]
        u_right = u_poles[m]
        candidates = [t for t in t_all_reg if u_left < t < u_right and t not in t_matched]
        if candidates:
            s = candidates[0]
            s_zeros.append({'m': m, 's': s, 'u_left': u_left, 'u_right': u_right})
            t_matched.add(s)

    print(f"\n--- 4. Hazard Rate: C_{{m+1}} / gap_R ---", flush=True)
    print(f"{'m':>4}  {'s_m':>12}  {'gap_R':>10}  {'C_m':>10}  {'C_{m+1}':>10}  "
          f"{'hazard':>10}  {'h<1':>5}  {'d·hazard':>10}  {'dh<1':>5}", flush=True)
    print('-' * 100, flush=True)

    hazard_results = []
    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        gap_R = entry['u_right'] - s_m
        d = s_m - w_min

        C_m = sum(b_abs[i] for i in range(m - 1, M))
        C_m1 = sum(b_abs[i] for i in range(m, M))

        hazard = C_m1 / gap_R if gap_R > 1e-50 else float('inf')
        d_hazard = d * hazard

        hazard_results.append({
            'm': m, 's': s_m, 'gap_R': gap_R, 'd': d,
            'C_m': C_m, 'C_m1': C_m1, 'hazard': hazard, 'd_hazard': d_hazard,
        })

        if m <= 8 or m >= M - 3 or d_hazard >= 0.5:
            h_ok = '✓' if hazard < 1 else '✗'
            dh_ok = '✓' if d_hazard < 1 else '✗'
            print(f"{m:4d}  {s_m:12.6e}  {gap_R:10.2e}  {C_m:10.4e}  {C_m1:10.4e}  "
                  f"{hazard:10.4e}  {h_ok:>5}  {d_hazard:10.4e}  {dh_ok:>5}", flush=True)
        elif m == 9:
            print("  ...", flush=True)

    # --- 5. Fisher-like Sensitivity ---
    print(f"\n--- 5. Fisher-like Sensitivity: (m_A^off)'(s_j) ---", flush=True)

    fisher_results = []
    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        gap_L = s_m - entry['u_left']
        gap_R = entry['u_right'] - s_m
        h_eff = gap_L * gap_R

        fisher = sum(alpha_a[a]**2 / (s_m - w_arr[a])**2
                     for a in noncl_A if abs(s_m - w_arr[a]) > 1e-30)

        local_bound = 1.0 / h_eff if h_eff > 0 else float('inf')
        ratio = fisher / local_bound if local_bound > 0 and local_bound < 1e30 else 0

        fisher_results.append({
            'm': m, 'fisher': fisher, 'local_bound': local_bound,
            'h_eff': h_eff, 'ratio': ratio,
        })

    if fisher_results:
        print(f"  min Fisher              = {min(r['fisher'] for r in fisher_results):.4e}", flush=True)
        print(f"  max Fisher              = {max(r['fisher'] for r in fisher_results):.4e}", flush=True)
        print(f"  Fisher/local_bound ratio: min={min(r['ratio'] for r in fisher_results):.6f}, "
              f"max={max(r['ratio'] for r in fisher_results):.6f}", flush=True)
        print(f"  (ratio ≤ 1 by C2aa.3 — Interlacing contracts non-local factors)", flush=True)

    # --- 6. Summary ---
    print(f"\n--- Summary λ={lam} ---", flush=True)
    print(f"n_off              = {n_off}", flush=True)
    print(f"H_off              = {H_off:.4f} nats", flush=True)
    print(f"N_eff              = {N_eff:.2f}  ({N_eff/n_off*100:.1f}% of n_off)", flush=True)
    print(f"50/90/99% mass     = {pct50}/{pct90}/{pct99} modes", flush=True)
    print(f"A_λ                = {A_lam:.6e}", flush=True)

    if hazard_results:
        max_hz = max(r['hazard'] for r in hazard_results)
        max_dhz = max(r['d_hazard'] for r in hazard_results)
        n_hz_ok = sum(1 for r in hazard_results if r['d_hazard'] < 1)
        print(f"max hazard         = {max_hz:.4e}", flush=True)
        print(f"max d·hazard       = {max_dhz:.4e}", flush=True)
        print(f"d·hazard < 1       = {n_hz_ok}/{len(hazard_results)}", flush=True)

    return {
        'H_off': H_off, 'N_eff': N_eff, 'n_off': n_off,
        'A_lam': A_lam, 'pct50': pct50, 'pct90': pct90,
        'hazard_results': hazard_results,
    }


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = info_theory_test(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    # Cross-λ summary
    print(f"\n{'='*120}", flush=True)
    print("Cross-λ: Information Theory Summary", flush=True)
    print(f"{'='*120}", flush=True)
    print(f"{'λ':>4}  {'n_off':>6}  {'H_off':>8}  {'N_eff':>8}  {'N/n':>6}  "
          f"{'50%':>4}  {'90%':>4}  {'99%':>4}  {'A_λ':>10}  "
          f"{'max d·hz':>10}  {'dh<1':>6}", flush=True)
    for lam in sorted(all_results):
        r = all_results[lam]
        max_dhz = max(h['d_hazard'] for h in r['hazard_results']) if r['hazard_results'] else 0
        n_ok = sum(1 for h in r['hazard_results'] if h['d_hazard'] < 1)
        n_t = len(r['hazard_results'])
        print(f"{lam:4.0f}  {r['n_off']:6d}  {r['H_off']:8.4f}  {r['N_eff']:8.2f}  "
              f"{r['N_eff']/r['n_off']:6.3f}  {r['pct50']:4d}  {r['pct90']:4d}  "
              f"{'—':>4}  {r['A_lam']:10.4e}  {max_dhz:10.4e}  {n_ok}/{n_t:>3}", flush=True)

    # Key question: does N_eff stay bounded?
    lams = sorted(all_results.keys())
    if len(lams) >= 2:
        n_effs = [all_results[l]['N_eff'] for l in lams]
        n_offs = [all_results[l]['n_off'] for l in lams]
        print(f"\n--- Key Question: N_eff Scaling ---", flush=True)
        print(f"  N_eff values: {', '.join(f'{x:.2f}' for x in n_effs)}", flush=True)
        print(f"  n_off values: {', '.join(f'{x}' for x in n_offs)}", flush=True)
        if n_effs[-1] > 0 and n_effs[0] > 0:
            growth = np.log(n_effs[-1] / n_effs[0]) / np.log(lams[-1] / lams[0])
            print(f"  N_eff growth exponent (λ^α): α ≈ {growth:.3f}", flush=True)
            print(f"  n_off growth exponent:       α ≈ {np.log(n_offs[-1]/n_offs[0])/np.log(lams[-1]/lams[0]):.3f}", flush=True)
            if growth < 0.5:
                print(f"  → N_eff wächst LANGSAMER als √λ — STARKE Massenkonzentration!", flush=True)
            elif growth < 1.0:
                print(f"  → N_eff wächst sublinear — moderate Konzentration", flush=True)
            else:
                print(f"  → N_eff wächst linear oder schneller — schwache Konzentration", flush=True)

    print("\nDone.", flush=True)
