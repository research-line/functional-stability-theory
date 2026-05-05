"""
c2ad_boundary_layer_test.py — C2AD: Near/Far Boundary-Layer Split

Verifiziert:
1. Shoaling-Parameter Σ_j = A_λ/√h_j^eff (C2ac)
2. J_near/J_far Split (C2ad.1)
3. Far-Sector: A_λ < ε_λ → |δ_m| < 1 (C2ad.2)
4. Near-Sector: Abel/PP-Bound (C2ad.3)
5. Combined: |δ_j| < 1 ∀ j ∈ J_reg (C2ad.4)

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


def boundary_layer_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*120}", flush=True)
    print(f"C2ad Boundary-Layer Split: λ={lam}, N={N}", flush=True)
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

    # --- A_λ = ‖R‖/α ---
    R_norm_sq = sum((c_arr[a] * (w_arr[a] - w_min))**2 for a in noncl_A)
    R_norm = np.sqrt(R_norm_sq)
    A_lam = R_norm / abs(alpha)

    # --- Interlacing: s_m in (u_{m-1}, u_m) ---
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

    # --- b_m coefficients for Abel bound ---
    b_coeffs = np.zeros(M + 1)
    b_coeffs[0] = A_cl
    for idx, a in enumerate(noncl_sorted):
        b_coeffs[idx + 1] = c_arr[a] * alpha_a[a] / alpha

    # --- Compute Shoaling + Abel for each matched mode ---
    results = []

    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        u_left = entry['u_left']
        u_right = entry['u_right']

        gap_L = s_m - u_left
        gap_R = u_right - s_m
        h_eff = gap_L * gap_R
        sigma_j = A_lam / np.sqrt(h_eff) if h_eff > 0 else float('inf')

        # δ_m actual
        d = s_m - w_min
        delta_actual = 0
        if abs(d) > 1e-30:
            delta_actual = -d * sum(c_arr[a] * alpha_a[a] / (alpha * (s_m - w_arr[a]))
                                    for a in noncl_A if abs(s_m - w_arr[a]) > 1e-30)

        # Abel/PP bound (C2ad.3)
        B_prev = sum(b_coeffs[n] for n in range(m))
        C_curr = sum(b_coeffs[n] for n in range(m, M + 1))
        abel_L = abs(B_prev) / gap_L if gap_L > 1e-50 else float('inf')
        abel_R = abs(C_curr) / gap_R if gap_R > 1e-50 else float('inf')
        abel_max = max(abel_L, abel_R)
        pp_bound = abs(d) * abel_max

        results.append({
            'm': m, 's': s_m, 'u_left': u_left, 'u_right': u_right,
            'gap_L': gap_L, 'gap_R': gap_R, 'h_eff': h_eff,
            'sigma': sigma_j, 'delta': delta_actual,
            'B_prev': B_prev, 'C_curr': C_curr,
            'abel_L': abel_L, 'abel_R': abel_R,
            'pp_bound': pp_bound,
        })

    # --- Near/Far split ---
    eps_lam = A_lam
    J_far = [r for r in results if r['h_eff'] > eps_lam**2]
    J_near = [r for r in results if r['h_eff'] <= eps_lam**2]

    print(f"\n--- Globale Parameter ---", flush=True)
    print(f"  A_λ = ‖R‖/α      = {A_lam:.6e}", flush=True)
    print(f"  A_cl              = {A_cl:.10f}", flush=True)
    print(f"  M (Off-Cluster)   = {M}", flush=True)
    print(f"  Matched modes     = {len(results)}", flush=True)
    print(f"  ε_λ = A_λ         = {eps_lam:.6e}", flush=True)
    print(f"  J_far             = {len(J_far)}", flush=True)
    print(f"  J_near            = {len(J_near)}", flush=True)

    # --- Far sector ---
    print(f"\n--- Far-Sector (C2ad.2): Shoaling Σ_j = A_λ/√h_eff ---", flush=True)
    far_ok = all(r['sigma'] < 1 for r in J_far)
    if J_far:
        max_sigma_far = max(r['sigma'] for r in J_far)
        print(f"  max Σ_j (far)     = {max_sigma_far:.6e}", flush=True)
        print(f"  Far bound < 1     = {'✓' if far_ok else '✗'}", flush=True)

    # --- Near sector ---
    print(f"\n--- Near-Sector (C2ad.3): Abel/PP Bound ---", flush=True)
    if J_near:
        print(f"\n{'m':>4}  {'s_m':>12}  {'h_eff':>10}  {'Σ_j':>10}  {'|δ|':>10}  "
              f"{'B_{m-1}':>10}  {'C_m':>10}  {'PP bnd':>10}  {'Σ<1':>5}  {'PP<1':>5}", flush=True)
        print('-' * 115, flush=True)
        for r in J_near:
            s_ok = '✓' if r['sigma'] < 1 else '✗'
            p_ok = '✓' if r['pp_bound'] < 1 else '✗'
            print(f"{r['m']:4d}  {r['s']:12.6e}  {r['h_eff']:10.2e}  {r['sigma']:10.4e}  "
                  f"{abs(r['delta']):10.4e}  {r['B_prev']:10.4e}  {r['C_curr']:10.4e}  "
                  f"{r['pp_bound']:10.4e}  {s_ok:>5}  {p_ok:>5}", flush=True)
    else:
        print(f"  (keine Near-Cluster-Moden)", flush=True)

    # --- Full Shoaling table (first/last + worst) ---
    print(f"\n--- Vollständige Shoaling-Tabelle (Auswahl) ---", flush=True)
    print(f"{'m':>4}  {'s_m':>12}  {'h_eff':>10}  {'Σ_j':>10}  {'|δ|':>10}  "
          f"{'PP bnd':>10}  {'Σ<1':>5}  {'PP<1':>5}  {'Sektor':>6}", flush=True)
    print('-' * 100, flush=True)

    for i, r in enumerate(results):
        s_ok = '✓' if r['sigma'] < 1 else '✗'
        p_ok = '✓' if r['pp_bound'] < 1 else '✗'
        sector = 'NEAR' if r['h_eff'] <= eps_lam**2 else 'FAR'
        if i < 10 or i >= len(results) - 5 or r['sigma'] >= 1 or r == max(results, key=lambda x: x['sigma']):
            print(f"{r['m']:4d}  {r['s']:12.6e}  {r['h_eff']:10.2e}  {r['sigma']:10.4e}  "
                  f"{abs(r['delta']):10.4e}  {r['pp_bound']:10.4e}  {s_ok:>5}  {p_ok:>5}  {sector:>6}", flush=True)
        elif i == 10:
            print("  ...", flush=True)

    # --- Summary ---
    print(f"\n--- Summary λ={lam} ---", flush=True)
    if results:
        max_sigma = max(r['sigma'] for r in results)
        max_delta = max(abs(r['delta']) for r in results)
        max_pp = max(r['pp_bound'] for r in results)
        n_sigma_ok = sum(1 for r in results if r['sigma'] < 1)
        n_pp_ok = sum(1 for r in results if r['pp_bound'] < 1)
        n_total = len(results)
        n_delta_ok = sum(1 for r in results if abs(r['delta']) < 1)

        # Combined: each mode passes if EITHER Σ_j < 1 OR pp_bound < 1
        n_combined = sum(1 for r in results if r['sigma'] < 1 or r['pp_bound'] < 1)

        worst_sigma = max(results, key=lambda r: r['sigma'])
        worst_pp = max(results, key=lambda r: r['pp_bound'])

        print(f"A_λ = ‖R‖/α         = {A_lam:.6e}", flush=True)
        print(f"max |δ|              = {max_delta:.6e}", flush=True)
        print(f"|δ| < 1              = {n_delta_ok}/{n_total}", flush=True)
        print(f"Shoaling Σ_j < 1     = {n_sigma_ok}/{n_total}  (max: {max_sigma:.4e})", flush=True)
        print(f"Abel PP bound < 1    = {n_pp_ok}/{n_total}  (max: {max_pp:.4e})", flush=True)
        print(f"COMBINED (Σ∨PP) < 1  = {n_combined}/{n_total}  ← C2ad.4", flush=True)
        print(f"Worst Σ: mode {worst_sigma['m']}, Σ={worst_sigma['sigma']:.4e}, "
              f"h_eff={worst_sigma['h_eff']:.2e}", flush=True)
        print(f"Worst PP: mode {worst_pp['m']}, PP={worst_pp['pp_bound']:.4e}", flush=True)

    return {
        'A_lam': A_lam, 'results': results,
        'J_near': J_near, 'J_far': J_far,
    }


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = boundary_layer_test(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    # Cross-λ summary
    print(f"\n{'='*120}", flush=True)
    print("Cross-λ: C2ad Boundary-Layer Summary", flush=True)
    print(f"{'='*120}", flush=True)
    print(f"{'λ':>4}  {'A_λ':>10}  {'|J_far|':>7}  {'|J_near|':>8}  {'max Σ':>10}  "
          f"{'max PP':>10}  {'max|δ|':>10}  {'Σ<1':>6}  {'PP<1':>6}  {'COMB':>6}", flush=True)
    for lam in sorted(all_results):
        r = all_results[lam]
        if r['results']:
            max_s = max(x['sigma'] for x in r['results'])
            max_p = max(x['pp_bound'] for x in r['results'])
            max_d = max(abs(x['delta']) for x in r['results'])
            n_s = sum(1 for x in r['results'] if x['sigma'] < 1)
            n_p = sum(1 for x in r['results'] if x['pp_bound'] < 1)
            n_c = sum(1 for x in r['results'] if x['sigma'] < 1 or x['pp_bound'] < 1)
            n_t = len(r['results'])
            print(f"{lam:4.0f}  {r['A_lam']:10.4e}  {len(r['J_far']):7d}  {len(r['J_near']):8d}  "
                  f"{max_s:10.4e}  {max_p:10.4e}  {max_d:10.4e}  "
                  f"{n_s}/{n_t:>3}  {n_p}/{n_t:>3}  {n_c}/{n_t:>3}", flush=True)

    print("\nDone.", flush=True)
