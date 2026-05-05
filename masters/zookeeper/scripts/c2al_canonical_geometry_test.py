"""
c2al_canonical_geometry_test.py — Canonical Local Geometry + Gap Balance Test

Testet C2al.1: d_m/g_R ≤ T_m^β / β₀
  und C2al.2: Θ^off ≤ ((1-β₀)/β₀) · M_off  (rechter Kanal GESCHLOSSEN)
  und C2am.1: β_{m-1}·g_R/g_L + β_m·g_L/g_R ≤ 1  (Slope-Inequality)
  und C2am.3: d_m/g_L ≤ T_m^β / (β₀ · β_{m-1})  (SCHÄRFER als C2al.3)
  und C2am.4: Λ^off ≤ M_off · T_m^β / (β₀ · β_{m-1})

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


def canonical_geometry_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*140}", flush=True)
    print(f"C2al Canonical Geometry: λ={lam}, N={N}", flush=True)
    print(f"{'='*140}", flush=True)

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
    j_reg = ops['j_reg']
    w_min = ops['w_min']

    A_cl = sum(alpha_a[a]**2 for a in cl_A)
    beta_0 = A_cl

    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    u_poles = [w_min]
    beta_arr = [beta_0]
    for a in noncl_sorted:
        u_poles.append(w_arr[a])
        beta_arr.append(alpha_a[a]**2)
    u_poles = np.array(u_poles)
    beta_arr = np.array(beta_arr)
    M = len(u_poles) - 1

    b_off = np.zeros(M)
    for idx, a in enumerate(noncl_sorted):
        b_off[idx] = c_arr[a] * alpha_a[a] / alpha
    M_off = np.sum(np.abs(b_off))

    T_beta = np.zeros(M + 1)
    T_beta[M] = beta_arr[M]
    for n in range(M - 1, -1, -1):
        T_beta[n] = T_beta[n + 1] + beta_arr[n]

    global_right_bound = (1 - beta_0) / beta_0 * M_off

    print(f"\n--- Globale Größen ---", flush=True)
    print(f"  β₀ = A_cl      = {beta_0:.10f}", flush=True)
    print(f"  1-β₀            = {1-beta_0:.6e}", flush=True)
    print(f"  M_off            = {M_off:.6e}", flush=True)
    print(f"  (1-β₀)/β₀·M_off = {global_right_bound:.6e}  (globaler Θ^off Bound)", flush=True)
    print(f"  Σβ_n             = {sum(beta_arr):.10f}  (sollte ≈ 1)", flush=True)

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

    print(f"\n{'m':>4}  {'d/gR':>10}  {'Tβ/β₀':>10}  {'al.1':>5}  "
          f"{'slope':>8}  {'am.1':>5}  "
          f"{'Tβ/β_{m-1}':>10}  {'Tβ/(β₀β)':>12}  {'am.3':>5}  "
          f"{'Λ_am4':>10}  {'Λ^off':>10}  {'<1':>4}", flush=True)
    print('-' * 130, flush=True)

    results = []
    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        d_m = s_m - w_min
        gap_L = s_m - entry['u_left']
        gap_R = entry['u_right'] - s_m

        T_m = T_beta[m]
        beta_m1 = beta_arr[m - 1] if m >= 1 else 0

        d_gR = d_m / gap_R if gap_R > 1e-50 else float('inf')
        d_gL = d_m / gap_L if gap_L > 1e-50 else float('inf')
        Tb_b0 = T_m / beta_0
        Tb_b0b = T_m / (beta_0 * beta_m1) if beta_m1 > 1e-50 else float('inf')
        Tb_ratio = T_m / beta_m1 if beta_m1 > 1e-50 else float('inf')

        beta_m_val = beta_arr[m] if m <= M else 0
        slope_lhs = beta_m1 * gap_R / gap_L + beta_m_val * gap_L / gap_R if gap_L > 1e-50 and gap_R > 1e-50 else float('inf')

        B_off_prev = np.sum(b_off[:m-1]) if m > 1 else 0
        C_off_curr = np.sum(b_off[m-1:])
        Lambda_off = abs(d_m) * abs(B_off_prev) / gap_L if gap_L > 1e-50 else float('inf')
        Theta_off = abs(d_m) * abs(C_off_curr) / gap_R if gap_R > 1e-50 else float('inf')

        Lambda_am4 = M_off * Tb_b0b

        ok_al1 = '✓' if d_gR <= Tb_b0 * 1.001 else '✗'
        ok_am1 = '✓' if slope_lhs <= 1.001 else '✗'
        ok_am3 = '✓' if d_gL <= Tb_b0b * 1.001 else '✗'
        lt1 = '✓' if max(Lambda_am4, global_right_bound) < 1 else '✗'

        results.append({
            'm': m, 'd_gR': d_gR, 'd_gL': d_gL,
            'Tb_b0': Tb_b0, 'Tb_b0b': Tb_b0b, 'Tb_ratio': Tb_ratio,
            'Lambda_off': Lambda_off, 'Theta_off': Theta_off,
            'Lambda_am4': Lambda_am4, 'slope_lhs': slope_lhs,
            'ok_al1': ok_al1, 'ok_am1': ok_am1, 'ok_am3': ok_am3,
            'beta_m1': beta_m1,
        })

        if m <= 8 or m >= M - 3 or Lambda_am4 >= 0.01 or Tb_ratio >= 10:
            print(f"{m:4d}  {d_gR:10.4e}  {Tb_b0:10.4e}  {ok_al1:>5}  "
                  f"{slope_lhs:8.4f}  {ok_am1:>5}  "
                  f"{Tb_ratio:10.4e}  {Tb_b0b:12.4e}  {ok_am3:>5}  "
                  f"{Lambda_am4:10.4e}  {Lambda_off:10.4e}  {lt1:>4}", flush=True)
        elif m == 9:
            print("  ...", flush=True)

    print(f"\n--- Summary λ={lam} ---", flush=True)
    n_total = len(results)
    n_al1 = sum(1 for r in results if r['ok_al1'] == '✓')
    n_am1 = sum(1 for r in results if r['ok_am1'] == '✓')
    n_am3 = sum(1 for r in results if r['ok_am3'] == '✓')
    max_Tb_b0 = max(r['Tb_b0'] for r in results)
    max_Tb_b0b = max(r['Tb_b0b'] for r in results)
    max_slope = max(r['slope_lhs'] for r in results)
    max_Lambda_am4 = max(r['Lambda_am4'] for r in results)
    max_Tb_ratio = max(r['Tb_b0b'] * beta_0 for r in results)
    n_full_lt1 = sum(1 for r in results if max(r['Lambda_am4'], global_right_bound) < 1)

    print(f"Total modes                          = {n_total}", flush=True)
    print(f"C2al.1 (d/gR ≤ Tβ/β₀)              = {n_al1}/{n_total}", flush=True)
    print(f"C2al.2 (Θ^off ≤ (1-β₀)/β₀·M_off)    = {global_right_bound:.6e}  (global, uniform)", flush=True)
    print(f"C2am.1 (slope ≤ 1)                   = {n_am1}/{n_total}  (max slope: {max_slope:.6f})", flush=True)
    print(f"C2am.3 (d/gL ≤ Tβ/(β₀β))            = {n_am3}/{n_total}", flush=True)
    print(f"max Tβ/β₀                            = {max_Tb_b0:.6e}", flush=True)
    print(f"max Tβ/(β₀β_{{m-1}})                  = {max_Tb_b0b:.6e}", flush=True)
    print(f"max Tβ/β_{{m-1}}                       = {max_Tb_ratio:.6e}  ← KRITISCHES VERHÄLTNIS", flush=True)
    print(f"max M_off·Tβ/(β₀β) (Λ^off am.4)     = {max_Lambda_am4:.6e}", flush=True)
    print(f"max(Λ_am4, global Θ) < 1             = {n_full_lt1}/{n_total}  ← ENTSCHEIDEND", flush=True)

    return {
        'n_total': n_total, 'n_al1': n_al1, 'n_am1': n_am1, 'n_am3': n_am3,
        'global_right': global_right_bound,
        'max_Tb_b0': max_Tb_b0, 'max_Tb_b0b': max_Tb_b0b,
        'max_slope': max_slope, 'max_Tb_ratio': max_Tb_ratio,
        'max_Lambda_am4': max_Lambda_am4,
        'n_full_lt1': n_full_lt1,
        'M_off': M_off, 'beta_0': beta_0,
    }


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = canonical_geometry_test(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    print(f"\n{'='*140}", flush=True)
    print("Cross-λ: C2al Canonical Geometry Summary", flush=True)
    print(f"{'='*140}", flush=True)
    print(f"{'λ':>4}  {'n':>4}  {'β₀':>8}  {'M_off':>10}  "
          f"{'gl.Θ':>10}  {'max Tβ/β₀':>10}  {'max Tβ/(β₀β)':>14}  {'max Tβ/β':>10}  "
          f"{'max Λ_am4':>10}  {'al.1':>5}  {'am.1':>5}  {'am.3':>5}  {'<1':>6}  {'SCHLUSS':>10}", flush=True)
    for lam in sorted(all_results):
        r = all_results[lam]
        ok = '✓ C2o.1' if r['n_full_lt1'] == r['n_total'] else '✗ OFFEN'
        print(f"{lam:4.0f}  {r['n_total']:4d}  {r['beta_0']:8.4f}  {r['M_off']:10.4e}  "
              f"{r['global_right']:10.4e}  {r['max_Tb_b0']:10.4e}  {r['max_Tb_b0b']:14.4e}  {r['max_Tb_ratio']:10.4e}  "
              f"{r['max_Lambda_am4']:10.4e}  {r['n_al1']:>3}/{r['n_total']}  {r['n_am1']:>3}/{r['n_total']}  "
              f"{r['n_am3']:>3}/{r['n_total']}  "
              f"{r['n_full_lt1']}/{r['n_total']:>3}  {ok:>10}", flush=True)

    print("\nDone.", flush=True)
