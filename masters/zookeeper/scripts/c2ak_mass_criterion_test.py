"""
c2ak_mass_criterion_test.py — Off-Cluster Total Mass Criterion

Testet C2ak.1: |δ_m| ≤ M_off · d_m / g_min
  und C2ak.4: M_off ≤ ‖R‖ / (α γ₁)
  und C2ak.5: |δ_m| ≤ (‖R‖/(α γ₁)) · d_m / g_min

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


def mass_criterion_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*130}", flush=True)
    print(f"C2ak Mass Criterion: λ={lam}, N={N}", flush=True)
    print(f"{'='*130}", flush=True)

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
    alpha_cl = sum(c_arr[a] * alpha_a[a] for a in cl_A)

    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    u_poles = [w_min]
    for a in noncl_sorted:
        u_poles.append(w_arr[a])
    u_poles = np.array(u_poles)
    M = len(u_poles) - 1

    b_off = np.zeros(M)
    for idx, a in enumerate(noncl_sorted):
        b_off[idx] = c_arr[a] * alpha_a[a] / alpha

    M_off = np.sum(b_off)
    M_off_abs = np.sum(np.abs(b_off))
    n_neg_b = np.sum(b_off < -1e-30)

    gamma1 = u_poles[1] - u_poles[0] if M > 0 else float('inf')

    R_components = np.zeros(len(w_arr))
    for a in range(len(w_arr)):
        R_components[a] = c_arr[a] * (w_arr[a] - w_min)
    R_norm = np.sqrt(np.sum(R_components**2))

    alpha_off = sum(c_arr[a] * alpha_a[a] for a in noncl_A)
    M_off_identity = alpha_off / alpha

    M_off_upper = R_norm / (alpha * gamma1)

    print(f"\n--- Globale Größen ---", flush=True)
    print(f"  α           = {alpha:.10f}", flush=True)
    print(f"  α_cl        = {alpha_cl:.10f}", flush=True)
    print(f"  α_off       = {alpha_off:.6e}", flush=True)
    print(f"  A_cl        = {A_cl:.10f}", flush=True)
    print(f"  M_off       = {M_off:.6e}  (signed)", flush=True)
    print(f"  |M_off|     = {M_off_abs:.6e}  (absolute)", flush=True)
    print(f"  α_off/α     = {M_off_identity:.6e}  (C2ak.3 Identity)", flush=True)
    print(f"  1-α_cl/α    = {1-alpha_cl/alpha:.6e}  (alternative)", flush=True)
    print(f"  b_n < 0     = {n_neg_b}/{M}  (Positivitätsverletzung)", flush=True)
    print(f"  γ₁          = {gamma1:.6e}  (erster Off-Cluster-Gap)", flush=True)
    print(f"  ‖R‖         = {R_norm:.6e}", flush=True)
    print(f"  ‖R‖/(αγ₁)   = {M_off_upper:.6e}  (C2ak.4 Upper)", flush=True)
    print(f"  M_off/Upper = {M_off_abs/M_off_upper:.4f}  (Tightness von C2ak.4)", flush=True)

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

    print(f"\n{'m':>4}  {'d_m':>10}  {'g_min':>10}  {'d/g_min':>10}  "
          f"{'M·d/g':>10}  {'R/(αγ)·d/g':>12}  {'|δ_m|':>10}  "
          f"{'C2ak.1':>7}  {'C2ak.5':>7}  {'<1':>4}", flush=True)
    print('-' * 120, flush=True)

    results = []
    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        d_m = s_m - w_min
        gap_L = s_m - entry['u_left']
        gap_R = entry['u_right'] - s_m
        g_min = min(gap_L, gap_R)

        delta_m = 0
        if abs(d_m) > 1e-30:
            delta_m = -d_m * sum(c_arr[a] * alpha_a[a] / (alpha * (s_m - w_arr[a]))
                                 for a in noncl_A if abs(s_m - w_arr[a]) > 1e-30)

        geo_ratio = d_m / g_min if g_min > 1e-50 else float('inf')
        c2ak1 = M_off_abs * geo_ratio
        c2ak5 = M_off_upper * geo_ratio

        ok_ak1 = '✓' if abs(delta_m) <= c2ak1 * 1.001 else '✗'
        ok_ak5 = '✓' if abs(delta_m) <= c2ak5 * 1.001 else '✗'
        lt1 = '✓' if c2ak1 < 1 else '✗'

        results.append({
            'm': m, 'd': d_m, 'g_min': g_min, 'geo': geo_ratio,
            'c2ak1': c2ak1, 'c2ak5': c2ak5, 'delta': delta_m,
        })

        if m <= 8 or m >= M - 3 or c2ak1 >= 0.01:
            print(f"{m:4d}  {d_m:10.4e}  {g_min:10.4e}  {geo_ratio:10.4e}  "
                  f"{c2ak1:10.4e}  {c2ak5:12.4e}  {abs(delta_m):10.4e}  "
                  f"{ok_ak1:>7}  {ok_ak5:>7}  {lt1:>4}", flush=True)
        elif m == 9:
            print("  ...", flush=True)

    print(f"\n--- Summary λ={lam} ---", flush=True)
    n_total = len(results)
    n_ak1_ok = sum(1 for r in results if abs(r['delta']) <= r['c2ak1'] * 1.001)
    n_ak5_ok = sum(1 for r in results if abs(r['delta']) <= r['c2ak5'] * 1.001)
    n_lt1 = sum(1 for r in results if r['c2ak1'] < 1)
    max_c2ak1 = max(r['c2ak1'] for r in results)
    max_c2ak5 = max(r['c2ak5'] for r in results)
    max_geo = max(r['geo'] for r in results)
    max_delta = max(abs(r['delta']) for r in results)
    worst = max(results, key=lambda r: r['c2ak1'])

    print(f"Total modes                = {n_total}", flush=True)
    print(f"|δ| ≤ M_off·d/g_min (C2ak.1) = {n_ak1_ok}/{n_total}", flush=True)
    print(f"|δ| ≤ ‖R‖/(αγ)·d/g (C2ak.5) = {n_ak5_ok}/{n_total}", flush=True)
    print(f"M_off·d/g_min < 1           = {n_lt1}/{n_total}  ← ENTSCHEIDEND", flush=True)
    print(f"max M_off·d/g_min           = {max_c2ak1:.6e}", flush=True)
    print(f"max ‖R‖/(αγ)·d/g            = {max_c2ak5:.6e}", flush=True)
    print(f"max d/g_min                 = {max_geo:.6e}", flush=True)
    print(f"max |δ|                     = {max_delta:.6e}", flush=True)
    print(f"Worst mode: m={worst['m']}, M_off·d/g={worst['c2ak1']:.4e}", flush=True)

    return {
        'n_total': n_total, 'n_lt1': n_lt1,
        'max_c2ak1': max_c2ak1, 'max_c2ak5': max_c2ak5,
        'max_geo': max_geo, 'max_delta': max_delta,
        'M_off': M_off_abs, 'M_off_upper': M_off_upper,
        'gamma1': gamma1, 'R_norm': R_norm,
        'results': results,
    }


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = mass_criterion_test(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    print(f"\n{'='*130}", flush=True)
    print("Cross-λ: C2ak Total Off-Cluster Mass Criterion", flush=True)
    print(f"{'='*130}", flush=True)
    print(f"{'λ':>4}  {'n':>4}  {'M_off':>10}  {'‖R‖/(αγ)':>10}  "
          f"{'max d/g':>10}  {'max C2ak.1':>10}  {'max C2ak.5':>10}  "
          f"{'max|δ|':>10}  {'<1':>6}  {'SCHLUSS':>10}", flush=True)
    for lam in sorted(all_results):
        r = all_results[lam]
        ok = '✓ C2o.1' if r['n_lt1'] == r['n_total'] else '✗ OFFEN'
        print(f"{lam:4.0f}  {r['n_total']:4d}  {r['M_off']:10.4e}  {r['M_off_upper']:10.4e}  "
              f"{r['max_geo']:10.4e}  {r['max_c2ak1']:10.4e}  {r['max_c2ak5']:10.4e}  "
              f"{r['max_delta']:10.4e}  {r['n_lt1']}/{r['n_total']:>3}  {ok:>10}", flush=True)

    print("\nDone.", flush=True)
