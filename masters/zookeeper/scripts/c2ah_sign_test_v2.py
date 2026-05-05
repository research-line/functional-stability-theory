"""
c2ah_sign_test_v2.py — Korrigierter Abel-Bound: NUR Off-Cluster-Massen

SCHLÜSSELKORREKTUR: δ_m summiert NUR über Off-Cluster-Pole (a ∉ Cl).
Die Abel-Schranken Λ^off, Θ^off müssen daher ebenfalls NUR Off-Cluster-
Massen b_n (n ≥ 1) enthalten, NICHT die Cluster-Masse b_0 = A_cl.

Wenn max(Λ^off, Θ^off) < 1 auf allen Moden → |δ_m| < 1 → C2o.1.
Keine Vorzeichenfrage, kein Near/Far-Split nötig.

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


def offcluster_abel_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*130}", flush=True)
    print(f"C2ah Off-Cluster Abel Bound: λ={lam}, N={N}", flush=True)
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

    # Compressed pole set
    u_poles = [w_min]
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    for a in noncl_sorted:
        u_poles.append(w_arr[a])
    u_poles = np.array(u_poles)
    M = len(u_poles) - 1

    # Off-cluster b_n (n = 1..M, NO cluster mass b_0)
    b_off = np.zeros(M)
    for idx, a in enumerate(noncl_sorted):
        b_off[idx] = c_arr[a] * alpha_a[a] / alpha

    total_off_mass = np.sum(np.abs(b_off))
    total_off_signed = np.sum(b_off)

    print(f"\n--- Off-Cluster Mass ---", flush=True)
    print(f"  A_cl (NICHT in Abel) = {A_cl:.10f}", flush=True)
    print(f"  M_off = Σ|b_n|      = {total_off_mass:.6e}", flush=True)
    print(f"  M_off_signed         = {total_off_signed:.6e}", flush=True)
    print(f"  M_off/A_cl           = {total_off_mass/A_cl:.6e}  (Verhältnis)", flush=True)

    # Interlacing
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

    print(f"\n{'m':>4}  {'s_m':>12}  {'d_m':>10}  {'δ_m':>12}  "
          f"{'Λ^off':>10}  {'Θ^off':>10}  {'max':>10}  {'|δ|<max':>8}  {'max<1':>6}", flush=True)
    print('-' * 120, flush=True)

    results = []

    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        u_left = entry['u_left']
        u_right = entry['u_right']
        d_m = s_m - w_min
        gap_L = s_m - u_left
        gap_R = u_right - s_m

        # δ_m (off-cluster only, same as before)
        delta_m = 0
        if abs(d_m) > 1e-30:
            delta_m = -d_m * sum(c_arr[a] * alpha_a[a] / (alpha * (s_m - w_arr[a]))
                                 for a in noncl_A if abs(s_m - w_arr[a]) > 1e-30)

        # OFF-CLUSTER Abel bounds (n = 1..M, NO b_0 = A_cl)
        # m is 1-indexed: poles are u_1..u_M, b_off is 0-indexed: b_off[0..M-1]
        # For mode m: left poles are u_1..u_{m-1} → b_off[0..m-2]
        #             right poles are u_m..u_M → b_off[m-1..M-1]
        B_off_prev = np.sum(b_off[:m-1]) if m > 1 else 0  # off-cluster left mass
        C_off_curr = np.sum(b_off[m-1:])  # off-cluster right mass (including u_m)

        Lambda_off = abs(d_m) * abs(B_off_prev) / gap_L if gap_L > 1e-50 else float('inf')
        Theta_off = abs(d_m) * abs(C_off_curr) / gap_R if gap_R > 1e-50 else float('inf')
        max_LT = max(Lambda_off, Theta_off)

        bound_ok = '✓' if abs(delta_m) <= max_LT * 1.001 else '✗'
        max_lt1 = '✓' if max_LT < 1 else '✗'

        results.append({
            'm': m, 's': s_m, 'd': d_m, 'delta': delta_m,
            'Lambda_off': Lambda_off, 'Theta_off': Theta_off,
            'max_LT': max_LT,
        })

        if m <= 8 or m >= M - 3 or max_LT >= 0.1:
            print(f"{m:4d}  {s_m:12.6e}  {d_m:10.4e}  {delta_m:12.4e}  "
                  f"{Lambda_off:10.4e}  {Theta_off:10.4e}  {max_LT:10.4e}  {bound_ok:>8}  {max_lt1:>6}", flush=True)
        elif m == 9:
            print("  ...", flush=True)

    # Summary
    print(f"\n--- Summary λ={lam} ---", flush=True)
    n_total = len(results)
    n_bound_ok = sum(1 for r in results if abs(r['delta']) <= r['max_LT'] * 1.001)
    n_max_lt1 = sum(1 for r in results if r['max_LT'] < 1)
    max_max = max(r['max_LT'] for r in results)
    max_delta = max(abs(r['delta']) for r in results)

    worst = max(results, key=lambda r: r['max_LT'])
    tightest = max(results, key=lambda r: abs(r['delta'])/r['max_LT'] if r['max_LT'] > 0 else 0)
    tight_ratio = abs(tightest['delta'])/tightest['max_LT'] if tightest['max_LT'] > 0 else 0

    print(f"Total modes              = {n_total}", flush=True)
    print(f"|δ| ≤ max(Λ^off, Θ^off) = {n_bound_ok}/{n_total}  (Abel-Bound korrekt)", flush=True)
    print(f"max(Λ^off, Θ^off) < 1   = {n_max_lt1}/{n_total}  ← ENTSCHEIDEND", flush=True)
    print(f"max max(Λ,Θ)             = {max_max:.6e}", flush=True)
    print(f"max |δ|                  = {max_delta:.6e}", flush=True)
    print(f"Worst mode: m={worst['m']}, max(Λ,Θ)={worst['max_LT']:.4e}", flush=True)
    print(f"Tightness (|δ|/bound)    = {tight_ratio:.4f}  (mode {tightest['m']})", flush=True)

    return {
        'n_total': n_total, 'n_max_lt1': n_max_lt1,
        'max_max': max_max, 'max_delta': max_delta,
        'total_off_mass': total_off_mass,
        'results': results,
    }


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = offcluster_abel_test(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    print(f"\n{'='*130}", flush=True)
    print("Cross-λ: Off-Cluster Abel Bound Summary", flush=True)
    print(f"{'='*130}", flush=True)
    print(f"{'λ':>4}  {'n':>4}  {'M_off':>10}  {'max(Λ,Θ)':>10}  {'max|δ|':>10}  "
          f"{'max<1':>6}  {'tight':>8}  {'SCHLUSS':>10}", flush=True)
    for lam in sorted(all_results):
        r = all_results[lam]
        tight = max(abs(x['delta'])/x['max_LT'] for x in r['results'] if x['max_LT'] > 0)
        ok = '✓ C2o.1' if r['n_max_lt1'] == r['n_total'] else '✗ OFFEN'
        print(f"{lam:4.0f}  {r['n_total']:4d}  {r['total_off_mass']:10.4e}  {r['max_max']:10.4e}  "
              f"{r['max_delta']:10.4e}  {r['n_max_lt1']}/{r['n_total']:>3}  {tight:8.4f}  {ok:>10}", flush=True)

    print("\nDone.", flush=True)
