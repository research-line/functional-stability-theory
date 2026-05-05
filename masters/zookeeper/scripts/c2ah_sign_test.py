"""
c2ah_sign_test.py — Vorzeichentest: δ_m ≥ 0 bzw. F(s_m) ≤ 0?

Prüft die entscheidende Restfrage aus C2AH:
Wenn δ_m ≥ 0 auf allen regulären Moden, dann reicht Θ_m < 1 allein.

Berechnet auch Λ_m und max(Λ_m, Θ_m) zur Verifikation von C2ah.2.

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


def sign_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*130}", flush=True)
    print(f"C2ah Sign Test: λ={lam}, N={N}", flush=True)
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

    # b_n coefficients (signed)
    b_signed = np.zeros(M + 1)
    b_signed[0] = A_cl
    for idx, a in enumerate(noncl_sorted):
        b_signed[idx + 1] = c_arr[a] * alpha_a[a] / alpha

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

    print(f"\n{'m':>4}  {'s_m':>12}  {'d_m':>10}  {'δ_m':>12}  {'sgn':>4}  "
          f"{'Λ_m':>10}  {'Θ_m':>10}  {'max(Λ,Θ)':>10}  {'|δ|≤max':>8}  "
          f"{'|δ|≤Θ':>7}  {'Θ<1':>5}  {'active':>7}", flush=True)
    print('-' * 140, flush=True)

    results = []
    n_pos = 0
    n_neg = 0
    n_theta_controls = 0

    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        u_left = entry['u_left']
        u_right = entry['u_right']
        d_m = s_m - w_min
        gap_L = s_m - u_left
        gap_R = u_right - s_m

        # δ_m with sign
        delta_m = 0
        if abs(d_m) > 1e-30:
            delta_m = -d_m * sum(c_arr[a] * alpha_a[a] / (alpha * (s_m - w_arr[a]))
                                 for a in noncl_A if abs(s_m - w_arr[a]) > 1e-30)

        # Λ_m and Θ_m
        B_prev = sum(b_signed[n] for n in range(m))
        C_curr = sum(b_signed[n] for n in range(m, M + 1))

        Lambda_m = d_m * abs(B_prev) / gap_L if gap_L > 1e-50 else float('inf')
        Theta_m = d_m * abs(C_curr) / gap_R if gap_R > 1e-50 else float('inf')
        max_LT = max(Lambda_m, Theta_m)

        sgn = '+' if delta_m >= 0 else '-'
        if delta_m >= 0:
            n_pos += 1
        else:
            n_neg += 1

        bound_ok = '✓' if abs(delta_m) <= max_LT * 1.001 else '✗'
        theta_ok = '✓' if abs(delta_m) <= Theta_m * 1.001 else '✗'
        theta_lt1 = '✓' if Theta_m < 1 else '✗'

        if delta_m >= 0 and Theta_m < 1:
            active = 'Θ'
            n_theta_controls += 1
        elif delta_m < 0 and Lambda_m < 1:
            active = 'Λ'
            n_theta_controls += 0
        else:
            active = '?'

        results.append({
            'm': m, 's': s_m, 'd': d_m, 'delta': delta_m,
            'Lambda': Lambda_m, 'Theta': Theta_m, 'max_LT': max_LT,
            'sgn': sgn, 'active': active,
        })

        if m <= 8 or m >= M - 3 or Theta_m >= 0.1 or Lambda_m >= 0.5:
            print(f"{m:4d}  {s_m:12.6e}  {d_m:10.4e}  {delta_m:12.4e}  {sgn:>4}  "
                  f"{Lambda_m:10.4e}  {Theta_m:10.4e}  {max_LT:10.4e}  {bound_ok:>8}  "
                  f"{theta_ok:>7}  {theta_lt1:>5}  {active:>7}", flush=True)
        elif m == 9:
            print("  ...", flush=True)

    # Summary
    print(f"\n--- Summary λ={lam} ---", flush=True)
    n_total = len(results)
    n_bound_ok = sum(1 for r in results if abs(r['delta']) <= r['max_LT'] * 1.001)
    n_theta_suff = sum(1 for r in results if abs(r['delta']) <= r['Theta'] * 1.001)
    n_theta_lt1 = sum(1 for r in results if r['Theta'] < 1)
    n_lambda_lt1 = sum(1 for r in results if r['Lambda'] < 1)
    max_theta = max(r['Theta'] for r in results)
    max_lambda = max(r['Lambda'] for r in results)

    print(f"Total modes          = {n_total}", flush=True)
    print(f"δ_m ≥ 0  (positive)  = {n_pos}/{n_total}  ← ENTSCHEIDEND für C2ah.3", flush=True)
    print(f"δ_m < 0  (negative)  = {n_neg}/{n_total}", flush=True)
    print(f"|δ| ≤ max(Λ,Θ)      = {n_bound_ok}/{n_total}  (C2ah.2)", flush=True)
    print(f"|δ| ≤ Θ             = {n_theta_suff}/{n_total}  (reine Tail-Hazard)", flush=True)
    print(f"Θ < 1               = {n_theta_lt1}/{n_total}  (max: {max_theta:.4e})", flush=True)
    print(f"Λ < 1               = {n_lambda_lt1}/{n_total}  (max: {max_lambda:.4e})", flush=True)
    print(f"Θ controls (δ≥0∧Θ<1)= {n_theta_controls}/{n_total}", flush=True)

    return {
        'n_pos': n_pos, 'n_neg': n_neg, 'n_total': n_total,
        'n_theta_suff': n_theta_suff, 'n_theta_lt1': n_theta_lt1,
        'max_theta': max_theta, 'max_lambda': max_lambda,
        'results': results,
    }


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = sign_test(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    print(f"\n{'='*130}", flush=True)
    print("Cross-λ: C2ah Sign-Resolved Hazard Summary", flush=True)
    print(f"{'='*130}", flush=True)
    print(f"{'λ':>4}  {'n':>4}  {'δ≥0':>5}  {'δ<0':>5}  {'|δ|≤Θ':>6}  "
          f"{'Θ<1':>5}  {'max Θ':>10}  {'max Λ':>10}  {'Θ allein reicht':>16}", flush=True)
    for lam in sorted(all_results):
        r = all_results[lam]
        theta_enough = r['n_pos'] == r['n_total'] and r['n_theta_lt1'] == r['n_total']
        print(f"{lam:4.0f}  {r['n_total']:4d}  {r['n_pos']:5d}  {r['n_neg']:5d}  "
              f"{r['n_theta_suff']:6d}  {r['n_theta_lt1']:5d}  {r['max_theta']:10.4e}  "
              f"{r['max_lambda']:10.4e}  {'✓ JA' if theta_enough else '✗ NEIN':>16}", flush=True)

    print("\nDone.", flush=True)
