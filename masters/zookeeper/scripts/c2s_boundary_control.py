"""
c2s_boundary_control.py — C2s Abel-Balance Lemma: Boundary-Term Kontrolle

Prüft die reduzierte Beweisaussage:
  |δ_j| < 1  ⟸  d · |P_m(t_j)| < 1

  wobei P_m = B_m/(t-w_m) + C_{m+1}/(w_{m+1}-t)
  und   B_m = Σ_{a≤m} b_a,  C_{m+1} = Σ_{a≥m+1} b_a

Schlüsselfrage: Ist die Massenverteilung b_a clusterkonzentriert?
Wenn ja, dann ist C_{m+1} klein für große m (Bulk-Moden),
und B_m klein für kleine m (Near-Cluster-Moden).

Autor: LG (Opus 4.6, Session 27b++)
Datum: 2026-04-20
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import mpmath as mp
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


def boundary_control(lam, N):
    mp.mp.dps = DPS
    print(f"\n{'='*95}", flush=True)
    print(f"C2s Boundary Control: λ={lam}, N={N}", flush=True)
    print(f"{'='*95}", flush=True)

    t0 = time.time()
    print("Building operators...", flush=True)
    ops = build_operators(lam, N)
    print(f"Build complete in {time.time()-t0:.1f}s", flush=True)

    alpha = ops['alpha']
    w_arr = ops['w_arr']
    c_arr = ops['c_arr']
    alpha_a = ops['alpha_a']
    noncl_A = ops['noncl_A']
    t_vals = ops['t_vals']
    f_arr = ops['f_arr']
    j_reg = ops['j_reg']
    w_min = ops['w_min']

    b_arr = np.array([c_arr[a] * alpha_a[a] / alpha for a in range(len(c_arr))])
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    w_sorted = np.array([w_arr[a] for a in noncl_sorted])
    b_sorted = np.array([b_arr[a] for a in noncl_sorted])
    n_off = len(noncl_sorted)

    B_left = np.cumsum(b_sorted)
    B_right = np.cumsum(b_sorted[::-1])[::-1]
    B_total = B_left[-1]

    # === Phase 1: Massenverteilung b_a ===
    print(f"\n--- Phase 1: Massenverteilung b_a ---", flush=True)
    print(f"B_total = {B_total:.6e}", flush=True)
    print(f"n_off   = {n_off}", flush=True)

    # Cumulative mass profile
    print(f"\n{'a':>4}  {'w_a':>12}  {'b_a':>12}  {'B_left':>12}  {'B_right':>12}  "
          f"{'B_L/B_tot':>9}  {'B_R/B_tot':>9}", flush=True)
    print('-' * 85, flush=True)
    for k in range(min(n_off, 10)):
        print(f"{k:4d}  {w_sorted[k]:12.6e}  {b_sorted[k]:12.4e}  "
              f"{B_left[k]:12.4e}  {B_right[k]:12.4e}  "
              f"{B_left[k]/B_total:9.4f}  {B_right[k]/B_total:9.4f}", flush=True)
    if n_off > 20:
        print("  ...", flush=True)
    for k in range(max(n_off-5, 10), n_off):
        print(f"{k:4d}  {w_sorted[k]:12.6e}  {b_sorted[k]:12.4e}  "
              f"{B_left[k]:12.4e}  {B_right[k]:12.4e}  "
              f"{B_left[k]/B_total:9.4f}  {B_right[k]/B_total:9.4f}", flush=True)

    # Mass concentration
    half_mass_idx = np.searchsorted(B_left, B_total / 2)
    q90_idx = np.searchsorted(B_left, 0.9 * B_total)
    q99_idx = np.searchsorted(B_left, 0.99 * B_total)
    print(f"\nMassenkonzentration:", flush=True)
    print(f"  50% der Masse in ersten {half_mass_idx+1}/{n_off} Eigenwerten", flush=True)
    print(f"  90% der Masse in ersten {q90_idx+1}/{n_off} Eigenwerten", flush=True)
    print(f"  99% der Masse in ersten {q99_idx+1}/{n_off} Eigenwerten", flush=True)

    # === Phase 2: P_m Kontrolle per Mode ===
    print(f"\n--- Phase 2: Boundary-Term |P_m| · d per Mode ---", flush=True)
    hdr = (f"{'j':>4}  {'m':>3}  {'gap_L':>9}  {'gap_R':>9}  "
           f"{'B_m·d/gL':>10}  {'C·d/gR':>10}  {'PP_bound':>10}  "
           f"{'delta':>10}  {'<1':>3}  {'which':>6}")
    print(hdr, flush=True)
    print('-' * len(hdr), flush=True)

    results = []
    for j in sorted(j_reg, key=lambda j: t_vals[j]):
        t = t_vals[j]
        d = t - w_min
        if d < 1e-30:
            continue

        m = np.searchsorted(w_sorted, t) - 1

        # Gaps
        gap_L = t - w_sorted[m] if m >= 0 else float('inf')
        gap_R = w_sorted[m+1] - t if m+1 < n_off else float('inf')

        # Boundary terms
        left_bdy = B_left[m] * d / gap_L if m >= 0 and gap_L > 1e-30 else 0
        right_bdy = B_right[m+1] * d / gap_R if m+1 < n_off and gap_R > 1e-30 else 0
        pp_bound = left_bdy + right_bdy

        # Actual delta
        F_direct = sum(b_sorted[k] / (t - w_sorted[k])
                       for k in range(n_off) if abs(t - w_sorted[k]) > 1e-30)
        delta = -d * F_direct

        ok = 'JA' if pp_bound < 1 else '!!'
        which = 'L' if left_bdy > right_bdy else 'R'

        results.append({
            'j': j, 'm': m, 'd': d,
            'gap_L': gap_L, 'gap_R': gap_R,
            'left_bdy': left_bdy, 'right_bdy': right_bdy,
            'pp_bound': pp_bound, 'delta': delta,
            'B_m': B_left[m] if m >= 0 else 0,
            'C_m1': B_right[m+1] if m+1 < n_off else 0,
        })

        print(f"{j:4d}  {m:3d}  {gap_L:9.2e}  {gap_R:9.2e}  "
              f"{left_bdy:10.3e}  {right_bdy:10.3e}  {pp_bound:10.3e}  "
              f"{delta:10.3e}  {ok:>3}  {which:>6}", flush=True)

    # === Phase 3: Summary ===
    if results:
        max_pp = max(r['pp_bound'] for r in results)
        max_left = max(r['left_bdy'] for r in results)
        max_right = max(r['right_bdy'] for r in results)
        max_delta = max(abs(r['delta']) for r in results)
        n_ok = sum(1 for r in results if r['pp_bound'] < 1)

        print(f"\n--- Summary λ={lam} ---", flush=True)
        print(f"PP bound < 1:    {n_ok}/{len(results)}", flush=True)
        print(f"max PP bound:    {max_pp:.6e}", flush=True)
        print(f"max left B·d/g:  {max_left:.6e}", flush=True)
        print(f"max right C·d/g: {max_right:.6e}", flush=True)
        print(f"max |δ|:         {max_delta:.6e}", flush=True)
        print(f"RESULT:          {'PP BOUND HOLDS' if n_ok == len(results) else 'PP BOUND FAILS'}", flush=True)

        worst = max(results, key=lambda r: r['pp_bound'])
        print(f"\nWorst mode: j={worst['j']}, m={worst['m']}", flush=True)
        print(f"  B_m={worst['B_m']:.4e}, C_{worst['m']+1}={worst['C_m1']:.4e}", flush=True)
        print(f"  gap_L={worst['gap_L']:.4e}, gap_R={worst['gap_R']:.4e}", flush=True)
        print(f"  left={worst['left_bdy']:.4e}, right={worst['right_bdy']:.4e}", flush=True)
        print(f"  PP={worst['pp_bound']:.6e}, |δ|={abs(worst['delta']):.4e}", flush=True)

    return results


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = boundary_control(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    # Cross-λ
    print(f"\n{'='*95}", flush=True)
    print("Cross-λ: C2s Boundary Control", flush=True)
    print(f"{'='*95}", flush=True)
    print(f"{'λ':>4}  {'PP<1':>7}  {'max_PP':>10}  {'max_L':>10}  {'max_R':>10}  "
          f"{'max|δ|':>10}  {'mass_50%':>8}  {'mass_90%':>8}", flush=True)
    for lam, res in sorted(all_results.items()):
        if res:
            mp_val = max(r['pp_bound'] for r in res)
            ml = max(r['left_bdy'] for r in res)
            mr = max(r['right_bdy'] for r in res)
            md = max(abs(r['delta']) for r in res)
            n_ok = sum(1 for r in res if r['pp_bound'] < 1)
            print(f"{lam:4.0f}  {n_ok}/{len(res):>3}  {mp_val:10.3e}  {ml:10.3e}  {mr:10.3e}  "
                  f"{md:10.3e}", flush=True)

    print("\nDone.", flush=True)
