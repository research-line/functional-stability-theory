"""
c2r_worst_mode_tracking.py — Worst-Mode Tracking + Principal-Pair Analyse

Verfolgt den schlimmsten Abel-Bound über λ=3,5,7 und zerlegt ihn in:
  (1) Principal-pair: die zwei flankierenden Pole w_m, w_{m+1}
  (2) Abel-rest: alle anderen Terme

GPT-Strategie (Session 27b++):
  - Principal-pair control: Der große Beitrag kommt vom kritischen Paar
  - Abel-rest small: Rest muss uniform subdominant bleiben
  - Worst-mode tracking: Immer derselbe geometrische Typ?

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


def principal_pair_analysis(lam, N):
    mp.mp.dps = DPS
    print(f"\n{'='*95}", flush=True)
    print(f"C2r Worst-Mode Tracking: λ={lam}, N={N}, dps={DPS}", flush=True)
    print(f"{'='*95}", flush=True)

    t0 = time.time()
    print("Building operators...", flush=True)
    ops = build_operators(lam, N)
    dt = time.time() - t0
    print(f"Build complete in {dt:.1f}s", flush=True)

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

    print(f"B_total = {B_total:.6e}, n_off = {n_off}", flush=True)
    print(f"b_a range: [{np.min(b_sorted):.4e}, {np.max(b_sorted):.4e}]", flush=True)
    neg_count = np.sum(b_sorted < 0)
    if neg_count > 0:
        print(f"WARNUNG: {neg_count} negative b_a (max |neg| = {np.max(np.abs(b_sorted[b_sorted < 0])):.4e})", flush=True)

    # Per-mode analysis with principal-pair decomposition
    hdr = (f"{'j':>4}  {'t_j':>10}  {'m':>3}  {'gap_L':>9}  {'gap_R':>9}  "
           f"{'PP_cont':>10}  {'Abel_rest':>10}  {'Abel_tot':>10}  "
           f"{'delta':>10}  {'<1':>3}  {'PP%':>5}")
    print(f"\n{hdr}", flush=True)
    print('-' * len(hdr), flush=True)

    results = []
    for j in sorted(j_reg, key=lambda j: t_vals[j]):
        t = t_vals[j]
        d = t - w_min
        if d < 1e-30:
            continue

        m = np.searchsorted(w_sorted, t) - 1

        # Direct Cauchy sum for actual delta
        F_direct = sum(b_sorted[k] / (t - w_sorted[k])
                       for k in range(n_off)
                       if abs(t - w_sorted[k]) > 1e-30)
        delta = -d * F_direct

        # === Principal Pair: the two flanking poles ===
        pp_left = b_sorted[m] / (t - w_sorted[m]) if m >= 0 and abs(t - w_sorted[m]) > 1e-30 else 0
        pp_right = b_sorted[m+1] / (t - w_sorted[m+1]) if m+1 < n_off and abs(t - w_sorted[m+1]) > 1e-30 else 0
        pp_sum = pp_left + pp_right
        pp_delta = abs(-d * pp_sum)

        # === Rest: everything except the principal pair ===
        rest_sum = F_direct - pp_sum if abs(pp_sum) > 1e-30 else F_direct
        rest_delta = abs(-d * rest_sum)

        # === Full Abel bound (same as c2r_abel_cancellation.py) ===
        abel_bound = compute_abel_bound(t, d, m, n_off, w_sorted, b_sorted, B_left, B_right)

        # === Principal-pair Abel bound (only boundary terms) ===
        pp_abel = 0.0
        if m >= 0 and abs(t - w_sorted[m]) > 1e-30:
            pp_abel += abs(B_left[m]) / abs(t - w_sorted[m]) * d
        if m + 1 < n_off and abs(t - w_sorted[m+1]) > 1e-30:
            pp_abel += abs(B_right[m+1]) / abs(t - w_sorted[m+1]) * d

        abel_rest = abel_bound - pp_abel

        gap_L = t - w_sorted[m] if m >= 0 else float('inf')
        gap_R = w_sorted[m+1] - t if m+1 < n_off else float('inf')
        ok = 'JA' if abel_bound < 1 else '!!'
        pp_pct = pp_abel / abel_bound * 100 if abel_bound > 1e-30 else 0

        results.append({
            'j': j, 't': t, 'm': m, 'd': d,
            'gap_L': gap_L, 'gap_R': gap_R,
            'pp_delta': pp_delta, 'rest_delta': rest_delta,
            'abel_bound': abel_bound, 'pp_abel': pp_abel, 'abel_rest': abel_rest,
            'delta': delta, 'pp_pct': pp_pct,
        })

        print(f"{j:4d}  {t:10.4e}  {m:3d}  {gap_L:9.2e}  {gap_R:9.2e}  "
              f"{pp_abel:10.3e}  {abel_rest:10.3e}  {abel_bound:10.3e}  "
              f"{delta:10.3e}  {ok:>3}  {pp_pct:5.1f}",
              flush=True)

    # === Summary ===
    if results:
        max_abel = max(r['abel_bound'] for r in results)
        max_delta = max(abs(r['delta']) for r in results)
        n_ok = sum(1 for r in results if r['abel_bound'] < 1)
        n_total = len(results)
        worst = max(results, key=lambda r: r['abel_bound'])

        print(f"\n--- Summary λ={lam} ---", flush=True)
        print(f"Modes:           {n_total}", flush=True)
        print(f"Abel bound < 1:  {n_ok}/{n_total}", flush=True)
        print(f"max Abel bound:  {max_abel:.6e}", flush=True)
        print(f"max |δ|:         {max_delta:.6e}", flush=True)
        print(f"Tightness:       {max_delta/max_abel:.6f}", flush=True)

        print(f"\n--- Worst mode detail ---", flush=True)
        w = worst
        print(f"  j={w['j']}, m={w['m']}, t={w['t']:.6e}", flush=True)
        print(f"  gap_L={w['gap_L']:.4e}, gap_R={w['gap_R']:.4e}", flush=True)
        print(f"  Abel total={w['abel_bound']:.6e}", flush=True)
        print(f"  PP (boundary)={w['pp_abel']:.6e} ({w['pp_pct']:.1f}%)", flush=True)
        print(f"  Abel rest={w['abel_rest']:.6e}", flush=True)
        print(f"  Actual |δ|={abs(w['delta']):.6e}", flush=True)
        print(f"  B_left[m]={B_left[w['m']]:.6e}" if w['m'] >= 0 else "", flush=True)
        if w['m']+1 < n_off:
            print(f"  B_right[m+1]={B_right[w['m']+1]:.6e}", flush=True)

        # Top 5 worst
        top5 = sorted(results, key=lambda r: r['abel_bound'], reverse=True)[:5]
        print(f"\n--- Top 5 worst modes ---", flush=True)
        for r in top5:
            print(f"  j={r['j']:3d}: abel={r['abel_bound']:.4e}, "
                  f"PP={r['pp_abel']:.4e}({r['pp_pct']:.0f}%), "
                  f"rest={r['abel_rest']:.4e}, "
                  f"|δ|={abs(r['delta']):.4e}, "
                  f"gaps=({r['gap_L']:.2e},{r['gap_R']:.2e})",
                  flush=True)

    return results


def compute_abel_bound(t, d, m, n_off, w_sorted, b_sorted, B_left, B_right):
    """Full Abel bound computation."""
    abel_bound = 0.0

    # Left side
    if m >= 0:
        if abs(t - w_sorted[m]) > 1e-30:
            abel_bound += abs(B_left[m]) / abs(t - w_sorted[m]) * d
        for k in range(m):
            gap = w_sorted[k+1] - w_sorted[k]
            d1 = t - w_sorted[k]
            d2 = t - w_sorted[k+1]
            if abs(d1) > 1e-30 and abs(d2) > 1e-30:
                abel_bound += abs(B_left[k]) * abs(gap / (d1 * d2)) * d

    # Right side
    if m + 1 < n_off:
        if abs(t - w_sorted[m+1]) > 1e-30:
            abel_bound += abs(B_right[m+1]) / abs(t - w_sorted[m+1]) * d
        for k in range(m+1, n_off - 1):
            gap = w_sorted[k+1] - w_sorted[k]
            d1 = t - w_sorted[k]
            d2 = t - w_sorted[k+1]
            if abs(d1) > 1e-30 and abs(d2) > 1e-30:
                abel_bound += abs(B_right[k+1]) * abs(gap / (d1 * d2)) * d

    return abel_bound


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = principal_pair_analysis(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    # Cross-λ tracking
    print(f"\n{'='*95}", flush=True)
    print("Cross-λ Worst-Mode Tracking", flush=True)
    print(f"{'='*95}", flush=True)
    print(f"{'λ':>4}  {'n_modes':>7}  {'Abel<1':>7}  {'max_abel':>10}  {'max|δ|':>10}  "
          f"{'tight':>8}  {'worst_j':>7}  {'PP%':>5}  {'gap_min':>9}", flush=True)
    for lam, res in sorted(all_results.items()):
        if res:
            ma = max(r['abel_bound'] for r in res)
            md = max(abs(r['delta']) for r in res)
            n_ok = sum(1 for r in res if r['abel_bound'] < 1)
            worst = max(res, key=lambda r: r['abel_bound'])
            gap_min = min(min(r['gap_L'], r['gap_R']) for r in res
                         if r['gap_L'] > 1e-30 and r['gap_R'] > 1e-30)
            print(f"{lam:4.0f}  {len(res):7d}  {n_ok}/{len(res):>3}  {ma:10.3e}  {md:10.3e}  "
                  f"{md/ma:8.5f}  j={worst['j']:>3d}  {worst['pp_pct']:5.1f}  {gap_min:9.2e}",
                  flush=True)

    # λ-scaling analysis
    lams = sorted(all_results.keys())
    if len(lams) >= 2:
        print(f"\n--- λ-Scaling ---", flush=True)
        for i in range(1, len(lams)):
            l1, l2 = lams[i-1], lams[i]
            r1 = [r for r in all_results[l1]]
            r2 = [r for r in all_results[l2]]
            ma1 = max(r['abel_bound'] for r in r1)
            ma2 = max(r['abel_bound'] for r in r2)
            md1 = max(abs(r['delta']) for r in r1)
            md2 = max(abs(r['delta']) for r in r2)
            ratio_abel = ma2 / ma1 if ma1 > 0 else 0
            ratio_delta = md2 / md1 if md1 > 0 else 0
            exp_abel = np.log(ratio_abel) / np.log(l2/l1) if ratio_abel > 0 else 0
            exp_delta = np.log(ratio_delta) / np.log(l2/l1) if ratio_delta > 0 else 0
            print(f"  λ={l1:.0f}→{l2:.0f}: Abel ratio={ratio_abel:.2f} (exp≈{exp_abel:.1f}), "
                  f"|δ| ratio={ratio_delta:.3f} (exp≈{exp_delta:.1f})", flush=True)

    print("\nDone.", flush=True)
