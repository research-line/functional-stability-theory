"""
c2q_stieltjes_bound.py — C2q' Near-Cluster Stieltjes Bound

Berechnet den GPT-Two-Lemma-Bound:
  |δ_j| ≤ (M_off/α) · (t_j - w_min) / min_a|t_j - w_a|

  Mass Lemma:   M_off/α klein
  Geometry Lemma: Interlacing-Ratio beschränkt

Autor: LG (Opus 4.6, Session 27b+)
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
]


def analyze_stieltjes_bound(lam, N):
    mp.mp.dps = DPS
    print(f"\n{'='*80}", flush=True)
    print(f"C2q' Stieltjes-Bound: λ={lam}, N={N}, dps={DPS}", flush=True)
    print(f"{'='*80}", flush=True)

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
    cl_A = ops['cl_A']
    t_vals = ops['t_vals']
    f_arr = ops['f_arr']
    g_arr = ops['g_arr']
    j_reg = ops['j_reg']
    w_min = ops['w_min']

    # === Mass Lemma: M_off ===
    M_off = sum(abs(c_arr[a] * alpha_a[a]) for a in noncl_A)
    alpha_cl = sum(c_arr[a] * alpha_a[a] for a in cl_A)
    mass_ratio = M_off / abs(alpha)

    print(f"\n--- Mass Lemma ---", flush=True)
    print(f"α         = {alpha:.6e}", flush=True)
    print(f"α_cl      = {alpha_cl:.6e}", flush=True)
    print(f"M_off     = {M_off:.6e}", flush=True)
    print(f"M_off/|α| = {mass_ratio:.6e}", flush=True)
    print(f"|noncl_A|  = {len(noncl_A)}", flush=True)

    # === Geometry Lemma: per-mode interlacing ratio ===
    print(f"\n--- Geometry Lemma (per mode) ---", flush=True)
    hdr = f"{'j':>4}  {'t_j':>12}  {'t-w_min':>10}  {'min|t-wa|':>10}  {'ratio':>10}  {'bound':>10}  {'delta':>10}  {'<1':>4}  {'tight':>6}"
    print(hdr, flush=True)
    print('-' * len(hdr), flush=True)

    results = []
    for j in sorted(j_reg, key=lambda j: t_vals[j]):
        d = t_vals[j] - w_min
        if d < 1e-30:
            continue

        min_gap = min(abs(t_vals[j] - w_arr[a]) for a in noncl_A)
        if min_gap < 1e-30:
            continue

        ratio = d / min_gap
        bound = mass_ratio * ratio

        # actual delta via spectral sum
        spec_sum = sum(
            c_arr[a] * alpha_a[a] / (t_vals[j] - w_arr[a])
            for a in noncl_A
            if abs(t_vals[j] - w_arr[a]) > 1e-30
        )
        delta_act = -(d / alpha) * spec_sum

        ok = 'JA' if bound < 1 else 'NEIN'
        tight = abs(delta_act) / bound if bound > 1e-30 else 0

        results.append({
            'j': j, 't_j': t_vals[j], 'd': d, 'min_gap': min_gap,
            'ratio': ratio, 'bound': bound, 'delta': delta_act, 'tight': tight,
        })

        print(f"{j:4d}  {t_vals[j]:12.6e}  {d:10.3e}  {min_gap:10.3e}  "
              f"{ratio:10.3e}  {bound:10.3e}  {delta_act:10.3e}  {ok:>4}  {tight:6.3f}",
              flush=True)

    # === Summary ===
    if results:
        max_ratio = max(r['ratio'] for r in results)
        max_bound = max(r['bound'] for r in results)
        max_delta = max(abs(r['delta']) for r in results)
        n_ok = sum(1 for r in results if r['bound'] < 1)
        n_total = len(results)

        print(f"\n--- Summary λ={lam} ---", flush=True)
        print(f"Modes analyzed:  {n_total}", flush=True)
        print(f"Bound < 1:       {n_ok}/{n_total}", flush=True)
        print(f"max ratio:       {max_ratio:.4f}", flush=True)
        print(f"max bound:       {max_bound:.6e}", flush=True)
        print(f"max |δ|:         {max_delta:.6e}", flush=True)
        print(f"M_off/|α|:       {mass_ratio:.6e}", flush=True)
        print(f"Product (M·R):   {mass_ratio * max_ratio:.6e}", flush=True)
        print(f"RESULT:          {'BOUND HOLDS' if max_bound < 1 else 'BOUND FAILS'}", flush=True)

    return results


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = analyze_stieltjes_bound(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    # === Cross-λ comparison ===
    print(f"\n{'='*80}", flush=True)
    print("Cross-λ Comparison", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"{'λ':>4}  {'M_off/α':>10}  {'max_ratio':>10}  {'product':>10}  {'max|δ|':>10}  {'bound<1':>8}",
          flush=True)
    for lam, res in sorted(all_results.items()):
        if res:
            mr = max(r['bound'] / max(r['ratio'], 1e-30) * max(r['ratio'], 1e-30) for r in res)
            mass_r = res[0]['bound'] / res[0]['ratio'] if res[0]['ratio'] > 0 else 0
            max_rat = max(r['ratio'] for r in res)
            max_bnd = max(r['bound'] for r in res)
            max_d = max(abs(r['delta']) for r in res)
            n_ok = sum(1 for r in res if r['bound'] < 1)
            print(f"{lam:4.0f}  {mass_r:10.3e}  {max_rat:10.3f}  {mass_r*max_rat:10.3e}  "
                  f"{max_d:10.3e}  {n_ok}/{len(res):>3}", flush=True)

    print("\nDone.", flush=True)
