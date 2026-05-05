"""
c2r_abel_cancellation.py — Abel/Stieltjes Cancellation Lemma (C2r)

Testet GPTs Abel-Summations-Ansatz für die Cauchy-Summe:

  F_j(t) = Σ_{a∉Cl} b_a / (t - w_a),   b_a := c_a α_a / α > 0

Statt |F_j| ≤ Σ |b_a|/|t-w_a| (Triangle Inequality, VERSAGT bei C2q')
nutze Abel-Transformation (Summation by parts):

  Ordne w_a aufsteigend. Für w_m < t_j < w_{m+1}:

  Linke Seite (a ≤ m):  Σ b_a/(t-w_a)  mit t-w_a > 0
  Rechte Seite (a > m):  Σ b_a/(t-w_a)  mit t-w_a < 0

  Abel-Transformation jeder Seite → Differenzkerne:
    (w_{a+1} - w_a) / ((t-w_a)(t-w_{a+1}))
  statt roher Polterme 1/(t-w_a).

  Die Differenzkerne sind POSITIV (bei richtiger Orientierung),
  daher keine Vorzeichenkompensation nötig → schärfere Bounds.

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
]


def abel_analysis(lam, N):
    """Abel-Summation der Cauchy-Summe F_j(t) = Σ b_a/(t-w_a)."""
    mp.mp.dps = DPS
    print(f"\n{'='*90}", flush=True)
    print(f"C2r Abel-Cancellation: λ={lam}, N={N}, dps={DPS}", flush=True)
    print(f"{'='*90}", flush=True)

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

    # b_a := c_a * alpha_a / alpha (normalized weights)
    b_arr = np.array([c_arr[a] * alpha_a[a] / alpha for a in range(len(c_arr))])

    # Sort non-cluster eigenvalues
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    w_sorted = np.array([w_arr[a] for a in noncl_sorted])
    b_sorted = np.array([b_arr[a] for a in noncl_sorted])
    n_off = len(noncl_sorted)

    # Cumulative masses
    B_left = np.cumsum(b_sorted)        # B_k = Σ_{a≤k} b_a
    B_right = np.cumsum(b_sorted[::-1])[::-1]  # C_k = Σ_{a≥k} b_a
    B_total = B_left[-1]  # = Σ b_a = (α - α_cl)/α

    print(f"\n--- Off-Cluster Gewichte ---", flush=True)
    print(f"B_total = Σ b_a = {B_total:.6e}", flush=True)
    print(f"(α-α_cl)/α      = {(alpha - ops['alpha_cl'])/alpha:.6e}", flush=True)
    print(f"Alle b_a > 0?    {np.all(b_sorted > 0)}", flush=True)
    print(f"min(b_a)         = {np.min(b_sorted):.4e}", flush=True)
    print(f"max(b_a)         = {np.max(b_sorted):.4e}", flush=True)

    # === Phase 1: Direct vs Abel per mode ===
    print(f"\n--- Phase 1: Cauchy-Summe direkt vs Abel-transformiert ---", flush=True)
    hdr = (f"{'j':>4}  {'t_j':>12}  {'m':>3}  {'F_direct':>12}  {'F_abel':>12}  "
           f"{'δ_j':>12}  {'|δ|<1':>5}  {'abel_bnd':>12}  {'a_bnd<1':>6}")
    print(hdr, flush=True)
    print('-' * len(hdr), flush=True)

    results = []
    for j in sorted(j_reg, key=lambda j: t_vals[j]):
        t = t_vals[j]
        d = t - w_min
        if d < 1e-30:
            continue

        # Direct Cauchy sum
        F_direct = sum(b_sorted[k] / (t - w_sorted[k])
                       for k in range(n_off)
                       if abs(t - w_sorted[k]) > 1e-30)

        # Find position: w_sorted[m] < t < w_sorted[m+1]
        m = np.searchsorted(w_sorted, t) - 1  # last index where w ≤ t

        # Abel transformation — LEFT side (a = 0..m, t - w_a > 0)
        # Σ_{a=0}^{m} b_a/(t-w_a) = B_0/(t-w_0)
        #   + Σ_{a=0}^{m-1} B_a · [(1/(t-w_{a+1})) - (1/(t-w_a))]
        # where B_a = Σ_{i=0}^{a} b_i (partial sums)
        #
        # Using Abel: Σ b_a f(a) = B_m f(m) - Σ_{a=0}^{m-1} B_a [f(a+1)-f(a)]
        # with f(a) = 1/(t-w_sorted[a])
        #
        # Kernel: f(a+1)-f(a) = 1/(t-w_{a+1}) - 1/(t-w_a)
        #       = (w_{a+1}-w_a) / ((t-w_a)(t-w_{a+1}))
        # For a ≤ m: t > w_a, t > w_{a+1} → kernel > 0 (same sign denom)
        # Wait: for a = m-1: w_m < t but w_{m+1} could be > t
        # Actually we only go to a = m, so w_{a+1} ≤ w_{m+1} which is fine for a < m

        abel_left = 0.0
        abel_left_bound = 0.0
        if m >= 0:
            # Boundary term: B_m / (t - w_sorted[m])
            if abs(t - w_sorted[m]) > 1e-30:
                abel_left = B_left[m] / (t - w_sorted[m])
                abel_left_bound = abs(B_left[m]) / abs(t - w_sorted[m])

            # Abel sum terms
            for k in range(m):
                gap = w_sorted[k+1] - w_sorted[k]
                denom1 = t - w_sorted[k]
                denom2 = t - w_sorted[k+1]
                if abs(denom1) > 1e-30 and abs(denom2) > 1e-30:
                    kernel = gap / (denom1 * denom2)
                    abel_left -= B_left[k] * kernel
                    abel_left_bound += abs(B_left[k]) * abs(kernel)

        # Abel transformation — RIGHT side (a = m+1..n-1, t - w_a < 0)
        # Σ_{a=m+1}^{n-1} b_a/(t-w_a)
        # Abel with f(a) = 1/(t-w_sorted[a]), going from right
        # Σ b_a f(a) = C_{m+1} f(m+1) - Σ_{a=m+1}^{n-2} C_{a+1} [f(a+1)-f(a)]
        # where C_a = Σ_{i=a}^{n-1} b_i (right partial sums)
        abel_right = 0.0
        abel_right_bound = 0.0
        if m + 1 < n_off:
            # Boundary term: C_{m+1} / (t - w_sorted[m+1])
            # Note: t < w_sorted[m+1], so t - w_sorted[m+1] < 0
            if abs(t - w_sorted[m+1]) > 1e-30:
                abel_right = B_right[m+1] / (t - w_sorted[m+1])
                abel_right_bound = abs(B_right[m+1]) / abs(t - w_sorted[m+1])

            # Abel sum terms (from right)
            for k in range(m+1, n_off - 1):
                gap = w_sorted[k+1] - w_sorted[k]
                denom1 = t - w_sorted[k]
                denom2 = t - w_sorted[k+1]
                if abs(denom1) > 1e-30 and abs(denom2) > 1e-30:
                    kernel = gap / (denom1 * denom2)
                    abel_right -= B_right[k+1] * kernel
                    abel_right_bound += abs(B_right[k+1]) * abs(kernel)

        F_abel = abel_left + abel_right
        abel_bound = abel_left_bound + abel_right_bound

        delta = -d * F_direct
        delta_abel_bound = d * abel_bound

        ok_delta = 'JA' if abs(delta) < 1 else 'NEIN'
        ok_abel = 'JA' if delta_abel_bound < 1 else 'NEIN'

        results.append({
            'j': j, 't': t, 'm': m,
            'F_direct': F_direct, 'F_abel': F_abel,
            'delta': delta, 'abel_bound': delta_abel_bound,
            'abel_left_bnd': abel_left_bound * d,
            'abel_right_bnd': abel_right_bound * d,
        })

        # Verify Abel = Direct
        rel_err = abs(F_direct - F_abel) / max(abs(F_direct), 1e-30) if abs(F_direct) > 1e-20 else 0

        print(f"{j:4d}  {t:12.6e}  {m:3d}  {F_direct:12.4e}  {F_abel:12.4e}  "
              f"{delta:12.4e}  {ok_delta:>5}  {delta_abel_bound:12.4e}  {ok_abel:>6}",
              flush=True)

    # === Phase 2: Summary ===
    if results:
        max_delta = max(abs(r['delta']) for r in results)
        max_abel = max(r['abel_bound'] for r in results)
        n_delta_ok = sum(1 for r in results if abs(r['delta']) < 1)
        n_abel_ok = sum(1 for r in results if r['abel_bound'] < 1)
        n_total = len(results)

        print(f"\n--- Summary λ={lam} ---", flush=True)
        print(f"Modes analyzed:    {n_total}", flush=True)
        print(f"|δ_j| < 1:         {n_delta_ok}/{n_total}", flush=True)
        print(f"Abel bound < 1:    {n_abel_ok}/{n_total}", flush=True)
        print(f"max |δ_j|:         {max_delta:.6e}", flush=True)
        print(f"max Abel bound:    {max_abel:.6e}", flush=True)
        print(f"Abel tightness:    {max_delta/max_abel:.4f}" if max_abel > 0 else "", flush=True)

        # Worst modes detail
        worst = sorted(results, key=lambda r: r['abel_bound'], reverse=True)[:5]
        print(f"\n--- Worst 5 modes (Abel bound) ---", flush=True)
        for r in worst:
            tight = abs(r['delta']) / r['abel_bound'] if r['abel_bound'] > 1e-30 else 0
            print(f"  j={r['j']:3d}: Abel bound={r['abel_bound']:.4e}, "
                  f"|δ|={abs(r['delta']):.4e}, tight={tight:.4f}, "
                  f"L={r['abel_left_bnd']:.4e}, R={r['abel_right_bnd']:.4e}",
                  flush=True)

    # === Phase 3: Boundary term analysis ===
    print(f"\n--- Phase 3: Boundary term vs Abel-Kerne ---", flush=True)
    print("Frage: Dominiert der Randterm B_m/(t-w_m) oder die Abel-Differenzkerne?", flush=True)
    for r in results[:10]:
        j = r['j']
        t = r['t']
        m = r['m']
        if m < 0 or m >= n_off:
            continue

        # Left boundary
        left_bdy = abs(B_left[m] / (t - w_sorted[m])) * (t - w_min) if abs(t - w_sorted[m]) > 1e-30 else 0
        # Right boundary
        right_bdy = abs(B_right[m+1] / (t - w_sorted[m+1])) * (t - w_min) if m+1 < n_off and abs(t - w_sorted[m+1]) > 1e-30 else 0
        total_bdy = left_bdy + right_bdy
        total_abel = r['abel_bound']
        bdy_frac = total_bdy / total_abel if total_abel > 1e-30 else 0

        print(f"  j={j:3d}: boundary={total_bdy:.4e}, total={total_abel:.4e}, "
              f"frac={bdy_frac:.3f}, B_m={B_left[m]:.4e}",
              flush=True)

    return results


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = abel_analysis(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    # Cross-λ
    print(f"\n{'='*90}", flush=True)
    print("Cross-λ Comparison (C2r Abel Cancellation)", flush=True)
    print(f"{'='*90}", flush=True)
    print(f"{'λ':>4}  {'max|δ|':>10}  {'max Abel':>10}  {'Abel<1':>8}  {'tightness':>10}", flush=True)
    for lam, res in sorted(all_results.items()):
        if res:
            md = max(abs(r['delta']) for r in res)
            ma = max(r['abel_bound'] for r in res)
            n_ok = sum(1 for r in res if r['abel_bound'] < 1)
            tight = md / ma if ma > 0 else 0
            print(f"{lam:4.0f}  {md:10.3e}  {ma:10.3e}  {n_ok}/{len(res):>3}  {tight:10.4f}", flush=True)

    print("\nDone.", flush=True)
