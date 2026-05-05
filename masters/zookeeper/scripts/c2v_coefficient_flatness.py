"""
c2v_coefficient_flatness.py — C2v: Rang-1-Null-Identität und Koeffizienten-Flatness

Verifiziert:
1. Exakte Null-Identität: Σ α_a²/(t_j - w_a) = 0
2. A_cl = Σ_{a∈Cl} α_a² ≤ 1
3. Formfaktor ϑ_a = c_a/(α·α_a) — Flatness-Profil
4. Zerlegung δ_j = ϑ* A_cl + Korrektur

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


def coefficient_flatness(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*100}", flush=True)
    print(f"C2v Coefficient Flatness: λ={lam}, N={N}", flush=True)
    print(f"{'='*100}", flush=True)

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

    dim_A = len(w_arr)
    cl_A = sorted(set(range(dim_A)) - set(noncl_A))

    # === Phase 1: Zero Identity Σ α_a²/(t_j - w_a) = 0 ===
    print(f"\n--- Phase 1: Zero Identity Σ α_a²/(t_j - w_a) = 0 ---", flush=True)

    A_cl = sum(alpha_a[a]**2 for a in cl_A)
    A_off = sum(alpha_a[a]**2 for a in noncl_A)
    A_total = A_cl + A_off
    print(f"A_cl   = Σ_{'{a∈Cl}'} α_a² = {A_cl:.10e}", flush=True)
    print(f"A_off  = Σ_{'{a∉Cl}'} α_a² = {A_off:.10e}", flush=True)
    print(f"A_total = {A_total:.10e} (should be 1.0)", flush=True)
    print(f"|Cl| = {len(cl_A)}, |Off| = {len(noncl_A)}, dim = {dim_A}", flush=True)

    zero_id_errors = []
    print(f"\n{'j':>4}  {'t_j':>12}  {'Σ_all':>12}  {'Σ_cl':>12}  {'Σ_off':>12}  "
          f"{'-(t-w)Σ_off':>12}  {'=A_cl?':>8}", flush=True)
    print('-' * 90, flush=True)

    for j in sorted(j_reg, key=lambda j: t_vals[j]):
        t = t_vals[j]

        S_all = sum(alpha_a[a]**2 / (t - w_arr[a])
                     for a in range(dim_A) if abs(t - w_arr[a]) > 1e-30)
        S_cl = sum(alpha_a[a]**2 / (t - w_min) for a in cl_A)
        S_off = sum(alpha_a[a]**2 / (t - w_arr[a])
                     for a in noncl_A if abs(t - w_arr[a]) > 1e-30)
        lhs = -(t - w_min) * S_off

        zero_id_errors.append(abs(S_all))

        if len(zero_id_errors) <= 10 or len(zero_id_errors) > len(j_reg) - 5:
            print(f"{j:4d}  {t:12.6e}  {S_all:12.4e}  {S_cl:12.4e}  {S_off:12.4e}  "
                  f"{lhs:12.6e}  {'✓' if abs(lhs - A_cl) < 1e-6 else '✗':>8}", flush=True)
        elif len(zero_id_errors) == 11:
            print("  ...", flush=True)

    print(f"\nZero Identity: max |Σ α_a²/(t_j-w_a)| = {max(zero_id_errors):.4e}", flush=True)
    print(f"               median                   = {np.median(zero_id_errors):.4e}", flush=True)

    # === Phase 2: Coefficient Profile ϑ_a = c_a/(α·α_a) ===
    print(f"\n--- Phase 2: Formfaktor ϑ_a = c_a/(α·α_a) ---", flush=True)

    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    theta_a = {}
    for a in noncl_sorted:
        if abs(alpha_a[a]) > 1e-30:
            theta_a[a] = c_arr[a] / (alpha * alpha_a[a])

    theta_vals = np.array([theta_a[a] for a in noncl_sorted if a in theta_a])
    alpha2_vals = np.array([alpha_a[a]**2 for a in noncl_sorted if a in theta_a])

    print(f"\n{'a':>4}  {'w_a':>12}  {'α_a':>12}  {'c_a':>12}  {'ϑ_a':>12}  "
          f"{'α_a²':>12}  {'b_a':>12}", flush=True)
    print('-' * 90, flush=True)
    for k, a in enumerate(noncl_sorted):
        if a not in theta_a:
            continue
        th = theta_a[a]
        ba = c_arr[a] * alpha_a[a] / alpha
        if k < 10 or k >= len(noncl_sorted) - 5:
            print(f"{k:4d}  {w_arr[a]:12.6e}  {alpha_a[a]:12.4e}  {c_arr[a]:12.4e}  "
                  f"{th:12.6e}  {alpha_a[a]**2:12.4e}  {ba:12.4e}", flush=True)
        elif k == 10:
            print("  ...", flush=True)

    print(f"\nϑ-Statistik (Off-Cluster):", flush=True)
    print(f"  n        = {len(theta_vals)}", flush=True)
    print(f"  mean     = {np.mean(theta_vals):.6e}", flush=True)
    print(f"  median   = {np.median(theta_vals):.6e}", flush=True)
    print(f"  std      = {np.std(theta_vals):.6e}", flush=True)
    print(f"  CV       = {np.std(theta_vals)/abs(np.mean(theta_vals)):.4f}", flush=True)
    print(f"  min      = {np.min(theta_vals):.6e}", flush=True)
    print(f"  max      = {np.max(theta_vals):.6e}", flush=True)
    print(f"  max/min  = {np.max(theta_vals)/np.min(theta_vals):.4f}" if np.min(theta_vals) > 0
          else f"  max/min  = SIGN CHANGE!", flush=True)

    # Weighted flatness (by α_a²)
    theta_weighted_mean = np.sum(theta_vals * alpha2_vals) / np.sum(alpha2_vals)
    theta_weighted_std = np.sqrt(np.sum(alpha2_vals * (theta_vals - theta_weighted_mean)**2) / np.sum(alpha2_vals))
    print(f"\nα²-gewichtete Statistik:", flush=True)
    print(f"  weighted mean = {theta_weighted_mean:.6e}", flush=True)
    print(f"  weighted std  = {theta_weighted_std:.6e}", flush=True)
    print(f"  weighted CV   = {theta_weighted_std/abs(theta_weighted_mean):.6f}", flush=True)

    # === Phase 3: Zerlegung δ_j = ϑ* A_cl + Korrektur ===
    print(f"\n--- Phase 3: δ_j = ϑ* A_cl + Korrektur ---", flush=True)

    theta_star = np.median(theta_vals)
    print(f"ϑ* (median) = {theta_star:.6e}", flush=True)
    print(f"ϑ* · A_cl   = {theta_star * A_cl:.6e}", flush=True)

    print(f"\n{'j':>4}  {'m':>3}  {'δ_actual':>12}  {'ϑ*·A_cl':>12}  {'correction':>12}  "
          f"{'|δ|<1':>6}  {'|corr|<1-ϑ*A':>14}", flush=True)
    print('-' * 85, flush=True)

    w_sorted = np.array([w_arr[a] for a in noncl_sorted])
    results = []

    for j in sorted(j_reg, key=lambda j: t_vals[j]):
        t = t_vals[j]
        d = t - w_min
        if d < 1e-30:
            continue

        # Actual δ
        F_actual = sum(c_arr[a] * alpha_a[a] / (alpha * (t - w_arr[a]))
                       for a in noncl_A if abs(t - w_arr[a]) > 1e-30)
        delta_actual = -d * F_actual

        # Canonical term: ϑ* A_cl
        canonical = theta_star * A_cl

        # Correction: -(t-w_min) Σ (ϑ_a - ϑ*) α_a² / (t - w_a)
        correction = 0
        for a in noncl_A:
            if a in theta_a and abs(t - w_arr[a]) > 1e-30:
                correction += -(t - w_min) * (theta_a[a] - theta_star) * alpha_a[a]**2 / (t - w_arr[a])

        ok_delta = '✓' if abs(delta_actual) < 1 else '✗'
        margin = 1 - theta_star * A_cl
        ok_corr = '✓' if abs(correction) < margin else '✗'

        m = np.searchsorted(w_sorted, t) - 1

        results.append({
            'j': j, 'm': m, 'delta': delta_actual,
            'canonical': canonical, 'correction': correction,
        })

        if len(results) <= 10 or len(results) > len(j_reg) - 5:
            print(f"{j:4d}  {m:3d}  {delta_actual:12.4e}  {canonical:12.6e}  "
                  f"{correction:12.4e}  {ok_delta:>6}  {ok_corr:>14}", flush=True)
        elif len(results) == 11:
            print("  ...", flush=True)

    # === Phase 4: Summary ===
    if results:
        print(f"\n--- Summary λ={lam} ---", flush=True)
        max_delta = max(abs(r['delta']) for r in results)
        max_corr = max(abs(r['correction']) for r in results)
        n_ok = sum(1 for r in results if abs(r['delta']) < 1)
        n_corr_ok = sum(1 for r in results if abs(r['correction']) < 1 - theta_star * A_cl)

        print(f"A_cl           = {A_cl:.10e}", flush=True)
        print(f"ϑ*             = {theta_star:.6e}", flush=True)
        print(f"ϑ* · A_cl      = {theta_star * A_cl:.6e}", flush=True)
        print(f"Margin (1-ϑ*A) = {1 - theta_star * A_cl:.6e}", flush=True)
        print(f"|δ| < 1        = {n_ok}/{len(results)}", flush=True)
        print(f"max |δ|        = {max_delta:.6e}", flush=True)
        print(f"|corr| < margin= {n_corr_ok}/{len(results)}", flush=True)
        print(f"max |corr|     = {max_corr:.6e}", flush=True)
        print(f"ϑ CV (unweighted) = {np.std(theta_vals)/abs(np.mean(theta_vals)):.4f}", flush=True)
        print(f"ϑ CV (α²-weighted)= {theta_weighted_std/abs(theta_weighted_mean):.6f}", flush=True)

    return {
        'A_cl': A_cl, 'theta_star': theta_star,
        'theta_vals': theta_vals, 'alpha2_vals': alpha2_vals,
        'theta_weighted_mean': theta_weighted_mean,
        'theta_weighted_std': theta_weighted_std,
        'results': results,
    }


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = coefficient_flatness(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    # Cross-λ summary
    print(f"\n{'='*100}", flush=True)
    print("Cross-λ: C2v Coefficient Flatness", flush=True)
    print(f"{'='*100}", flush=True)
    print(f"{'λ':>4}  {'A_cl':>12}  {'ϑ*':>12}  {'ϑ*·A_cl':>12}  {'margin':>12}  "
          f"{'CV_raw':>8}  {'CV_wt':>10}  {'max|δ|':>10}  {'max|corr|':>10}", flush=True)
    for lam in sorted(all_results):
        r = all_results[lam]
        if r['results']:
            max_d = max(abs(x['delta']) for x in r['results'])
            max_c = max(abs(x['correction']) for x in r['results'])
            cv_raw = np.std(r['theta_vals']) / abs(np.mean(r['theta_vals']))
            cv_wt = r['theta_weighted_std'] / abs(r['theta_weighted_mean'])
            margin = 1 - r['theta_star'] * r['A_cl']
            print(f"{lam:4.0f}  {r['A_cl']:12.6e}  {r['theta_star']:12.6e}  "
                  f"{r['theta_star']*r['A_cl']:12.6e}  {margin:12.6e}  "
                  f"{cv_raw:8.4f}  {cv_wt:10.6f}  {max_d:10.4e}  {max_c:10.4e}", flush=True)

    print("\nDone.", flush=True)
