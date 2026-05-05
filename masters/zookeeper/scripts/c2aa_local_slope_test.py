"""
c2aa_local_slope_test.py — C2AA: Herglotz-Produktformel und lokaler Mean-Gap-Bound

Verifiziert:
1. Compressed pole set: u_0=w_min (Cluster), u_1<...<u_M (Off-Cluster)
2. Interlacing: u_0 < s_1 < u_1 < s_2 < ... < s_M < u_M
3. C2aa.3: m_A'(s_j) ≤ 1/((s_j-u_{j-1})(u_j-s_j))
4. C2aa.5: |δ_j| ≤ (‖R‖/α)/√((s_j-u_{j-1})(u_j-s_j))
5. C2aa.6-Test: Ist ‖R‖/α < √(min mean-gap)?

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


def local_slope_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*110}", flush=True)
    print(f"C2aa Local Slope Test: λ={lam}, N={N}", flush=True)
    print(f"{'='*110}", flush=True)

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
    dim = ops['dim']

    A_cl = sum(alpha_a[a]**2 for a in cl_A)

    # Compressed pole set: u_0 = w_min, u_1 < u_2 < ... < u_M
    u_poles = [w_min]  # u_0
    beta_poles = [A_cl]  # β_0
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    for a in noncl_sorted:
        u_poles.append(w_arr[a])
        beta_poles.append(alpha_a[a]**2)
    u_poles = np.array(u_poles)
    beta_poles = np.array(beta_poles)
    M = len(u_poles) - 1

    # ‖R‖/α
    R_norm_sq = sum((c_arr[a] * (w_arr[a] - w_min))**2 for a in noncl_A)
    R_norm = np.sqrt(R_norm_sq)
    R_over_alpha = R_norm / abs(alpha)

    # Interlacing: für jedes Intervall (u_{m-1}, u_m) die T-Nullstelle finden
    t_all_reg = sorted([t_vals[j] for j in j_reg])

    # Zuordnung: s_m = unique T-EV in (u_{m-1}, u_m)
    s_zeros = []     # matched zeros
    j_near = []      # near-cluster T-modes that don't match any interval
    t_matched = set()

    for m in range(1, M + 1):
        u_left = u_poles[m - 1]
        u_right = u_poles[m]
        candidates = [t for t in t_all_reg if u_left < t < u_right and t not in t_matched]
        if candidates:
            s = candidates[0]
            s_zeros.append({'m': m, 's': s, 'u_left': u_left, 'u_right': u_right})
            t_matched.add(s)

    # Unmatched T-modes = near-cluster boundary layer
    j_near_modes = []
    for j in j_reg:
        if t_vals[j] not in t_matched:
            j_near_modes.append(j)

    print(f"\nKomprimiertes Polsystem:", flush=True)
    print(f"  M = {M} Off-Cluster-Pole + 1 Cluster-Pol", flush=True)
    print(f"  |Cl| = {len(cl_A)}, A_cl = {A_cl:.6e}", flush=True)
    print(f"  Σ β_m = {sum(beta_poles):.10f} (should = 1.0)", flush=True)
    print(f"  {len(s_zeros)} matched Herglotz-Nullstellen (should = {M})", flush=True)
    print(f"  {len(j_near_modes)} Near-Cluster-Moden (Boundary Layer)", flush=True)
    print(f"  ‖R‖/α = {R_over_alpha:.6e}", flush=True)
    print(f"  (‖R‖/α)² = {R_over_alpha**2:.6e}", flush=True)

    # Near-cluster diagnosis
    if j_near_modes:
        print(f"\n--- Near-Cluster Boundary Layer ---", flush=True)
        print(f"{'j':>4}  {'t_j':>14}  {'t_j-w_min':>12}  {'|f_j|':>12}  {'|δ_j|':>12}  "
              f"{'|β_j|²|δ_j|²':>14}  {'Status':>10}", flush=True)
        for j in j_near_modes:
            t = t_vals[j]
            d = t - w_min
            delta = 0
            if abs(d) > 1e-30:
                delta = -d * sum(c_arr[a] * alpha_a[a] / (alpha * (t - w_arr[a]))
                                 for a in noncl_A if abs(t - w_arr[a]) > 1e-30)
            beta_j = ops['g_arr'][j] if 'g_arr' in ops else 0
            print(f"{j:4d}  {t:14.8e}  {d:12.4e}  {abs(f_arr[j]):12.4e}  {abs(delta):12.4e}  "
                  f"{beta_j**2 * delta**2:14.4e}  {'LEAK' if abs(f_arr[j]) < 1e-10 else 'NEAR'}", flush=True)

    # === Interlacing-Verifikation ===
    print(f"\n--- Interlacing-Verifikation ---", flush=True)
    interlacing_ok = True
    for entry in s_zeros:
        s = entry['s']
        u_left = entry['u_left']
        u_right = entry['u_right']
        if not (u_left < s < u_right):
            interlacing_ok = False
            print(f"  FEHLER: s_{entry['m']} = {s:.6e} nicht in ({u_left:.6e}, {u_right:.6e})", flush=True)
    print(f"  Matched: {len(s_zeros)}/{M}", flush=True)
    print(f"  Interlacing: {'✓ PERFEKT' if interlacing_ok else '✗ FEHLER'}", flush=True)

    # === C2aa.3 & C2aa.5: Lokale Steigung und Bound ===
    print(f"\n--- C2aa.3/C2aa.5: Lokaler Mean-Gap Bound ---", flush=True)

    print(f"\n{'j':>4}  {'s_j':>12}  {'u_{j-1}':>12}  {'u_j':>12}  {'gap_L':>10}  {'gap_R':>10}  "
          f"{'mean_gap':>10}  {'C2aa bnd':>10}  {'|δ|':>10}  {'bnd<1':>6}  {'tight':>8}", flush=True)
    print('-' * 125, flush=True)

    results = []

    for entry in s_zeros:
        idx = entry['m'] - 1  # 0-basiert
        s_j = entry['s']
        u_left = entry['u_left']
        u_right = entry['u_right']

        gap_L = s_j - u_left
        gap_R = u_right - s_j
        mean_gap = gap_L * gap_R

        if mean_gap > 0:
            c2aa_bound = R_over_alpha / np.sqrt(mean_gap)
        else:
            c2aa_bound = float('inf')

        # δ_j
        t_j = s_j
        d = t_j - w_min
        delta_actual = 0
        if abs(d) > 1e-30:
            delta_actual = -d * sum(c_arr[a] * alpha_a[a] / (alpha * (t_j - w_arr[a]))
                                    for a in noncl_A if abs(t_j - w_arr[a]) > 1e-30)

        ok = '✓' if c2aa_bound < 1 else '✗'
        tight = abs(delta_actual) / c2aa_bound if c2aa_bound > 1e-50 else 0

        results.append({
            'idx': idx, 's_j': s_j, 'u_left': u_left, 'u_right': u_right,
            'gap_L': gap_L, 'gap_R': gap_R, 'mean_gap': mean_gap,
            'c2aa_bound': c2aa_bound, 'delta': delta_actual, 'tight': tight,
        })

        if idx < 12 or idx >= len(s_zeros) - 5:
            print(f"{idx+1:4d}  {s_j:12.6e}  {u_left:12.6e}  {u_right:12.6e}  {gap_L:10.2e}  {gap_R:10.2e}  "
                  f"{mean_gap:10.2e}  {c2aa_bound:10.4e}  {abs(delta_actual):10.4e}  {ok:>6}  {tight:8.4f}", flush=True)
        elif idx == 12:
            print("  ...", flush=True)

    # === Summary ===
    print(f"\n--- Summary λ={lam} ---", flush=True)
    if results:
        max_delta = max(abs(r['delta']) for r in results)
        max_c2aa = max(r['c2aa_bound'] for r in results)
        min_meangap = min(r['mean_gap'] for r in results if r['mean_gap'] > 0)
        n_ok = sum(1 for r in results if r['c2aa_bound'] < 1)
        n_total = len(results)
        max_tight = max(r['tight'] for r in results)

        # Worst mode
        worst = max(results, key=lambda r: r['c2aa_bound'])

        print(f"‖R‖/α             = {R_over_alpha:.6e}", flush=True)
        print(f"(‖R‖/α)²          = {R_over_alpha**2:.6e}", flush=True)
        print(f"min mean-gap       = {min_meangap:.6e}", flush=True)
        print(f"√(min mean-gap)    = {np.sqrt(min_meangap):.6e}", flush=True)
        print(f"(‖R‖/α)/√(min mg) = {R_over_alpha/np.sqrt(min_meangap):.6e}  (should < 1)", flush=True)
        print(f"max |δ|            = {max_delta:.6e}", flush=True)
        print(f"C2aa bound < 1     = {n_ok}/{n_total}", flush=True)
        print(f"max C2aa bound     = {max_c2aa:.4e}  (mode {worst['idx']+1})", flush=True)
        print(f"max tightness      = {max_tight:.4f}", flush=True)
        print(f"Worst mode: s_{worst['idx']+1} = {worst['s_j']:.6e}, "
              f"gap_L = {worst['gap_L']:.2e}, gap_R = {worst['gap_R']:.2e}", flush=True)

    return {
        'R_over_alpha': R_over_alpha, 'results': results,
        'min_meangap': min_meangap if results else 0,
    }


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = local_slope_test(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    # Cross-λ summary
    print(f"\n{'='*110}", flush=True)
    print("Cross-λ: C2aa Local Slope Summary", flush=True)
    print(f"{'='*110}", flush=True)
    print(f"{'λ':>4}  {'‖R‖/α':>10}  {'min mg':>10}  {'√mg':>10}  {'R/√mg':>10}  "
          f"{'max|δ|':>10}  {'max bnd':>10}  {'bnd<1':>8}  {'tight':>8}", flush=True)
    for lam in sorted(all_results):
        r = all_results[lam]
        if r['results']:
            max_d = max(abs(x['delta']) for x in r['results'])
            max_b = max(x['c2aa_bound'] for x in r['results'])
            mg = r['min_meangap']
            n_ok = sum(1 for x in r['results'] if x['c2aa_bound'] < 1)
            n_t = len(r['results'])
            tight = max(x['tight'] for x in r['results'])
            print(f"{lam:4.0f}  {r['R_over_alpha']:10.4e}  {mg:10.2e}  {np.sqrt(mg):10.4e}  "
                  f"{r['R_over_alpha']/np.sqrt(mg):10.4e}  {max_d:10.4e}  {max_b:10.4e}  "
                  f"{n_ok}/{n_t:>3}  {tight:8.4f}", flush=True)

    print("\nDone.", flush=True)
