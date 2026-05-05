"""
c2ao_matched_sidemass_test.py — Matched Side-Mass Lemma Diagnostic

Zeigt für jeden regulären Mode:
  - Welche Seite hat den kleineren Gap? (g_L oder g_R)
  - Wie groß ist die zugehörige Seitenmasse? (B_{m-1}^off oder C_m^off)
  - Wie groß ist das Produkt (d/gap_danger) × (Seitenmasse_danger)?
  - Vergleich: C2AK (M_off × d/g_min) vs C2AJ (matched side-mass)

Kernhypothese: Der gefährliche kleine Gap paart sich IMMER mit der kleinen Seitenmasse.

Autor: LG (Opus 4.6, Session 27b+++)
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


def matched_sidemass_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*150}", flush=True)
    print(f"C2ao Matched Side-Mass: λ={lam}, N={N}", flush=True)
    print(f"{'='*150}", flush=True)

    t0 = time.time()
    ops = build_operators(lam, N)
    print(f"Build: {time.time()-t0:.1f}s", flush=True)

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
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    u_poles = [w_min] + [w_arr[a] for a in noncl_sorted]
    M = len(noncl_sorted)

    b_off = np.zeros(M)
    for idx, a in enumerate(noncl_sorted):
        b_off[idx] = c_arr[a] * alpha_a[a] / alpha
    M_off = np.sum(np.abs(b_off))

    B_cum = np.zeros(M + 1)
    for n in range(M):
        B_cum[n + 1] = B_cum[n] + b_off[n]

    C_cum = np.zeros(M + 1)
    C_cum[M] = 0
    for n in range(M - 1, -1, -1):
        C_cum[n] = C_cum[n + 1] + b_off[n]

    t_all_reg = sorted([t_vals[j] for j in j_reg])
    s_zeros = []
    t_matched = set()
    for m in range(1, M + 1):
        u_left = u_poles[m - 1]
        u_right = u_poles[m]
        cands = [t for t in t_all_reg if u_left < t < u_right and t not in t_matched]
        if cands:
            s = cands[0]
            s_zeros.append({'m': m, 's': s, 'u_left': u_left, 'u_right': u_right})
            t_matched.add(s)

    print(f"\n  β₀={A_cl:.6f}  M_off={M_off:.4e}  M={M} off-cluster poles", flush=True)

    print(f"\n{'m':>3}  {'gL':>10}  {'gR':>10}  {'danger':>6}  "
          f"{'B_left':>10}  {'C_right':>10}  {'side_mass':>10}  "
          f"{'Λ^off':>10}  {'Θ^off':>10}  {'max(Λ,Θ)':>10}  "
          f"{'C2AK':>10}  {'ratio':>8}  {'<1':>3}", flush=True)
    print('-' * 145, flush=True)

    results = []
    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        d_m = s_m - w_min
        gL = s_m - entry['u_left']
        gR = entry['u_right'] - s_m
        g_min = min(gL, gR)

        B_left = abs(B_cum[m - 1])
        C_right = abs(C_cum[m])

        Lambda_off = d_m * B_left / gL if gL > 1e-50 else float('inf')
        Theta_off = d_m * C_right / gR if gR > 1e-50 else float('inf')
        c2aj = max(Lambda_off, Theta_off)

        c2ak = M_off * d_m / g_min if g_min > 1e-50 else float('inf')

        danger = 'LEFT' if gL <= gR else 'RIGHT'
        side_mass = B_left if gL <= gR else C_right
        ratio = c2aj / c2ak if c2ak > 1e-50 else 0

        ok = '✓' if c2aj < 1 else '✗'

        results.append({
            'm': m, 'gL': gL, 'gR': gR, 'danger': danger,
            'B_left': B_left, 'C_right': C_right, 'side_mass': side_mass,
            'Lambda': Lambda_off, 'Theta': Theta_off,
            'c2aj': c2aj, 'c2ak': c2ak, 'ratio': ratio,
        })

        if m <= 6 or m >= M - 3 or c2aj >= 0.01 or c2ak >= 1:
            print(f"{m:3d}  {gL:10.2e}  {gR:10.2e}  {danger:>6}  "
                  f"{B_left:10.4e}  {C_right:10.4e}  {side_mass:10.4e}  "
                  f"{Lambda_off:10.4e}  {Theta_off:10.4e}  {c2aj:10.4e}  "
                  f"{c2ak:10.4e}  {ratio:8.2e}  {ok:>3}", flush=True)
        elif m == 7:
            print("  ...", flush=True)

    print(f"\n--- Matched Side-Mass Analyse λ={lam} ---", flush=True)
    n_total = len(results)

    n_danger_left = sum(1 for r in results if r['danger'] == 'LEFT')
    n_danger_right = sum(1 for r in results if r['danger'] == 'RIGHT')

    max_c2aj = max(r['c2aj'] for r in results)
    max_c2ak = max(r['c2ak'] for r in results)

    max_side_mass_danger = max(r['side_mass'] for r in results)
    mean_ratio = np.mean([r['ratio'] for r in results])
    max_ratio = max(r['ratio'] for r in results)

    n_side_smaller = sum(1 for r in results
                         if (r['danger'] == 'LEFT' and r['B_left'] <= r['C_right'])
                         or (r['danger'] == 'RIGHT' and r['C_right'] <= r['B_left']))

    print(f"  Total modes              = {n_total}", flush=True)
    print(f"  Danger LEFT / RIGHT      = {n_danger_left} / {n_danger_right}", flush=True)
    print(f"  Side-mass < other side   = {n_side_smaller}/{n_total}  ← MATCHED PAIRING", flush=True)
    print(f"  max C2AJ                 = {max_c2aj:.6e}  (Λ^off oder Θ^off)", flush=True)
    print(f"  max C2AK                 = {max_c2ak:.6e}  (M_off · d/g_min)", flush=True)
    print(f"  max C2AJ/C2AK ratio      = {max_ratio:.4e}  (matched/raw)", flush=True)
    print(f"  mean C2AJ/C2AK ratio     = {mean_ratio:.4e}", flush=True)
    print(f"  max danger-side mass     = {max_side_mass_danger:.6e}", flush=True)

    front_mass_50 = 0
    for n in range(min(3, M)):
        front_mass_50 += abs(b_off[n])
    front_frac = front_mass_50 / M_off if M_off > 1e-50 else 0
    print(f"  First 3 modes mass frac  = {front_frac:.4f}  (front-loading)", flush=True)

    return {
        'n_total': n_total, 'max_c2aj': max_c2aj, 'max_c2ak': max_c2ak,
        'max_ratio': max_ratio, 'mean_ratio': mean_ratio,
        'n_side_smaller': n_side_smaller, 'M_off': M_off,
        'n_danger_left': n_danger_left, 'n_danger_right': n_danger_right,
        'front_frac': front_frac,
    }


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = matched_sidemass_test(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    print(f"\n{'='*150}", flush=True)
    print("Cross-λ: Matched Side-Mass Summary", flush=True)
    print(f"{'='*150}", flush=True)
    print(f"{'λ':>4}  {'n':>4}  {'M_off':>10}  {'max C2AJ':>10}  {'max C2AK':>10}  "
          f"{'max ratio':>10}  {'mean ratio':>10}  {'side<other':>10}  "
          f"{'dL/dR':>8}  {'front3':>8}  {'SCHLUSS':>10}", flush=True)
    for lam in sorted(all_results):
        r = all_results[lam]
        ok = '✓ MATCHED' if r['max_c2aj'] < 1 else '✗'
        print(f"{lam:4.0f}  {r['n_total']:4d}  {r['M_off']:10.4e}  {r['max_c2aj']:10.4e}  {r['max_c2ak']:10.4e}  "
              f"{r['max_ratio']:10.4e}  {r['mean_ratio']:10.4e}  "
              f"{r['n_side_smaller']:>3}/{r['n_total']}  "
              f"{r['n_danger_left']}/{r['n_danger_right']:>3}  "
              f"{r['front_frac']:8.4f}  {ok:>10}", flush=True)

    print("\nDone.", flush=True)
