"""
c2as3_canonical_channel_test.py — Test ob C2as.3 gilt

Berechnet die KANONISCHEN Kanäle Λ_m^can, Θ_m^can und vergleicht
mit den REALEN Kanälen Λ_m^off, Θ_m^off.

Kernfrage: Ist max(Λ^can, Θ^can) < 1?
Falls nicht, ist C2as.3 FALSCH und die Beweisarchitektur braucht Revision.
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
]


def canonical_channel_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*120}", flush=True)
    print(f"C2as.3 Canonical Channel Test: λ={lam}, N={N}", flush=True)
    print(f"{'='*120}", flush=True)

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

    beta = np.zeros(M)
    Delta = np.zeros(M)
    a_can = np.zeros(M)
    b_off = np.zeros(M)

    for idx, a in enumerate(noncl_sorted):
        beta[idx] = alpha_a[a]**2
        Delta[idx] = w_arr[a] - w_min
        a_can[idx] = beta[idx] / Delta[idx]
        b_off[idx] = c_arr[a] * alpha_a[a] / alpha

    print(f"\n  β₀ = A_cl = {A_cl:.6f}", flush=True)
    print(f"  Total canonical mass Σa_n = {np.sum(a_can):.6e}", flush=True)
    print(f"  Total real mass Σ|b_n|   = {np.sum(np.abs(b_off)):.6e}", flush=True)
    print(f"  Ratio real/can           = {np.sum(np.abs(b_off))/np.sum(a_can):.6e}", flush=True)

    A_can_L = np.zeros(M + 1)
    for n in range(M):
        A_can_L[n + 1] = A_can_L[n] + a_can[n]

    A_can_R = np.zeros(M + 1)
    A_can_R[M] = 0
    for n in range(M - 1, -1, -1):
        A_can_R[n] = A_can_R[n + 1] + a_can[n]

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

    print(f"\n{'m':>3}  {'gL':>10}  {'gR':>10}  "
          f"{'Λ^can':>12}  {'Θ^can':>12}  {'max_can':>12}  "
          f"{'Λ^off':>12}  {'Θ^off':>12}  {'max_off':>12}  "
          f"{'can/off':>10}  {'can<1':>5}", flush=True)
    print('-' * 130, flush=True)

    results = []
    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        d_m = s_m - w_min
        gL = s_m - entry['u_left']
        gR = entry['u_right'] - s_m

        Lam_can = d_m * abs(A_can_L[m - 1]) / gL if gL > 1e-50 else float('inf')
        The_can = d_m * abs(A_can_R[m]) / gR if gR > 1e-50 else float('inf')
        max_can = max(Lam_can, The_can)

        Lam_off = d_m * abs(B_cum[m - 1]) / gL if gL > 1e-50 else float('inf')
        The_off = d_m * abs(C_cum[m]) / gR if gR > 1e-50 else float('inf')
        max_off = max(Lam_off, The_off)

        ratio = max_can / max_off if max_off > 1e-50 else float('inf')
        ok = '✓' if max_can < 1 else '✗'

        results.append({
            'm': m, 'Lam_can': Lam_can, 'The_can': The_can, 'max_can': max_can,
            'Lam_off': Lam_off, 'The_off': The_off, 'max_off': max_off,
            'ratio': ratio, 'gL': gL, 'gR': gR, 'd': d_m,
        })

        if m <= 5 or m >= M - 2 or max_can >= 1 or max_off >= 0.01:
            print(f"{m:3d}  {gL:10.2e}  {gR:10.2e}  "
                  f"{Lam_can:12.4e}  {The_can:12.4e}  {max_can:12.4e}  "
                  f"{Lam_off:12.4e}  {The_off:12.4e}  {max_off:12.4e}  "
                  f"{ratio:10.2e}  {ok:>5}", flush=True)
        elif m == 6:
            print("  ...", flush=True)

    n_can_ok = sum(1 for r in results if r['max_can'] < 1)
    n_can_fail = sum(1 for r in results if r['max_can'] >= 1)
    max_can_all = max(r['max_can'] for r in results)
    max_off_all = max(r['max_off'] for r in results)
    max_ratio = max(r['ratio'] for r in results)

    print(f"\n--- Ergebnis λ={lam} ---", flush=True)
    print(f"  Canonical < 1:  {n_can_ok}/{len(results)}", flush=True)
    print(f"  Canonical >= 1: {n_can_fail}/{len(results)}  ← C2as.3 VERSAGT hier", flush=True)
    print(f"  max(Λ^can,Θ^can) = {max_can_all:.6e}", flush=True)
    print(f"  max(Λ^off,Θ^off) = {max_off_all:.6e}", flush=True)
    print(f"  max Ratio can/off = {max_ratio:.6e}", flush=True)

    return {
        'n_can_ok': n_can_ok, 'n_can_fail': n_can_fail,
        'max_can': max_can_all, 'max_off': max_off_all,
        'max_ratio': max_ratio, 'n_total': len(results),
    }


if __name__ == "__main__":
    for cfg in CONFIGS:
        canonical_channel_test(cfg["lam"], cfg["N"])
    print("\nDone.", flush=True)
