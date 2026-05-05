"""
c2ay_peeled_moment_test.py — Peeled Moment Bound

Kernidee: C_m^off = Σ_{n≥m} b_n ≤ (1/Δ_m)(μ/α - Σ_{n<m} Δ_n b_n)

Besser als roher Markov (μ/α)/Δ_m, weil die bereits verbrauchten
frühen Momente Δ_n b_n abgezogen werden.

Ziel: Schließe m=2 (einziger Markov-Ausreißer).
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from c2q_first_shell_dominance import build_operators

DPS = int(os.environ.get("DPS", 50))


def peeled_moment_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*120}", flush=True)
    print(f"C2AY Peeled Moment Bound: λ={lam}, N={N}", flush=True)
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

    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    M = len(noncl_sorted)

    b = np.zeros(M)
    Delta = np.zeros(M)
    beta = np.zeros(M)
    for idx, a in enumerate(noncl_sorted):
        b[idx] = c_arr[a] * alpha_a[a] / alpha
        Delta[idx] = w_arr[a] - w_min
        beta[idx] = alpha_a[a]**2

    mu_alpha = np.sum(Delta * b)

    # Kumulative Moment-Summe: S_m = Σ_{n<m} Δ_n b_n
    S_cum = np.zeros(M + 1)
    for n in range(M):
        S_cum[n + 1] = S_cum[n] + Delta[n] * b[n]

    # Tail-Summe: C_m = Σ_{n≥m} b_n
    C_cum = np.zeros(M + 1)
    C_cum[M] = 0
    for n in range(M - 1, -1, -1):
        C_cum[n] = C_cum[n + 1] + b[n]

    # B-Summe: B_m = Σ_{n<m} b_n
    B_cum = np.zeros(M + 1)
    for n in range(M):
        B_cum[n + 1] = B_cum[n] + b[n]

    print(f"\n--- Moment-Budget ---", flush=True)
    print(f"  μ/α = Σ Δ_n b_n = {mu_alpha:.6e}", flush=True)
    print(f"", flush=True)
    print(f"  {'n':>3}  {'Δ_n':>12}  {'b_n':>12}  {'Δ_n·b_n':>12}  "
          f"{'S_n (kum)':>12}  {'S_n/μα':>8}  {'Rest':>12}", flush=True)
    print(f"  {'-'*85}", flush=True)
    for n in range(min(10, M)):
        rest = mu_alpha - S_cum[n + 1]
        frac = S_cum[n + 1] / mu_alpha if abs(mu_alpha) > 1e-50 else 0
        print(f"  {n:3d}  {Delta[n]:12.4e}  {b[n]:+12.4e}  {Delta[n]*b[n]:+12.4e}  "
              f"{S_cum[n+1]:12.4e}  {frac:8.4f}  {rest:+12.4e}", flush=True)

    # --- Interlacing ---
    u_poles = [w_min] + [w_arr[a] for a in noncl_sorted]
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

    # --- Haupt-Tabelle: Peeled vs Markov vs Actual ---
    print(f"\n{'='*140}", flush=True)
    print(f"{'m':>3}  {'d_m':>10}  {'gR':>10}  "
          f"{'|C_m| act':>12}  {'Markov':>12}  {'Peeled':>12}  "
          f"{'Θ_act':>12}  {'Θ_Markov':>12}  {'Θ_Peeled':>12}  "
          f"{'Gain':>8}  {'Θ_P<1':>5}", flush=True)
    print('-' * 140, flush=True)

    n_markov_ok = 0
    n_peeled_ok = 0
    n_actual_ok = 0
    worst_markov = 0
    worst_peeled = 0
    worst_actual = 0

    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        d_m = s_m - w_min
        gR = entry['u_right'] - s_m
        gL = s_m - entry['u_left']
        idx = m - 1

        C_m_act = abs(C_cum[idx])

        markov = abs(mu_alpha) / Delta[idx] if Delta[idx] > 1e-50 else float('inf')

        remaining_moment = mu_alpha - S_cum[idx]
        peeled = max(remaining_moment, 0) / Delta[idx] if Delta[idx] > 1e-50 else float('inf')

        theta_act = d_m * C_m_act / gR if gR > 1e-50 else float('inf')
        theta_markov = d_m * markov / gR if gR > 1e-50 else float('inf')
        theta_peeled = d_m * peeled / gR if gR > 1e-50 else float('inf')

        gain = markov / peeled if peeled > 1e-50 else float('inf')

        if theta_markov < 1:
            n_markov_ok += 1
        if theta_peeled < 1:
            n_peeled_ok += 1
        if theta_act < 1:
            n_actual_ok += 1
        worst_markov = max(worst_markov, theta_markov)
        worst_peeled = max(worst_peeled, theta_peeled)
        worst_actual = max(worst_actual, theta_act)

        show = (m <= 5 or m >= M - 2 or theta_markov >= 0.5
                or theta_peeled >= 0.1 or gain > 5)
        if show:
            print(f"{m:3d}  {d_m:10.2e}  {gR:10.2e}  "
                  f"{C_m_act:12.4e}  {markov:12.4e}  {peeled:12.4e}  "
                  f"{theta_act:12.4e}  {theta_markov:12.4e}  {theta_peeled:12.4e}  "
                  f"{gain:8.1f}×  {'✓' if theta_peeled < 1 else '✗':>5}", flush=True)
        elif m == 6:
            print("  ...", flush=True)

    n_total = len(s_zeros)

    # --- Auch linken Kanal mit Peeled prüfen ---
    print(f"\n--- Linker Kanal Λ_m (Peeled Moment von rechts) ---", flush=True)
    # Für B_m = Σ_{n<m} b_n brauchen wir den umgekehrten Peeled:
    # B_m = Σ_{n<m} b_n. Moment: Σ_{n<m} Δ_n b_n = S_m.
    # Markov (von oben): B_m ≤ S_m / Δ_{m-1}? Nein, das geht nicht direkt.
    # B_m hat keinen einfachen Moment-Bound weil Δ_n WÄCHST.
    # Aber B_m ist direkt klein (front-loading + wenige Terme).
    worst_lambda = 0
    n_lambda_ok = 0
    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        d_m = s_m - w_min
        gL = s_m - entry['u_left']
        B_m = abs(B_cum[m - 1])
        lam_act = d_m * B_m / gL if gL > 1e-50 else float('inf')
        if lam_act < 1:
            n_lambda_ok += 1
        worst_lambda = max(worst_lambda, lam_act)
    print(f"  Λ actual: {n_lambda_ok}/{n_total} OK, worst = {worst_lambda:.4e}", flush=True)

    # --- Fokus: m=2 Detailanalyse ---
    print(f"\n{'='*80}", flush=True)
    print(f"--- FOKUS: m=2 Detailanalyse ---", flush=True)
    print(f"{'='*80}", flush=True)

    m2_entry = next((e for e in s_zeros if e['m'] == 2), None)
    if m2_entry:
        m = 2
        idx = 1
        s_m = m2_entry['s']
        d_m = s_m - w_min
        gR = m2_entry['u_right'] - s_m
        gL = s_m - m2_entry['u_left']

        print(f"  s_2       = w_min + {d_m:.6e}", flush=True)
        print(f"  u_1 (Pol) = w_min + {Delta[0]:.6e}", flush=True)
        print(f"  u_2 (Pol) = w_min + {Delta[1]:.6e}", flush=True)
        print(f"  g_L = s_2 - u_1 = {gL:.6e}", flush=True)
        print(f"  g_R = u_2 - s_2 = {gR:.6e}", flush=True)
        print(f"  d_2/g_R         = {d_m/gR:.4f}", flush=True)
        print(f"", flush=True)

        print(f"  C_2 = Σ_{{n≥2}} b_n = {C_cum[1]:.6e}  (actual)", flush=True)
        print(f"  b_1              = {b[0]:+.6e}", flush=True)
        print(f"  b_2              = {b[1]:+.6e}", flush=True)
        print(f"  C_3 = Σ_{{n≥3}} b_n = {C_cum[2]:.6e}", flush=True)
        print(f"  Also: C_2 = b_2 + C_3 = {b[1]:.6e} + {C_cum[2]:.6e} = {b[1]+C_cum[2]:.6e}", flush=True)
        print(f"", flush=True)

        print(f"  Moment-Budget:", flush=True)
        print(f"    μ/α                = {mu_alpha:.6e}", flush=True)
        print(f"    Δ_1·b_1           = {Delta[0]*b[0]:.6e}  ({Delta[0]*b[0]/mu_alpha*100:.4f}% des Moments)", flush=True)
        print(f"    Rest nach n=1     = {mu_alpha - S_cum[1]:.6e}", flush=True)
        rest_1 = mu_alpha - S_cum[1]
        print(f"    Peeled: C_2 ≤ Rest/Δ_2 = {rest_1:.6e}/{Delta[1]:.6e} = {rest_1/Delta[1]:.6e}", flush=True)
        print(f"    Markov: C_2 ≤ μα/Δ_2   = {abs(mu_alpha)/Delta[1]:.6e}", flush=True)
        print(f"    Actual: C_2             = {abs(C_cum[1]):.6e}", flush=True)
        print(f"    Gain Peeled/Markov      = {abs(mu_alpha)/max(rest_1,1e-50):.2f}×", flush=True)
        print(f"", flush=True)

        theta_markov_2 = d_m * abs(mu_alpha) / (Delta[1] * gR)
        theta_peeled_2 = d_m * max(rest_1, 0) / (Delta[1] * gR)
        theta_act_2 = d_m * abs(C_cum[1]) / gR

        print(f"  Θ_2:", flush=True)
        print(f"    Θ Markov  = {theta_markov_2:.6f}  {'✓' if theta_markov_2 < 1 else '✗ > 1'}", flush=True)
        print(f"    Θ Peeled  = {theta_peeled_2:.6f}  {'✓' if theta_peeled_2 < 1 else '✗ > 1'}", flush=True)
        print(f"    Θ Actual  = {theta_act_2:.6e}", flush=True)
        print(f"    Tightness Peeled/Actual = {theta_peeled_2/theta_act_2:.2e}", flush=True)
        print(f"", flush=True)

        # Split C_2 = b_2 + C_3, Bound C_3 per peeled ab m=3
        rest_2 = mu_alpha - S_cum[2]
        c3_peeled = max(rest_2, 0) / Delta[2] if Delta[2] > 1e-50 else float('inf')
        c2_split = abs(b[1]) + c3_peeled
        theta_split = d_m * c2_split / gR

        print(f"  Split-Bound: C_2 = |b_2| + C_3^peeled", flush=True)
        print(f"    |b_2|             = {abs(b[1]):.6e}", flush=True)
        print(f"    Rest nach n=1,2   = {rest_2:.6e}", flush=True)
        print(f"    C_3^peeled        = {c3_peeled:.6e}", flush=True)
        print(f"    |b_2| + C_3^peel  = {c2_split:.6e}", flush=True)
        print(f"    Θ Split           = {theta_split:.6f}  {'✓' if theta_split < 1 else '✗'}", flush=True)

    # --- Zusammenfassung ---
    print(f"\n{'='*80}", flush=True)
    print(f"--- ZUSAMMENFASSUNG λ={lam}, N={N} ---", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"  {'Methode':<35} {'ok':>6} {'worst':>12} {'< 1?':>6}", flush=True)
    print(f"  {'Θ actual':<35} {n_actual_ok:>3}/{n_total:>2} {worst_actual:12.4e} {'✓' if worst_actual < 1 else '✗':>6}", flush=True)
    print(f"  {'Θ Markov (μ/α / Δ_m)':<35} {n_markov_ok:>3}/{n_total:>2} {worst_markov:12.4e} {'✓' if worst_markov < 1 else '✗':>6}", flush=True)
    print(f"  {'Θ Peeled (Rest-Moment / Δ_m)':<35} {n_peeled_ok:>3}/{n_total:>2} {worst_peeled:12.4e} {'✓' if worst_peeled < 1 else '✗':>6}", flush=True)
    print(f"  {'Λ actual (links)':<35} {n_lambda_ok:>3}/{n_total:>2} {worst_lambda:12.4e} {'✓' if worst_lambda < 1 else '✗':>6}", flush=True)

    if n_peeled_ok == n_total and n_lambda_ok == n_total:
        print(f"\n  ★ HYBRID GESCHLOSSEN: max(Λ^off, Θ^peeled) < 1 auf ALLEN {n_total} Moden!", flush=True)
    else:
        fails = n_total - n_peeled_ok
        print(f"\n  Noch {fails} Θ-Moden offen.", flush=True)


if __name__ == "__main__":
    peeled_moment_test(3.0, 30)
    print("\nDone.", flush=True)
