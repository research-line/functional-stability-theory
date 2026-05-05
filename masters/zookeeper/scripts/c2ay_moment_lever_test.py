"""
c2ay_moment_lever_test.py — Test des Momentenhebels aus C2AY

Kernidentität: μ/α = Σ_{n≥1} Δ_n b_n  (First-Moment Identity)
Markov-Bound:  C_m^off ≤ (μ/α)/Δ_m    (unter b_n ≥ 0)
Tail-β-Bound:  C_m^off ≤ (||R||/(α Δ_m)) √(Σ_{n≥m} β_n)

Testet ob der rechte Abel-Kanal Θ_m^off damit < 1 wird.
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from c2q_first_shell_dominance import build_operators

DPS = int(os.environ.get("DPS", 50))


def moment_lever_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*120}", flush=True)
    print(f"C2AY Moment Lever Test: λ={lam}, N={N}", flush=True)
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

    # --- b_n, Δ_n, β_n ---
    b = np.zeros(M)
    Delta = np.zeros(M)
    beta = np.zeros(M)
    for idx, a in enumerate(noncl_sorted):
        b[idx] = c_arr[a] * alpha_a[a] / alpha
        Delta[idx] = w_arr[a] - w_min
        beta[idx] = alpha_a[a]**2

    # --- C2ay.1: First-Moment Identity ---
    mu_over_alpha_moment = np.sum(Delta * b)

    # Residual R = (A - w_min)k, off-cluster part
    residual_sq = sum((w_arr[a] - w_min)**2 * c_arr[a]**2 for a in noncl_sorted)
    residual_norm = np.sqrt(residual_sq)

    # μ = ⟨ũ, R⟩ = Σ c_a Δ_a α_a = α Σ Δ_n b_n
    mu_direct = alpha * mu_over_alpha_moment

    print(f"\n--- First-Moment Identity (C2ay.1) ---", flush=True)
    print(f"  μ/α = Σ Δ_n b_n      = {mu_over_alpha_moment:.6e}", flush=True)
    print(f"  |μ/α|                 = {abs(mu_over_alpha_moment):.6e}", flush=True)
    print(f"  ||R||                 = {residual_norm:.6e}", flush=True)
    print(f"  α                     = {alpha:.6f}", flush=True)

    # --- Vorzeichen der b_n ---
    n_pos = np.sum(b > 0)
    n_neg = np.sum(b < 0)
    n_zero = np.sum(b == 0)
    print(f"\n--- Vorzeichen der b_n ---", flush=True)
    print(f"  b_n > 0: {n_pos}/{M}", flush=True)
    print(f"  b_n < 0: {n_neg}/{M}", flush=True)
    print(f"  b_n = 0: {n_zero}/{M}", flush=True)
    print(f"  ALLE b_n ≥ 0?  {'JA' if n_neg == 0 else 'NEIN'}", flush=True)

    if n_neg > 0:
        neg_idx = np.where(b < 0)[0]
        print(f"  Negative b_n bei Indizes: {neg_idx[:10]}...", flush=True)
        print(f"  |b_neg|/|b_pos| Ratio:  {np.sum(np.abs(b[b<0]))/np.sum(np.abs(b[b>0])):.4e}", flush=True)

    # --- Erste 10 b_n ---
    print(f"\n  Erste 10 b_n:", flush=True)
    for n in range(min(10, M)):
        print(f"    n={n:2d}  Δ_n={Delta[n]:.4e}  b_n={b[n]:+.4e}  "
              f"Δ_n·b_n={Delta[n]*b[n]:+.4e}  β_n={beta[n]:.4e}", flush=True)

    # --- Interlacing: s_m (T-Eigenwerte zwischen A-Eigenwerten) ---
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

    # --- Tail-Summen ---
    # C_m^off = Σ_{n≥m} b_n  (right tail sum)
    # B_m^off = Σ_{n<m} b_n   (left partial sum)
    B_cum = np.zeros(M + 1)
    for n in range(M):
        B_cum[n + 1] = B_cum[n] + b[n]

    C_cum = np.zeros(M + 1)
    C_cum[M] = 0
    for n in range(M - 1, -1, -1):
        C_cum[n] = C_cum[n + 1] + b[n]

    # Tail-β-Summe: Σ_{n≥m} β_n
    beta_tail = np.zeros(M + 1)
    beta_tail[M] = 0
    for n in range(M - 1, -1, -1):
        beta_tail[n] = beta_tail[n + 1] + beta[n]

    # --- Mode-für-Mode: Bounds vs. Aktuelle Werte ---
    print(f"\n{'m':>3}  {'d_m':>10}  {'gR':>10}  "
          f"{'|C_m|':>12}  {'Markov':>12}  {'Tail-β':>12}  "
          f"{'Θ_actual':>12}  {'Θ_Markov':>12}  {'Θ_tail-β':>12}  "
          f"{'Θ<1':>5}", flush=True)
    print('-' * 140, flush=True)

    n_theta_markov_ok = 0
    n_theta_tailbeta_ok = 0
    n_theta_actual_ok = 0
    worst_theta_actual = 0
    worst_theta_markov = 0
    worst_theta_tailbeta = 0
    results = []

    abs_mu_alpha = abs(mu_over_alpha_moment)

    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        d_m = s_m - w_min
        gR = entry['u_right'] - s_m
        gL = s_m - entry['u_left']

        C_m_actual = abs(C_cum[m - 1])
        markov_bound = abs_mu_alpha / Delta[m - 1] if Delta[m - 1] > 1e-50 else float('inf')
        tailbeta_bound = (residual_norm / (alpha * Delta[m - 1])) * np.sqrt(beta_tail[m - 1]) if Delta[m - 1] > 1e-50 else float('inf')

        theta_actual = d_m * C_m_actual / gR if gR > 1e-50 else float('inf')
        theta_markov = d_m * markov_bound / gR if gR > 1e-50 else float('inf')
        theta_tailbeta = d_m * tailbeta_bound / gR if gR > 1e-50 else float('inf')

        if theta_markov < 1:
            n_theta_markov_ok += 1
        if theta_tailbeta < 1:
            n_theta_tailbeta_ok += 1
        if theta_actual < 1:
            n_theta_actual_ok += 1
        worst_theta_actual = max(worst_theta_actual, theta_actual)
        worst_theta_markov = max(worst_theta_markov, theta_markov)
        worst_theta_tailbeta = max(worst_theta_tailbeta, theta_tailbeta)

        results.append({
            'm': m, 'theta_actual': theta_actual, 'theta_markov': theta_markov,
            'theta_tailbeta': theta_tailbeta, 'C_m': C_m_actual,
            'markov': markov_bound, 'tailbeta': tailbeta_bound,
            'd': d_m, 'gR': gR, 'gL': gL,
        })

        if m <= 5 or m >= M - 2 or theta_markov >= 0.5:
            print(f"{m:3d}  {d_m:10.2e}  {gR:10.2e}  "
                  f"{C_m_actual:12.4e}  {markov_bound:12.4e}  {tailbeta_bound:12.4e}  "
                  f"{theta_actual:12.4e}  {theta_markov:12.4e}  {theta_tailbeta:12.4e}  "
                  f"{'✓' if theta_actual < 1 else '✗':>5}", flush=True)
        elif m == 6:
            print("  ...", flush=True)

    n_total = len(s_zeros)

    # --- Auch Λ-Kanal (links) testen ---
    print(f"\n--- Linker Kanal Λ_m^off (zum Vergleich) ---", flush=True)
    n_lambda_ok = 0
    worst_lambda = 0
    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        d_m = s_m - w_min
        gL = s_m - entry['u_left']
        B_m = abs(B_cum[m - 1])
        lam_actual = d_m * B_m / gL if gL > 1e-50 else float('inf')
        if lam_actual < 1:
            n_lambda_ok += 1
        worst_lambda = max(worst_lambda, lam_actual)

    # --- Zusammenfassung ---
    print(f"\n{'='*80}", flush=True)
    print(f"--- ZUSAMMENFASSUNG λ={lam}, N={N} ---", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"  First-Moment μ/α     = {mu_over_alpha_moment:.6e}", flush=True)
    print(f"  ||R||                = {residual_norm:.6e}", flush=True)
    print(f"  b_n ≥ 0?            = {'JA' if n_neg == 0 else f'NEIN ({n_neg} negativ)'}", flush=True)
    print(f"", flush=True)
    print(f"  {'Methode':<30} {'ok':>6} {'worst':>12} {'< 1?':>6}", flush=True)
    print(f"  {'Θ actual':<30} {n_theta_actual_ok:>3}/{n_total:>2} {worst_theta_actual:12.4e} {'✓' if worst_theta_actual < 1 else '✗':>6}", flush=True)
    print(f"  {'Θ Markov (μ/α/Δ_m)':<30} {n_theta_markov_ok:>3}/{n_total:>2} {worst_theta_markov:12.4e} {'✓' if worst_theta_markov < 1 else '✗':>6}", flush=True)
    print(f"  {'Θ Tail-β':<30} {n_theta_tailbeta_ok:>3}/{n_total:>2} {worst_theta_tailbeta:12.4e} {'✓' if worst_theta_tailbeta < 1 else '✗':>6}", flush=True)
    print(f"  {'Λ actual (links)':<30} {n_lambda_ok:>3}/{n_total:>2} {worst_lambda:12.4e} {'✓' if worst_lambda < 1 else '✗':>6}", flush=True)
    print(f"", flush=True)

    # Ist der Markov-Bound OHNE Annahme b_n≥0 nutzbar?
    # Cauchy-Schwarz: |C_m| = |Σ_{n≥m} b_n| ≤ √M · √(Σ b_n²)
    # Besser: |C_m| ≤ Σ_{n≥m} |b_n|
    abs_tail = np.zeros(M + 1)
    abs_tail[M] = 0
    for n in range(M - 1, -1, -1):
        abs_tail[n] = abs_tail[n + 1] + abs(b[n])

    # Absoluter Markov-Bound: Σ|b_n| ≤ Σ(Δ_n/Δ_m)|b_n| ≤ (Σ Δ_n |b_n|)/Δ_m
    abs_moment = np.sum(Delta * np.abs(b))
    print(f"  Σ Δ_n |b_n|         = {abs_moment:.6e}  (absolutes Moment)", flush=True)
    print(f"  Σ Δ_n b_n           = {mu_over_alpha_moment:.6e}  (signed Moment = μ/α)", flush=True)
    print(f"  Ratio |signed|/abs  = {abs(mu_over_alpha_moment)/abs_moment:.6e}", flush=True)
    print(f"", flush=True)

    # Worst Θ mit absolutem Markov
    worst_theta_abs = 0
    n_abs_ok = 0
    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        d_m = s_m - w_min
        gR = entry['u_right'] - s_m
        abs_markov = abs_moment / Delta[m - 1] if Delta[m - 1] > 1e-50 else float('inf')
        theta_abs = d_m * abs_markov / gR if gR > 1e-50 else float('inf')
        if theta_abs < 1:
            n_abs_ok += 1
        worst_theta_abs = max(worst_theta_abs, theta_abs)

    print(f"  {'Θ abs-Markov (Σ Δ|b|/Δ_m)':<30} {n_abs_ok:>3}/{n_total:>2} {worst_theta_abs:12.4e} {'✓' if worst_theta_abs < 1 else '✗':>6}", flush=True)

    # --- Circularity Check ---
    print(f"\n--- Circularitäts-Check ---", flush=True)
    print(f"  Frage: Ist μ/α klein WEIL MS2 gilt (||k_⊥|| klein)?", flush=True)
    k_perp_sq = sum(c_arr[a]**2 for a in noncl_sorted)
    k_perp = np.sqrt(k_perp_sq)
    print(f"  ||k_⊥||             = {k_perp:.6e}", flush=True)
    print(f"  |μ/α|               = {abs_mu_alpha:.6e}", flush=True)
    print(f"  Triviale Schranke: |μ/α| ≤ ||R||/α = {residual_norm/alpha:.6e}", flush=True)
    print(f"  ||R|| = ||(A-w_min)k||_off ≤ max(Δ_n)·||k_⊥||", flush=True)
    max_Delta = np.max(Delta)
    print(f"  max(Δ_n)·||k_⊥||   = {max_Delta * k_perp:.6e}", flush=True)
    print(f"  Also: |μ/α| ≤ max(Δ)·||k_⊥||/α = {max_Delta * k_perp / alpha:.6e}", flush=True)
    print(f"  → μ/α klein ⟺ ||k_⊥|| klein ⟺ MS2", flush=True)

    return results


if __name__ == "__main__":
    moment_lever_test(3.0, 30)
    print("\nDone.", flush=True)
