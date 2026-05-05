"""
c2_variational_bound_test.py — Test ob der Variational/Extremal-Ansatz
|δ_m| < 1 erzwingen kann OHNE MS2.

Kernidee: k minimiert Rayleigh-Quotient => (A-ρ)k = R (klein).
Quasimode-Bound: ||k_⊥|| ≤ ||(A-ρ)k|| / gap.
Dann CS: |δ_m| ≤ (d_m/α) ||k_⊥|| √(Σ α_a²/(s_m-w_a)²).

ODER: Direkter Extremal-Bound: δ_m = β₀ - (d_m/α)⟨φ_m, k_⊥⟩
mit φ_m = Σ_{off} α_a q_a/(s_m - w_a).
Trivial: |⟨φ_m, k_⊥⟩| ≤ ||φ_m|| · ||k_⊥|| (CS).
Ohne MS2: ||k_⊥|| ≤ √(1-α²). Reicht das?
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from c2q_first_shell_dominance import build_operators

DPS = int(os.environ.get("DPS", 50))


def variational_test(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*100}", flush=True)
    print(f"Variational Bound Test: λ={lam}, N={N}", flush=True)
    print(f"{'='*100}", flush=True)

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
    beta0 = A_cl
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    M = len(noncl_sorted)

    k_perp_sq = sum(c_arr[a]**2 for a in noncl_sorted)
    k_perp = np.sqrt(k_perp_sq)

    Delta = np.array([w_arr[a] - w_min for a in noncl_sorted])
    gap = Delta[0]

    c_off = np.array([c_arr[a] for a in noncl_sorted])
    alpha_off = np.array([alpha_a[a] for a in noncl_sorted])

    residual_sq = sum((w_arr[a] - w_min)**2 * c_arr[a]**2 for a in noncl_sorted)
    residual = np.sqrt(residual_sq)

    quasimode_bound = residual / gap

    print(f"\n--- Variational Größen ---", flush=True)
    print(f"  α = ⟨ũ,k⟩           = {alpha:.6f}", flush=True)
    print(f"  β₀ = A_cl           = {beta0:.6f}", flush=True)
    print(f"  ||k_⊥|| (actual)    = {k_perp:.6e}", flush=True)
    print(f"  ||k_⊥|| (trivial)   = {np.sqrt(1 - alpha**2):.6e}  (√(1-α²))", flush=True)
    print(f"  gap (Δ₁)            = {gap:.6e}", flush=True)
    print(f"  ||(A-ρ)k||_off      = {residual:.6e}", flush=True)
    print(f"  Quasimode bound     = {quasimode_bound:.6e}  (residual/gap)", flush=True)
    print(f"  Quasimode / actual  = {quasimode_bound/k_perp:.2e}  ← Überschätzungsfaktor", flush=True)

    # --- Test 1: Proportionaler Fall ---
    # Wenn c_a = γ α_a, dann δ_m = (γ/α) β₀ für ALLE m
    # Worst case: |γ| = √(1-α²)/√(1-A_cl)
    gamma_max = np.sqrt(1 - alpha**2) / np.sqrt(1 - A_cl)
    delta_proportional = gamma_max / alpha * beta0

    print(f"\n--- Proportionaler Worst-Case (c_a = γ α_a) ---", flush=True)
    print(f"  |γ|_max             = {gamma_max:.6f}", flush=True)
    print(f"  |δ_m| = |γ/α| β₀   = {delta_proportional:.4f}", flush=True)
    print(f"  |δ_m| < 1?          = {'✓' if delta_proportional < 1 else '✗ VERSAGT'}", flush=True)

    # --- Test 2: Allgemeiner CS-Bound mit ||k_⊥|| = √(1-α²) ---
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
            s_zeros.append({'m': m, 's': s})
            t_matched.add(s)

    trivial_kperp = np.sqrt(1 - alpha**2)

    n_cs_trivial_ok = 0
    n_cs_actual_ok = 0
    n_cs_quasimode_ok = 0
    worst_cs_trivial = 0
    worst_cs_actual = 0
    actual_delta_max = 0

    for entry in s_zeros:
        m = entry['m']
        s_m = entry['s']
        d_m = s_m - w_min

        herglotz_slope = sum(alpha_off[n]**2 / (s_m - w_arr[noncl_sorted[n]])**2
                             for n in range(M))
        sqrt_slope = np.sqrt(herglotz_slope)

        cs_trivial = (d_m / alpha) * trivial_kperp * sqrt_slope
        cs_actual = (d_m / alpha) * k_perp * sqrt_slope
        cs_quasi = (d_m / alpha) * quasimode_bound * sqrt_slope

        actual_delta = abs(d_m / alpha * sum(
            c_off[n] * alpha_off[n] / (s_m - w_arr[noncl_sorted[n]])
            for n in range(M)))

        if cs_trivial < 1:
            n_cs_trivial_ok += 1
        if cs_actual < 1:
            n_cs_actual_ok += 1
        if cs_quasi < 1:
            n_cs_quasimode_ok += 1
        worst_cs_trivial = max(worst_cs_trivial, cs_trivial)
        worst_cs_actual = max(worst_cs_actual, cs_actual)
        actual_delta_max = max(actual_delta_max, actual_delta)

    n_total = len(s_zeros)
    print(f"\n--- CS-Bounds auf |δ_m| ---", flush=True)
    print(f"  {'Methode':<30} {'ok':>6} {'max bound':>12} {'< 1?':>6}", flush=True)
    print(f"  {'CS(||k_⊥||=√(1-α²))':<30} {n_cs_trivial_ok:>3}/{n_total:>2} {worst_cs_trivial:12.4e} {'✓' if worst_cs_trivial < 1 else '✗':>6}", flush=True)
    print(f"  {'CS(quasimode bound)':<30} {n_cs_quasimode_ok:>3}/{n_total:>2} {worst_cs_actual * quasimode_bound/k_perp:12.4e} {'✓' if n_cs_quasimode_ok == n_total else '✗':>6}", flush=True)
    print(f"  {'CS(actual ||k_⊥||)':<30} {n_cs_actual_ok:>3}/{n_total:>2} {worst_cs_actual:12.4e} {'✓' if worst_cs_actual < 1 else '✗':>6}", flush=True)
    print(f"  {'Actual |δ_m|':<30} {n_total:>3}/{n_total:>2} {actual_delta_max:12.4e} {'✓':>6}", flush=True)

    print(f"\n--- Fazit ---", flush=True)
    print(f"  Variational (trivial ||k_⊥||): {'VERSAGT' if worst_cs_trivial >= 1 else 'OK'} (bound = {worst_cs_trivial:.2e})", flush=True)
    print(f"  Variational (quasimode):       {'VERSAGT' if n_cs_quasimode_ok < n_total else 'OK'}", flush=True)
    print(f"  Mit MS2 (actual ||k_⊥||):      {'OK' if worst_cs_actual < 1 else 'VERSAGT'} (bound = {worst_cs_actual:.2e})", flush=True)
    print(f"  Lücke trivial/actual:          {worst_cs_trivial/actual_delta_max:.2e}×", flush=True)


if __name__ == "__main__":
    variational_test(3.0, 30)
    print("\nDone.", flush=True)
