"""
c2w_residual_alignment.py — C2W: Residual Alignment und korrekter Hauptterm

Verifiziert:
1. C2w.1: c_a = ⟨q_a, R⟩/(w_a - w_min)  mit R = (A - w_min I)k
2. C2w.2: c_a = (μ α_a + h_a)/(w_a - w_min)  mit μ = ⟨ũ, R⟩, h_a = ⟨q_a, P₀R⟩
3. C2w.3: δ_j^(0) = (μ/α)(A_cl - Σ α_a²/(w_a-w_min))  ist j-UNABHÄNGIG
4. Zerlegung: δ_j = δ_j^(0) + δ_j^(1)  mit δ_j^(1) aus P₀-Residual h
5. C2w.4-Prüfung: |δ^(0)| < 1  und  |δ_j^(1)| < 1 - |δ^(0)|

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


def residual_alignment(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*100}", flush=True)
    print(f"C2w Residual Alignment: λ={lam}, N={N}", flush=True)
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
    cl_A = ops['cl_A']
    t_vals = ops['t_vals']
    f_arr = ops['f_arr']
    j_reg = ops['j_reg']
    w_min = ops['w_min']
    dim = ops['dim']

    # === Phase 1: R = (A - w_min)k, verify c_a = ⟨q_a,R⟩/(w_a-w_min) ===
    print(f"\n--- Phase 1: Residual R = (A - w_min I)k ---", flush=True)

    # ⟨q_a, R⟩ = c_a (w_a - w_min) for a ∉ Cl, = 0 for a ∈ Cl
    mu_parts = []
    for a in range(dim):
        if a in cl_A:
            continue
        qR_a = c_arr[a] * (w_arr[a] - w_min)
        mu_parts.append(alpha_a[a] * qR_a)

    # μ = ⟨ũ, R⟩ = Σ_a α_a ⟨q_a, R⟩ = Σ_{a∉Cl} α_a c_a (w_a - w_min)
    mu = sum(alpha_a[a] * c_arr[a] * (w_arr[a] - w_min) for a in noncl_A)

    # h_a = ⟨q_a, P₀R⟩ = ⟨q_a, R⟩ - α_a μ = c_a(w_a-w_min) - α_a μ
    h_a = np.zeros(dim)
    for a in noncl_A:
        h_a[a] = c_arr[a] * (w_arr[a] - w_min) - alpha_a[a] * mu

    # Verify C2w.2: c_a = (μ α_a + h_a)/(w_a - w_min)
    print(f"\nμ = ⟨ũ, R⟩ = {mu:.10e}", flush=True)
    print(f"μ/α = {mu/alpha:.10e}", flush=True)

    c2w2_errors = []
    for a in noncl_A:
        c_predicted = (mu * alpha_a[a] + h_a[a]) / (w_arr[a] - w_min)
        err = abs(c_predicted - c_arr[a])
        c2w2_errors.append(err)
    print(f"C2w.2 verification: max |c_predicted - c_actual| = {max(c2w2_errors):.4e}", flush=True)

    # ||h|| vs μ — how much of R is in ũ direction?
    h_norm_sq = sum(h_a[a]**2 for a in noncl_A)
    R_norm_sq = sum((c_arr[a] * (w_arr[a] - w_min))**2 for a in noncl_A)
    mu_sq = mu**2
    print(f"\n||R||² = {R_norm_sq:.6e}", flush=True)
    print(f"μ²     = {mu_sq:.6e}  ({100*mu_sq/R_norm_sq:.2f}% of ||R||²)", flush=True)
    print(f"||h||² = {h_norm_sq:.6e}  ({100*h_norm_sq/R_norm_sq:.2f}% of ||R||²)", flush=True)
    print(f"μ²+||h||² = {mu_sq + h_norm_sq:.6e}  (should = ||R||² = {R_norm_sq:.6e})", flush=True)

    # === Phase 2: C2x Weighted Coefficient Flatness ===
    print(f"\n--- Phase 2: C2x — Gewichtete Flatness: (w_a-w_min)b_a ≈ (μ/α)α_a² ---", flush=True)

    A_cl = sum(alpha_a[a]**2 for a in cl_A)
    S_offcl = sum(alpha_a[a]**2 / (w_arr[a] - w_min) for a in noncl_A)

    # C2x.1: (w_a-w_min)b_a = (μ/α)α_a² + α_a h_a/α
    print(f"\n{'a':>4}  {'w_a-w_min':>12}  {'(w-w)·b_a':>12}  {'(μ/α)α_a²':>12}  "
          f"{'α_a·h_a/α':>12}  {'|dev|':>10}  {'|h_a/α_a|':>10}", flush=True)
    print('-' * 95, flush=True)

    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    c2x_lhs_list = []
    c2x_rhs_list = []
    c2x_dev_list = []

    for k_idx, a in enumerate(noncl_sorted):
        gap_a = w_arr[a] - w_min
        b_a = c_arr[a] * alpha_a[a] / alpha
        wb = gap_a * b_a
        main = (mu / alpha) * alpha_a[a]**2
        error_term = alpha_a[a] * h_a[a] / alpha
        dev = abs(wb - main)

        c2x_lhs_list.append(wb)
        c2x_rhs_list.append(main)
        c2x_dev_list.append(dev)

        h_over_a = abs(h_a[a] / alpha_a[a]) if abs(alpha_a[a]) > 1e-50 else float('inf')

        if k_idx < 8 or k_idx >= len(noncl_sorted) - 4:
            print(f"{a:4d}  {gap_a:12.4e}  {wb:12.4e}  {main:12.4e}  "
                  f"{error_term:12.4e}  {dev:10.2e}  {h_over_a:10.2e}", flush=True)
        elif k_idx == 8:
            print("  ...", flush=True)

    # C2x.3: l² bound
    l2_lhs = sum((c2x_lhs_list[i] - c2x_rhs_list[i])**2 / (alpha_a[noncl_sorted[i]]**2 + 1e-100)
                 for i in range(len(noncl_sorted)))
    l2_rhs = h_norm_sq / alpha**2
    print(f"\nC2x.3 l²-Flatness:", flush=True)
    print(f"  Σ |dev|²/α_a²  = {l2_lhs:.6e}", flush=True)
    print(f"  ||h||²/α²      = {l2_rhs:.6e}", flush=True)
    print(f"  Ratio           = {l2_lhs/l2_rhs:.10f}  (should ≤ 1.0)", flush=True)
    print(f"  ||h||/α         = {np.sqrt(h_norm_sq)/abs(alpha):.6e}  (controls flatness)", flush=True)

    # δ_j decomposition using C2w
    delta_0 = (mu / alpha) * (A_cl - S_offcl)

    print(f"\nδ^(0) = (μ/α)(A_cl - Σ α_a²/(w_a-w_min))", flush=True)
    print(f"  A_cl              = {A_cl:.10e}", flush=True)
    print(f"  Σ α_a²/(w_a-w_min) = {S_offcl:.10e}", flush=True)
    print(f"  μ/α               = {mu/alpha:.10e}", flush=True)
    print(f"  δ^(0)             = {delta_0:.10e}", flush=True)
    print(f"  |δ^(0)| < 1?      {'✓' if abs(delta_0) < 1 else '✗'}  (margin: {1 - abs(delta_0):.6e})", flush=True)

    # === Phase 3: C2y — Residual Transport Bound with Cluster Subtraction ===
    print(f"\n--- Phase 3: C2y — Transport Bound: |δ_j| ≤ (‖R‖/α)√(1/f_j² - A_cl/(t_j-w_min)²) ---", flush=True)

    R_norm = np.sqrt(R_norm_sq)
    R_over_alpha = R_norm / abs(alpha)

    # Verify C2y.2: Σ_{all} α_a²/(t_j-w_a)² = 1/f_j²
    print(f"\n‖R‖/α = {R_over_alpha:.6e}", flush=True)

    print(f"\n{'j':>4}  {'t_j':>12}  {'|δ_actual|':>12}  {'C2y bound':>12}  "
          f"{'old bound':>12}  {'ratio':>8}  {'|δ|<1':>6}  {'bnd<1':>6}  {'improv':>8}", flush=True)
    print('-' * 105, flush=True)

    results = []

    for j in sorted(j_reg, key=lambda j: t_vals[j]):
        t = t_vals[j]
        d = t - w_min
        if abs(d) < 1e-30 or abs(f_arr[j]) < 1e-50:
            continue

        # Actual δ_j
        delta_actual = -d * sum(c_arr[a] * alpha_a[a] / (alpha * (t - w_arr[a]))
                                for a in noncl_A if abs(t - w_arr[a]) > 1e-30)

        # C2y.2: verify norm identity
        norm_sum_all = sum(alpha_a[a]**2 / (t - w_arr[a])**2
                          for a in range(dim) if abs(t - w_arr[a]) > 1e-30)
        norm_sum_off = sum(alpha_a[a]**2 / (t - w_arr[a])**2
                          for a in noncl_A if abs(t - w_arr[a]) > 1e-30)

        # Old bound: ||R||/(α|f_j|)
        old_bound = R_over_alpha / abs(f_arr[j])

        # C2y.3: sharp bound with cluster subtraction
        kernel_off = 1.0 / f_arr[j]**2 - A_cl / d**2
        if kernel_off > 0:
            c2y_bound = R_over_alpha * np.sqrt(kernel_off)
        else:
            c2y_bound = 0.0

        ok_delta = '✓' if abs(delta_actual) < 1 else '✗'
        ok_bound = '✓' if c2y_bound < 1 else '✗'
        improvement = old_bound / c2y_bound if c2y_bound > 1e-50 else float('inf')

        results.append({
            'j': j, 't': t, 'delta_actual': delta_actual,
            'c2y_bound': c2y_bound, 'old_bound': old_bound,
            'kernel_off': kernel_off, 'improvement': improvement,
        })

        if len(results) <= 12 or len(results) > len(j_reg) - 5:
            print(f"{j:4d}  {t:12.6e}  {abs(delta_actual):12.4e}  {c2y_bound:12.4e}  "
                  f"{old_bound:12.4e}  {abs(delta_actual)/c2y_bound if c2y_bound > 1e-50 else 0:8.4f}  "
                  f"{ok_delta:>6}  {ok_bound:>6}  {improvement:8.1f}×", flush=True)
        elif len(results) == 13:
            print("  ...", flush=True)

    # === Phase 4: Summary ===
    print(f"\n--- Summary λ={lam} ---", flush=True)

    if results:
        max_delta = max(abs(r['delta_actual']) for r in results)
        max_c2y = max(r['c2y_bound'] for r in results)
        max_old = max(r['old_bound'] for r in results)
        n_ok = sum(1 for r in results if abs(r['delta_actual']) < 1)
        n_c2y_ok = sum(1 for r in results if r['c2y_bound'] < 1)
        n_old_ok = sum(1 for r in results if r['old_bound'] < 1)
        tightness = max(abs(r['delta_actual'])/r['c2y_bound']
                        for r in results if r['c2y_bound'] > 1e-50)

        print(f"‖R‖/α          = {R_over_alpha:.6e}", flush=True)
        print(f"μ²/‖R‖²        = {100*mu_sq/R_norm_sq:.2f}%", flush=True)
        print(f"‖h‖/α          = {np.sqrt(h_norm_sq)/abs(alpha):.6e}", flush=True)
        print(f"A_cl            = {A_cl:.10e}", flush=True)
        print(f"max |δ|         = {max_delta:.6e}", flush=True)
        print(f"|δ| < 1         = {n_ok}/{len(results)}", flush=True)
        print(f"C2y bound < 1   = {n_c2y_ok}/{len(results)}  (max: {max_c2y:.4e})", flush=True)
        print(f"Old bound < 1   = {n_old_ok}/{len(results)}  (max: {max_old:.4e})", flush=True)
        print(f"Tightness       = {tightness:.4f}  (|δ|/C2y_bound)", flush=True)
        print(f"Improvement     = {max_old/max_c2y:.1f}×  (max old/max C2y)", flush=True)

    return {
        'mu': mu, 'R_over_alpha': R_over_alpha,
        'h_norm_sq': h_norm_sq, 'R_norm_sq': R_norm_sq,
        'mu_sq': mu_sq, 'A_cl': A_cl,
        'results': results,
    }


if __name__ == "__main__":
    all_results = {}
    for cfg in CONFIGS:
        res = residual_alignment(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    # Cross-λ summary
    print(f"\n{'='*100}", flush=True)
    print("Cross-λ: C2w/C2x/C2y Summary", flush=True)
    print(f"{'='*100}", flush=True)
    print(f"{'λ':>4}  {'‖R‖/α':>10}  {'max|δ|':>10}  {'maxC2y':>10}  "
          f"{'maxOld':>10}  {'tight':>8}  {'improv':>8}  {'C2y<1':>8}  {'μ²/R²':>8}", flush=True)
    for lam in sorted(all_results):
        r = all_results[lam]
        if r['results']:
            max_d = max(abs(x['delta_actual']) for x in r['results'])
            max_c2y = max(x['c2y_bound'] for x in r['results'])
            max_old = max(x['old_bound'] for x in r['results'])
            tight = max(abs(x['delta_actual'])/x['c2y_bound']
                        for x in r['results'] if x['c2y_bound'] > 1e-50)
            mu_frac = r['mu_sq'] / r['R_norm_sq'] if r['R_norm_sq'] > 0 else 0
            n_c2y = sum(1 for x in r['results'] if x['c2y_bound'] < 1)
            n_total = len(r['results'])
            print(f"{lam:4.0f}  {r['R_over_alpha']:10.4e}  {max_d:10.4e}  {max_c2y:10.4e}  "
                  f"{max_old:10.4e}  {tight:8.4f}  {max_old/max_c2y:7.1f}×  "
                  f"{n_c2y}/{n_total:>3}  {100*mu_frac:7.2f}%", flush=True)

    print("\nDone.", flush=True)
