"""
c2br_contraction_test.py

Kontraktionsargument fuer die Abel-Schließung.

Schlüsselidee:
  δ_j = (1 - 2 A_cl) + correction(k_perp)

  Der kanonische Teil δ^can = 1 - 2A_cl ist NICHT-ZIRKULAR (folgt aus C2v).
  |δ^can| < 1 mit Marge M = 1 - |1-2A_cl| = min(2A_cl, 2-2A_cl).

  Die Korrektur kommt aus k_perp und ist via CS beschränkt:
  |correction_j| <= (d_j/α) × ||k_perp|| × sqrt(Σ_{off} α_a²/(w_a-t_j)²)

  Definiere: C_j = (d_j/α) × sqrt(1/f_j² - A_cl/d_j²)
  Dann: |correction_j| <= C_j × ||k_perp||

  Falls C_max := max_j C_j < M, dann:
  ||k_perp|| <= 1 (trivial) => |δ_j| <= |δ^can| + C_max < 1
  => Sign Coherence => C2h => MS2 => ||k_perp|| -> 0

  NICHT-ZIRKULÄRER BEWEIS! (wenn C_max < M)
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import mpmath as mp
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, chi_trivial,
)
from c2_approximation_test import (
    h_educated_guess, k_lambda_value, norm,
)
from c2_poisson_decomposition import (
    h_hat_analytical, project_to_fourier,
)

DPS = int(os.environ.get("DPS", 25))


def to_np(M, dim):
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = float(M[i, j])
    return A


def compute_contraction(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    A_np = to_np(Aq, dim)
    H_np = to_np(Ah, dim)
    W_np = A_np - H_np

    col0 = W_np[:, 0]
    cn = np.linalg.norm(col0)
    ut = col0 / cn

    P_perp = np.eye(dim) - np.outer(ut, ut)

    ws_A, vs_A = np.linalg.eigh(A_np)
    w0 = ws_A[0]

    # T = P_perp H P_perp
    T_np = P_perp @ H_np @ P_perp
    ts_all, vs_T_all = np.linalg.eigh(T_np)

    # Cluster in A
    cluster_tol = 1e-6
    cl_mask = np.abs(ws_A - w0) < cluster_tol
    off_mask = ~cl_mask
    cl_idx = np.where(cl_mask)[0]
    off_idx = np.where(off_mask)[0]

    A_cl = sum(float(vs_A[:, i] @ ut)**2 for i in cl_idx)

    # Off-cluster A-Eigenwerte
    w_off = ws_A[off_idx]
    v_off = vs_A[:, off_idx]
    off_sort = np.argsort(w_off)
    w_off = w_off[off_sort]
    v_off = v_off[:, off_sort]
    n_off = len(w_off)

    alpha_a = np.array([float(v_off[:, a] @ ut) for a in range(n_off)])

    # T-Eigenwerte (nicht-trivial)
    t_mask = np.abs(ts_all) > 1e-10
    t_idx = np.where(t_mask)[0]
    ts = ts_all[t_idx]
    vs_T = vs_T_all[:, t_idx]
    t_sort = np.argsort(ts)
    ts = ts[t_sort]
    vs_T = vs_T[:, t_sort]
    n_T = len(ts)

    # f_j berechnen
    uHu = float(ut @ H_np @ ut)
    f_vec = H_np @ ut - uHu * ut
    f_j = np.array([float(vs_T[:, j] @ f_vec) for j in range(n_T)])

    # Poisson-Vektor
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
    alpha = float(ut @ kn)

    # k_perp
    k_perp = kn - alpha * ut
    # Projiziere auf off-cluster A-Eigenraum
    c_a_off = np.array([float(v_off[:, a] @ kn) for a in range(n_off)])
    k_perp_sq = np.sum(c_a_off**2)
    k_perp_norm = np.sqrt(k_perp_sq)

    # mu/alpha
    mu_alpha = float(ut @ (A_np - w0 * np.eye(dim)) @ kn) / alpha

    # === KONTRAKTIONSANALYSE ===

    # Kanonisches delta
    delta_can = 1.0 - 2.0 * A_cl
    margin = 1.0 - abs(delta_can)

    # C_j fuer jeden T-Eigenwert
    C_vals = np.zeros(n_T)
    delta_actual = np.zeros(n_T)
    correction_actual = np.zeros(n_T)

    for j in range(n_T):
        t_j = ts[j]
        d_j = t_j - w0

        if abs(d_j) < 1e-15 or abs(f_j[j]) < 1e-30:
            continue

        # Off-cluster spektrale Summe: Σ alpha_a^2/(w_a-t_j)^2
        off_spec_sum = np.sum(alpha_a**2 / (w_off - t_j)**2)

        # Cluster-Beitrag: A_cl/d_j^2
        cl_spec = A_cl / d_j**2

        # Off-cluster-Beitrag (sollte = 1/f_j^2 - A_cl/d_j^2)
        total_spec = 1.0 / f_j[j]**2
        off_spec_from_f = total_spec - cl_spec

        # C_j = (d_j/alpha) * sqrt(off_spec_sum)
        if off_spec_sum > 0:
            C_vals[j] = (d_j / alpha) * np.sqrt(off_spec_sum)

        # Tatsächliches delta (direkt berechnet)
        cauchy_sum_off = np.sum(c_a_off * alpha_a / (w_off - t_j))
        eta_j_from_off = -(d_j / alpha) * cauchy_sum_off

        # alpha_cl = Σ_{Cl} c_a alpha_a
        alpha_cl_actual = sum(float(vs_A[:, i] @ kn) * float(vs_A[:, i] @ ut) for i in cl_idx)
        eta_j_cl = alpha_cl_actual / alpha

        delta_actual[j] = 1.0 - (eta_j_cl + eta_j_from_off)
        correction_actual[j] = delta_actual[j] - delta_can

    return {
        'lam': lam, 'N': N, 'n_off': n_off, 'n_T': n_T,
        'w0': w0, 'alpha': alpha, 'A_cl': A_cl,
        'mu_alpha': mu_alpha, 'k_perp_norm': k_perp_norm,
        'delta_can': delta_can, 'margin': margin,
        'C_vals': C_vals, 'C_max': np.max(C_vals),
        'delta_actual': delta_actual,
        'correction_actual': correction_actual,
        'ts': ts, 'f_j': f_j,
    }


def main():
    print("C2br: Kontraktionsargument")
    print(f"DPS={DPS}")
    print()

    all_results = []
    for lam in [3.0, 5.0, 7.0, 9.0]:
        for N in [15, 20, 25]:
            t0 = time.time()
            d = compute_contraction(lam, N)
            elapsed = time.time() - t0
            all_results.append(d)

            print(f"lam={lam}, N={N} ({elapsed:.0f}s):")
            print(f"  A_cl={d['A_cl']:.6f}, alpha={d['alpha']:.6f}")
            print(f"  delta^can = 1-2A_cl = {d['delta_can']:+.6f}")
            print(f"  |delta^can| = {abs(d['delta_can']):.6f}")
            print(f"  Marge M = 1-|d^can| = {d['margin']:.6f}")
            print(f"  C_max = {d['C_max']:.6f}")
            print(f"  C_max < M? {'JA => KONTRAKTION!' if d['C_max'] < d['margin'] else 'NEIN'}")
            print(f"  ||k_perp|| = {d['k_perp_norm']:.6e}")
            print(f"  C_max * ||k_perp|| = {d['C_max'] * d['k_perp_norm']:.6e}")
            print(f"  |delta^can| + C_max*||k_perp|| = {abs(d['delta_can']) + d['C_max']*d['k_perp_norm']:.6f}")
            print(f"  max |delta_actual| = {np.max(np.abs(d['delta_actual'])):.6e}")
            print(f"  max |correction| = {np.max(np.abs(d['correction_actual'])):.6e}")

            # Modeweise Analyse (Top-3 C_j)
            top3 = np.argsort(d['C_vals'])[-3:][::-1]
            cvals = d['C_vals']
            parts = [f"j={j}: C={cvals[j]:.4f}" for j in top3]
            print(f"  Top-3 C_j: {', '.join(parts)}")
            print()

    # Zusammenfassung
    print("=" * 90)
    print("ZUSAMMENFASSUNG: KONTRAKTIONSARGUMENT")
    print("=" * 90)
    print(f"{'lam':>5} {'N':>4} {'A_cl':>8} {'|d^can|':>8} {'Marge':>8} {'C_max':>8} {'C<M?':>6} {'eps_0':>10}")
    for d in all_results:
        ok = "JA" if d['C_max'] < d['margin'] else "NEIN"
        eps0 = d['margin'] / d['C_max'] if d['C_max'] > 1e-10 else float('inf')
        print(f"{d['lam']:5.1f} {d['N']:4d} {d['A_cl']:8.4f} {abs(d['delta_can']):8.4f} "
              f"{d['margin']:8.4f} {d['C_max']:8.4f} {ok:>6} {eps0:10.4f}")

    print()
    print("Legende:")
    print("  delta^can = 1 - 2*A_cl (kanonisch, aus C2v, NICHT-ZIRKULAER)")
    print("  Marge = 1 - |delta^can| (verfügbarer Spielraum)")
    print("  C_max = max_j (d_j/alpha) * sqrt(off-cluster spectral sum)")
    print("  C_max < M => Triviale Kontraktion: ||k_perp||<=1 => |delta|<1 => MS2")
    print("  eps_0 = M/C_max (kritischer Startpunkt)")
    print("  Wenn eps_0 > 1: NICHT-ZIRKULAERER BEWEIS (triviale Schranke genuegt)")
    print("  Wenn eps_0 < 1 aber > ||k_perp||: Kontraktion mit Anfangsbedingung")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
