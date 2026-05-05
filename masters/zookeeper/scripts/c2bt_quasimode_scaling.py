"""
c2bt_quasimode_scaling.py

Quantifiziert die Poisson Quasimode / Ground-State Proximity:

(1) Rayleigh-Quotient: rho_k = <k,Ak> vs w0 = min eig(A)
(2) Quasimode-Residual: ||(A - rho_k)k|| als Funktion von lambda, N
(3) Ground-State Overlap: |<k, q0>|^2, ||P_cl k||^2
(4) Quasimode-to-Projector: ||(I-P_cl)k|| <= ||(A-rho)k|| / gap
(5) Power-Law-Fits fuer alle Groessen
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
    k_lambda_value, norm,
)
from c2_poisson_decomposition import (
    project_to_fourier,
)

DPS = int(os.environ.get("DPS", 25))


def to_np(M, dim):
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = float(M[i, j])
    return A


def compute_quasimode(lam, N):
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

    # Full eigendecomposition of A
    ws, vs = np.linalg.eigh(A_np)
    w0 = ws[0]
    q0 = vs[:, 0]

    # Cluster identification
    cluster_tol = 1e-6
    cl_mask = np.abs(ws - w0) < cluster_tol
    n_cl = np.sum(cl_mask)
    cl_idx = np.where(cl_mask)[0]
    off_idx = np.where(~cl_mask)[0]

    # Cluster gap
    if len(off_idx) > 0:
        gap = np.min(np.abs(ws[off_idx] - w0))
    else:
        gap = 0.0

    # P_cl projector
    V_cl = vs[:, cl_idx]
    P_cl_k = V_cl @ (V_cl.T)  # cluster projector matrix

    # Poisson kernel k
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
    alpha = float(ut @ kn)

    # === (1) RAYLEIGH QUOTIENT ===
    Ak = A_np @ kn
    rho_k = float(kn @ Ak)  # <k, Ak> since ||k||=1
    rayleigh_gap = rho_k - w0

    # === (2) QUASIMODE RESIDUAL ===
    res_rho = Ak - rho_k * kn  # (A - rho_k)k
    res_w0 = Ak - w0 * kn      # (A - w0)k
    norm_res_rho = np.linalg.norm(res_rho)
    norm_res_w0 = np.linalg.norm(res_w0)

    # === (3) GROUND-STATE OVERLAP ===
    overlap_q0 = abs(float(kn @ q0))**2
    k_cl = P_cl_k @ kn
    overlap_cl = float(kn @ k_cl)  # ||P_cl k||^2
    leakage = 1.0 - overlap_cl
    leakage_norm = np.sqrt(max(leakage, 0))

    # Per-eigenvector overlaps (sorted by eigenvalue)
    overlaps = np.array([float(kn @ vs[:, i])**2 for i in range(dim)])

    # === (4) QUASIMODE-TO-PROJECTOR CHECK ===
    # ||(I-P_cl)k|| <= ||(A-rho)k|| / gap
    qtp_bound = norm_res_rho / gap if gap > 1e-15 else float('inf')
    qtp_actual = leakage_norm
    qtp_tight = qtp_actual / qtp_bound if qtp_bound > 1e-30 else 0

    # Also with w0 instead of rho
    qtp_bound_w0 = norm_res_w0 / gap if gap > 1e-15 else float('inf')

    # H-residual for comparison
    Hk = H_np @ kn
    res_H = Hk - w0 * kn
    norm_res_H = np.linalg.norm(res_H)

    # mu/alpha
    mu_alpha = float(ut @ (Ak - w0 * kn)) / alpha

    return {
        'lam': lam, 'N': N, 'dim': dim,
        'w0': w0, 'rho_k': rho_k,
        'rayleigh_gap': rayleigh_gap,
        'norm_res_rho': norm_res_rho,
        'norm_res_w0': norm_res_w0,
        'norm_res_H': norm_res_H,
        'overlap_q0': overlap_q0,
        'overlap_cl': overlap_cl,
        'leakage': leakage,
        'leakage_norm': leakage_norm,
        'n_cl': n_cl,
        'gap': gap,
        'qtp_bound': qtp_bound,
        'qtp_bound_w0': qtp_bound_w0,
        'qtp_actual': qtp_actual,
        'qtp_tight': qtp_tight,
        'alpha': alpha,
        'mu_alpha': mu_alpha,
        'C_tilde': cn,
        'top5_overlaps': np.sort(overlaps)[::-1][:5],
    }


def power_law_fit(x, y):
    """Fit y = C * x^alpha, return (alpha, C, r^2)."""
    mask = (np.array(x) > 0) & (np.array(y) > 0)
    if np.sum(mask) < 2:
        return 0, 0, 0
    lx = np.log(np.array(x)[mask])
    ly = np.log(np.array(y)[mask])
    n = len(lx)
    sx = np.sum(lx)
    sy = np.sum(ly)
    sxx = np.sum(lx**2)
    sxy = np.sum(lx * ly)
    denom = n * sxx - sx**2
    if abs(denom) < 1e-30:
        return 0, 0, 0
    alpha = (n * sxy - sx * sy) / denom
    logC = (sy - alpha * sx) / n
    C = np.exp(logC)
    ss_res = np.sum((ly - alpha * lx - logC)**2)
    ss_tot = np.sum((ly - np.mean(ly))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 1e-30 else 0
    return alpha, C, r2


def main():
    print("C2bt: Poisson Quasimode / Ground-State Proximity")
    print(f"DPS={DPS}")
    print()

    # Scan over lambda and N
    lambdas = [3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]
    Ns = [15, 20, 25]

    all_results = []
    for N in Ns:
        print(f"{'='*90}")
        print(f"  N = {N}")
        print(f"{'='*90}")

        for lam in lambdas:
            t0 = time.time()
            d = compute_quasimode(lam, N)
            elapsed = time.time() - t0
            all_results.append(d)

            print(f"\n  lam={lam}, N={N} ({elapsed:.0f}s):")
            print(f"    w0 = {d['w0']:.8f}, rho_k = {d['rho_k']:.8f}")
            print(f"    Rayleigh gap = rho_k - w0 = {d['rayleigh_gap']:+.6e}")
            print(f"    Cluster: {d['n_cl']} EVs, gap = {d['gap']:.6e}")
            print(f"    --- QUASIMODE RESIDUAL ---")
            print(f"    ||(A-rho)k|| = {d['norm_res_rho']:.6e}")
            print(f"    ||(A-w0)k||  = {d['norm_res_w0']:.6e}")
            print(f"    ||(H-w0)k||  = {d['norm_res_H']:.6e}")
            print(f"    --- GROUND-STATE PROXIMITY ---")
            print(f"    |<k,q0>|^2   = {d['overlap_q0']:.10f}")
            print(f"    ||P_cl k||^2 = {d['overlap_cl']:.10f}")
            print(f"    leakage      = {d['leakage']:.6e}")
            print(f"    ||k_perp||   = {d['leakage_norm']:.6e}")
            print(f"    --- QUASIMODE-TO-PROJECTOR ---")
            print(f"    Bound: ||k_perp|| <= ||(A-rho)k||/gap = {d['qtp_bound']:.6e}")
            print(f"    Actual ||k_perp|| = {d['qtp_actual']:.6e}")
            print(f"    Tightness = actual/bound = {d['qtp_tight']:.6f}")
            print(f"    mu/alpha = {d['mu_alpha']:.6e}")

    # Power-law fits per N
    print()
    print("=" * 90)
    print("POWER-LAW-SKALIERUNG")
    print("=" * 90)

    for N in Ns:
        subset = [d for d in all_results if d['N'] == N]
        lams = [d['lam'] for d in subset]

        quantities = [
            ('||(A-rho)k||', [d['norm_res_rho'] for d in subset]),
            ('||(A-w0)k||', [d['norm_res_w0'] for d in subset]),
            ('Rayleigh gap', [abs(d['rayleigh_gap']) for d in subset]),
            ('leakage ||k_perp||', [d['leakage_norm'] for d in subset]),
            ('|mu/alpha|', [abs(d['mu_alpha']) for d in subset]),
            ('gap', [d['gap'] for d in subset]),
        ]

        print(f"\n  N={N}:")
        print(f"  {'Groesse':30s} {'Exponent':>10} {'R^2':>8} {'Werte (lam=3..15)':>60}")
        for name, vals in quantities:
            alpha_fit, C_fit, r2 = power_law_fit(lams, vals)
            vals_str = "  ".join(f"{v:.2e}" for v in vals)
            print(f"  {name:30s} {alpha_fit:+10.3f} {r2:8.4f} {vals_str}")

    # Summary table
    print()
    print("=" * 90)
    print("ZUSAMMENFASSUNG: QUASIMODE-QUALITAET")
    print("=" * 90)
    print(f"{'lam':>5} {'N':>4} {'||(A-rho)k||':>14} {'Rayleigh gap':>14} "
          f"{'||k_perp||':>12} {'gap':>12} {'QTP bound':>12} {'QTP tight':>10}")
    for d in all_results:
        print(f"{d['lam']:5.1f} {d['N']:4d} {d['norm_res_rho']:14.6e} {d['rayleigh_gap']:+14.6e} "
              f"{d['leakage_norm']:12.6e} {d['gap']:12.6e} {d['qtp_bound']:12.6e} {d['qtp_tight']:10.6f}")

    print()
    print("LEGENDE:")
    print("  ||(A-rho)k||: Quasimode-Residual (Rayleigh)")
    print("  Rayleigh gap: rho_k - w0 (positiv = k hoeherer Rayleigh als GS)")
    print("  ||k_perp||: Off-Cluster-Leckage von k")
    print("  gap: Cluster-Gap (kleinster Abstand Cluster<->Off-Cluster)")
    print("  QTP bound: ||(A-rho)k||/gap (Quasimode-to-Projector)")
    print("  QTP tight: actual/bound (naeher an 1 = tighter)")
    print()
    print("  Wenn ||(A-rho)k|| ~ lam^{-sigma} mit sigma > 0: k ist Quasimode")
    print("  Wenn gap ~ lam^{-beta} mit beta < sigma: QTP bound -> 0")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
