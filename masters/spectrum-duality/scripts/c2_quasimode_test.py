"""
c2_quasimode_test.py — Quasimode-Argument fuer MS2

ChatGPT-Audit-Route:
  ||(1 - P_cl)k_lambda|| <= ||(A_lambda - rho_lambda)k_lambda|| / g_lambda

Numerische Verifikation:
  1. Residual r_lambda = ||(QW - rho)k|| fuer verschiedene rho-Definitionen
  2. Gap g_lambda = dist(rho, spec \ cluster) fuer verschiedene Cluster-Definitionen
  3. Quasimode-Bound vs. tatsaechliches ||k_perp||
  4. Lambda-Skalierung aller drei Groessen

Schluesselinsight: Der Cluster ist nahezu degeneriert (Gap ~10^{-10}).
Frage: Kann man den Cluster "breiter" definieren und einen groesseren Gap nutzen?

Autor: LG (Opus 4.7 [1M], Session 22)
Datum: 2026-04-19
"""

import os
import sys
import time
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp,
    project_to_parity,
    diagonalize_mp,
    chi_trivial,
)
from c2_approximation_test import (
    k_lambda_value,
    inner_product,
    norm,
)

DPS = int(os.environ.get("DPS", 50))
N = int(os.environ.get("N", 30))
LAM = float(os.environ.get("LAM", 3.0))


def project_k_onto_even_basis(lam, N_val, dps_val):
    mp.mp.dps = dps_val
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    inv_sqrt_L = 1 / mp.sqrt(L_mp)
    dim = N_val + 1

    def k_of_x(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(u, lam)

    coeffs = mp.matrix(dim, 1)
    def integrand_0(x):
        return k_of_x(x) * inv_sqrt_L
    coeffs[0, 0] = mp.quad(integrand_0, [0, L_mp], maxdegree=8)

    for m in range(1, dim):
        def integrand_m(x, m_val=m):
            return k_of_x(x) * mp.sqrt(mp.mpf(2)) * inv_sqrt_L * mp.cos(2 * mp.pi * m_val * x / L_mp)
        coeffs[m, 0] = mp.quad(integrand_m, [0, L_mp], maxdegree=8)

    return coeffs, L_mp


def main():
    mp.mp.dps = DPS
    lam = LAM
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    dim = N + 1

    print("=" * 80)
    print(f"QUASIMODE-ANALYSE: lambda={lam}, N={N}, dim={dim}, DPS={DPS}")
    print("=" * 80)

    # QW aufbauen + diagonalisieren
    print("\n--- Phase 1: QW + Eigensystem ---")
    M_full, _ = build_QW_mp(N, lam, chi_func=chi_trivial, q_mod=1,
                            include_W02=True, sign_WR=-1, verbose=False)
    M_even, U_even = project_to_parity(M_full, N, parity="even")
    w_sorted, Q_sorted = diagonalize_mp(M_even, verbose=False)
    w_min = w_sorted[0]

    # Vollstaendiges Eigenwertspektrum
    print(f"\n  Vollstaendiges Spektrum (alle {dim} Eigenwerte):")
    for k in range(dim):
        gap_from_min = float(w_sorted[k] - w_min)
        print(f"    k={k:3d}: w={mp.nstr(w_sorted[k], 12):>20s}  gap_from_min={gap_from_min:.4e}")

    # k_lambda
    print("\n--- Phase 2: k_lambda Fourier-Projektion ---")
    k_coeffs, _ = project_k_onto_even_basis(lam, N, DPS)
    k_norm_val = norm(k_coeffs, dim)
    k_normed = mp.matrix(dim, 1)
    for i in range(dim):
        k_normed[i, 0] = k_coeffs[i, 0] / k_norm_val

    # Overlaps mit ALLEN Eigenvektoren
    print(f"\n--- Phase 3: k_lambda in Eigenbasis ---")
    overlaps = []
    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]
        ov = inner_product(k_normed, qk, dim)
        overlaps.append(ov)

    print(f"  {'k':>4s}  {'|c_k|':>12s}  {'|c_k|^2':>12s}  {'w_k':>15s}  {'w_k - w_min':>12s}")
    cumul_sq = mp.mpf(0)
    for k in range(min(20, dim)):
        ck = abs(overlaps[k])
        ck_sq = ck**2
        cumul_sq += ck_sq
        print(f"  {k:4d}  {float(ck):12.6e}  {float(ck_sq):12.6e}  {mp.nstr(w_sorted[k], 10):>15s}  {float(w_sorted[k]-w_min):12.4e}")
    print(f"  ... (kumulative |c|^2 bis k=19: {float(cumul_sq):.12f})")

    # --- Phase 4: Quasimode-Analyse mit verschiedenen Cluster-Definitionen ---
    print(f"\n--- Phase 4: Quasimode-Bound fuer verschiedene Cluster-Toleranzen ---")
    print(f"\n  Idee: P_cl = Projektor auf alle EVs mit |w_k - w_min| < tol")
    print(f"  r = ||(QW - rho)k||, g = dist(rho, spec \\ cluster)")
    print(f"  Bound: ||P_perp k|| <= r / g")

    tolerances = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 0.5, 1.0, 2.0]

    print(f"\n  {'tol':>10s}  {'Cl.Size':>8s}  {'||P_cl k||':>12s}  {'||P_perp k||':>12s}  {'rho':>12s}  {'r=||(QW-rho)k||':>16s}  {'g=Gap':>12s}  {'r/g (Bound)':>12s}  {'Tight?':>8s}")

    for tol in tolerances:
        # Cluster fuer diese Toleranz
        cl_idx = [k for k in range(dim) if abs(float(w_sorted[k] - w_min)) < tol]
        cl_size = len(cl_idx)

        if cl_size == 0 or cl_size == dim:
            continue

        # P_cl k
        k_cl = mp.matrix(dim, 1)
        cl_proj_sq = mp.mpf(0)
        for k in cl_idx:
            ov = overlaps[k]
            cl_proj_sq += ov**2
            qk = mp.matrix(dim, 1)
            for i in range(dim):
                qk[i, 0] = Q_sorted[i, k]
            for i in range(dim):
                k_cl[i, 0] += ov * qk[i, 0]

        p_cl_norm = mp.sqrt(cl_proj_sq)
        p_perp_norm = mp.sqrt(1 - cl_proj_sq) if cl_proj_sq < 1 else mp.mpf(0)

        # rho = Rayleigh quotient auf Cluster (gewichteter Mittelwert)
        rho_weight = mp.mpf(0)
        for k in cl_idx:
            rho_weight += overlaps[k]**2 * w_sorted[k]
        rho = rho_weight / cl_proj_sq if cl_proj_sq > 0 else w_min

        # Residual r = ||(QW - rho)k||
        QWk = M_even * k_normed
        res = mp.matrix(dim, 1)
        for i in range(dim):
            res[i, 0] = QWk[i, 0] - rho * k_normed[i, 0]
        r = norm(res, dim)

        # Gap g = min distance from rho to non-cluster eigenvalues
        non_cl_idx = [k for k in range(dim) if k not in cl_idx]
        if non_cl_idx:
            g = min(abs(float(w_sorted[k] - rho)) for k in non_cl_idx)
        else:
            g = float('inf')

        # Bound
        bound = float(r) / g if g > 0 else float('inf')
        actual = float(p_perp_norm)
        tight = actual / bound if bound > 0 and bound < float('inf') else 0

        print(f"  {tol:10.0e}  {cl_size:8d}  {float(p_cl_norm):12.8f}  {actual:12.4e}  {float(rho):12.6f}  {float(r):16.6e}  {g:12.4e}  {bound:12.4e}  {tight:8.4f}")

    # --- Phase 5: Spektrale Zerlegung des Residuals ---
    print(f"\n--- Phase 5: Spektrale Zerlegung des Residuals ---")
    print(f"  ||(QW - rho)k||^2 = sum_k |c_k|^2 (w_k - rho)^2")
    print(f"  Verwende rho = Rayleigh-Quotient (global)")

    QWk = M_even * k_normed
    rho_global = inner_product(k_normed, QWk, dim)
    res_global = mp.matrix(dim, 1)
    for i in range(dim):
        res_global[i, 0] = QWk[i, 0] - rho_global * k_normed[i, 0]
    r_global = norm(res_global, dim)

    print(f"\n  rho_global = {mp.nstr(rho_global, 12)}")
    print(f"  ||res||    = {mp.nstr(r_global, 6)}")

    print(f"\n  Beitraege zum Residual (|c_k|^2 * (w_k - rho)^2):")
    print(f"  {'k':>4s}  {'|c_k|^2':>12s}  {'(w_k-rho)^2':>14s}  {'Beitrag':>14s}  {'kum.%':>8s}")
    res_sq_total = float(r_global)**2
    cumul_res = 0
    for k in range(min(25, dim)):
        ck_sq = float(overlaps[k])**2
        dw_sq = float(w_sorted[k] - rho_global)**2
        contrib = ck_sq * dw_sq
        cumul_res += contrib
        pct = 100 * cumul_res / res_sq_total if res_sq_total > 0 else 0
        if contrib > 1e-30 or k < 10:
            print(f"  {k:4d}  {ck_sq:12.4e}  {dw_sq:14.4e}  {contrib:14.4e}  {pct:8.2f}%")

    # --- Phase 6: Alternative — Cluster PLUS Uebergangszone ---
    print(f"\n--- Phase 6: Optimaler Cluster + Gap ---")
    print(f"  Suche die Cluster-Definition die den Bound minimiert:")

    best_bound = float('inf')
    best_tol = 0
    best_info = {}

    for log_tol in np.arange(-12, 1, 0.5):
        tol = 10**log_tol
        cl_idx = [k for k in range(dim) if abs(float(w_sorted[k] - w_min)) < tol]
        cl_size = len(cl_idx)
        if cl_size == 0 or cl_size == dim:
            continue

        cl_proj_sq = sum(overlaps[k]**2 for k in cl_idx)
        rho_w = sum(overlaps[k]**2 * w_sorted[k] for k in cl_idx)
        rho = rho_w / cl_proj_sq if cl_proj_sq > 0 else w_min

        # Residual analytisch: r^2 = sum |c_k|^2 (w_k - rho)^2
        r_sq = sum(float(overlaps[k])**2 * float(w_sorted[k] - rho)**2 for k in range(dim))
        r = np.sqrt(r_sq)

        non_cl = [k for k in range(dim) if k not in cl_idx]
        if non_cl:
            g = min(abs(float(w_sorted[k] - rho)) for k in non_cl)
        else:
            continue

        if g > 0:
            bound = r / g
        else:
            continue

        p_perp = float(mp.sqrt(1 - cl_proj_sq)) if float(cl_proj_sq) < 1 else 0

        if bound < best_bound and bound > 0:
            best_bound = bound
            best_tol = tol
            best_info = {
                "tol": tol, "cl_size": cl_size, "p_perp": p_perp,
                "r": r, "g": g, "bound": bound, "tight": p_perp / bound if bound > 0 else 0,
            }

    if best_info:
        bi = best_info
        print(f"  Optimaler Cluster: tol={bi['tol']:.2e}, size={bi['cl_size']}")
        print(f"    r = {bi['r']:.6e}")
        print(f"    g = {bi['g']:.6e}")
        print(f"    Bound r/g = {bi['bound']:.6e}")
        print(f"    Actual ||P_perp k|| = {bi['p_perp']:.6e}")
        print(f"    Tightness = {bi['tight']:.4f}")


if __name__ == "__main__":
    main()
