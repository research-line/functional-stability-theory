"""
c2_n_scaling.py — N-Skalierung der C2/MS2-Approximation

Fixiert lambda=3, variiert N=15,30,45,60.
Frage: Verbessert sich der Overlap k_lambda <-> xi_lambda bei groesserem N?
Bricht die Cluster-Degenereszenz auf?

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
    RIEMANN_ZEROS,
)
from c2_approximation_test import (
    k_lambda_value,
    h_educated_guess,
    inner_product,
    norm,
    evaluate_F_even,
)

DPS = int(os.environ.get("DPS", 50))
LAM = float(os.environ.get("LAM", 3.0))
N_VALUES = [15, 30, 45, 60]


def project_k_onto_even_basis(lam, N_val, dps_val):
    """Fourier-Projektion von k_lambda auf Even-Basis der Dimension N+1."""
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


def run_for_N(N_val, lam, dps_val):
    mp.mp.dps = dps_val
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    dim = N_val + 1

    print(f"\n{'='*80}")
    print(f"  N = {N_val}, lambda = {lam}, dim = {dim}")
    print(f"{'='*80}")

    # QW + Diagonalisierung
    print(f"  Phase A: QW-Matrix...")
    t0 = time.time()
    M_full, _ = build_QW_mp(N_val, lam, chi_func=chi_trivial, q_mod=1,
                            include_W02=True, sign_WR=-1, verbose=False)
    M_even, U_even = project_to_parity(M_full, N_val, parity="even")
    w_sorted, Q_sorted = diagonalize_mp(M_even, verbose=False)
    dt_qw = time.time() - t0
    print(f"    QW in {dt_qw:.1f}s")

    w_min = w_sorted[0]

    # Cluster
    cluster_tol = mp.mpf("1e-10")
    cluster_idx = [k for k in range(dim) if abs(w_sorted[k] - w_min) < cluster_tol]
    cluster_size = len(cluster_idx)
    cluster_edge = cluster_idx[-1]

    # Cluster-Gap: Abstand des ersten Nicht-Cluster-EW zu w_min
    if cluster_size < dim:
        gap = w_sorted[cluster_size] - w_min
    else:
        gap = mp.mpf(0)

    xi_edge = mp.matrix(dim, 1)
    for i in range(dim):
        xi_edge[i, 0] = Q_sorted[i, cluster_edge]

    print(f"    w_min = {mp.nstr(w_min, 10)}, Cluster = {cluster_size}, Gap = {mp.nstr(gap, 4)}")

    # k_lambda Fourier-Projektion
    print(f"  Phase B: k_lambda Projektion...")
    t0 = time.time()
    k_coeffs, _ = project_k_onto_even_basis(lam, N_val, dps_val)
    dt_k = time.time() - t0
    print(f"    Projektion in {dt_k:.1f}s")

    k_norm_val = norm(k_coeffs, dim)
    if k_norm_val < mp.mpf("1e-30"):
        print("    WARNUNG: k_lambda ~ 0")
        return None

    k_normed = mp.matrix(dim, 1)
    for i in range(dim):
        k_normed[i, 0] = k_coeffs[i, 0] / k_norm_val

    # Overlap mit jedem Cluster-EV
    best_k = 0
    best_ov = mp.mpf(0)
    cluster_proj_sq = mp.mpf(0)
    for k in cluster_idx:
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]
        ov = abs(inner_product(k_normed, qk, dim))
        cluster_proj_sq += ov**2
        if ov > best_ov:
            best_ov = ov
            best_k = k

    cluster_proj = mp.sqrt(cluster_proj_sq)
    print(f"    Cluster-Projektion ||P_cluster k|| = {mp.nstr(cluster_proj, 8)}")
    print(f"    Best single EV overlap: k={best_k}, |ov| = {mp.nstr(best_ov, 8)}")

    # Rayleigh
    QWk = M_even * k_normed
    eps_rayleigh = inner_product(k_normed, QWk, dim)

    # Residuum
    res = mp.matrix(dim, 1)
    for i in range(dim):
        res[i, 0] = QWk[i, 0] - eps_rayleigh * k_normed[i, 0]
    res_norm = norm(res, dim)

    print(f"    Rayleigh eps = {mp.nstr(eps_rayleigh, 10)}, Diff = {mp.nstr(eps_rayleigh - w_min, 4)}")
    print(f"    Residuum ||QW*k - eps*k|| = {mp.nstr(res_norm, 4)}")

    # F an Riemann-Nullstellen
    print(f"\n    F(gamma) an Riemann-Nullstellen:")
    f_k_zeros = []
    f_xi_zeros = []
    for gamma in RIEMANN_ZEROS:
        F_k = evaluate_F_even(k_coeffs, N_val, L_mp, gamma)
        F_xi = evaluate_F_even(xi_edge, N_val, L_mp, gamma)
        f_k_zeros.append(abs(float(F_k)))
        f_xi_zeros.append(abs(float(F_xi)))
        print(f"      gamma={gamma:.6f}: |F_k|={abs(float(F_k)):.3e}, |F_xi|={abs(float(F_xi)):.3e}")

    # Optimalen Cluster-Vektor fuer k finden: Projektion von k auf Cluster
    # und F davon testen
    print(f"\n    F der Cluster-Projektion von k:")
    k_cluster = mp.matrix(dim, 1)
    for k in cluster_idx:
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]
        ov = inner_product(k_normed, qk, dim)
        for i in range(dim):
            k_cluster[i, 0] += ov * qk[i, 0]
    # k_cluster sollte ≈ k_normed sein (da 100% im Cluster)
    k_cl_norm = norm(k_cluster, dim)
    # Unnormiert lassen, da k_coeffs auch unnormiert in F eingeht
    k_cl_coeffs = mp.matrix(dim, 1)
    for i in range(dim):
        k_cl_coeffs[i, 0] = k_cluster[i, 0] * k_norm_val
    for gamma in RIEMANN_ZEROS[:3]:
        F_cl = evaluate_F_even(k_cl_coeffs, N_val, L_mp, gamma)
        print(f"      gamma={gamma:.6f}: |F_cluster_proj|={abs(float(F_cl)):.3e}")

    return {
        "N": N_val,
        "w_min": float(w_min),
        "cluster_size": cluster_size,
        "cluster_gap": float(gap),
        "cluster_proj": float(cluster_proj),
        "best_k": best_k,
        "best_overlap": float(best_ov),
        "rayleigh": float(eps_rayleigh),
        "residual": float(res_norm),
        "F_k_gamma1": f_k_zeros[0] if f_k_zeros else 0,
        "F_xi_gamma1": f_xi_zeros[0] if f_xi_zeros else 0,
    }


def main():
    mp.mp.dps = DPS
    print("=" * 80)
    print(f"C2/MS2: N-SKALIERUNG bei lambda={LAM}")
    print("=" * 80)
    print(f"  N-Werte: {N_VALUES}")
    print(f"  DPS = {DPS}")

    results = []
    for N_val in N_VALUES:
        r = run_for_N(N_val, LAM, DPS)
        if r:
            results.append(r)

    print(f"\n{'='*80}")
    print(f"  N-SKALIERUNG ZUSAMMENFASSUNG (lambda={LAM})")
    print(f"{'='*80}")
    print(f"  {'N':>4s}  {'Cluster':>8s}  {'Gap':>10s}  {'Cl.Proj':>8s}  {'BestOv':>8s}  {'best_k':>7s}  {'|F_k(g1)|':>12s}  {'|F_xi(g1)|':>12s}  {'Ratio':>10s}")
    for r in results:
        ratio = r["F_k_gamma1"] / r["F_xi_gamma1"] if r["F_xi_gamma1"] > 0 else float("inf")
        print(f"  {r['N']:4d}  {r['cluster_size']:8d}  {r['cluster_gap']:10.2e}  {r['cluster_proj']:8.5f}  {r['best_overlap']:8.5f}  {r['best_k']:7d}  {r['F_k_gamma1']:12.3e}  {r['F_xi_gamma1']:12.3e}  {ratio:10.1e}")


if __name__ == "__main__":
    main()
