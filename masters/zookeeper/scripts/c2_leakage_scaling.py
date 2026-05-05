"""
c2_leakage_scaling.py — Lambda-Skalierung der Out-of-Cluster-Leckage

Kernfrage: Verschwindet die Nicht-Cluster-Komponente von k_lambda im Limes lambda->inf?

Messgroessen:
  - epsilon = 1 - ||P_cluster k||^2  (Out-of-Cluster-Anteil)
  - |F_raw(gamma)| vs |F_proj(gamma)|  (F mit/ohne Cluster-Projektion)
  - Cluster-Groesse vs. Dimension

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
    inner_product,
    norm,
    evaluate_F_even,
)

DPS = int(os.environ.get("DPS", 50))
N = int(os.environ.get("N", 30))
LAMBDAS = [3.0, 5.0, 7.0, 10.0]


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


def run_leakage_analysis(lam, N_val, dps_val):
    mp.mp.dps = dps_val
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    dim = N_val + 1

    print(f"\n{'='*80}")
    print(f"  LAMBDA = {lam}, N = {N_val}, dim = {dim}")
    print(f"{'='*80}")

    # QW
    t0 = time.time()
    M_full, _ = build_QW_mp(N_val, lam, chi_func=chi_trivial, q_mod=1,
                            include_W02=True, sign_WR=-1, verbose=False)
    M_even, U_even = project_to_parity(M_full, N_val, parity="even")
    w_sorted, Q_sorted = diagonalize_mp(M_even, verbose=False)
    dt_qw = time.time() - t0

    w_min = w_sorted[0]

    # Cluster
    cluster_tol = mp.mpf("1e-10")
    cluster_idx = [k for k in range(dim) if abs(w_sorted[k] - w_min) < cluster_tol]
    cluster_size = len(cluster_idx)

    if cluster_size < dim:
        gap = w_sorted[cluster_size] - w_min
    else:
        gap = mp.mpf(0)

    print(f"  QW in {dt_qw:.1f}s, w_min={mp.nstr(w_min, 8)}, Cluster={cluster_size}/{dim}, Gap={mp.nstr(gap, 3)}")

    # k_lambda
    t0 = time.time()
    k_coeffs, _ = project_k_onto_even_basis(lam, N_val, dps_val)
    dt_k = time.time() - t0

    k_norm_val = norm(k_coeffs, dim)
    if k_norm_val < mp.mpf("1e-30"):
        print("  WARNUNG: k_lambda ~ 0")
        return None

    k_normed = mp.matrix(dim, 1)
    for i in range(dim):
        k_normed[i, 0] = k_coeffs[i, 0] / k_norm_val

    print(f"  k_lambda projiziert in {dt_k:.1f}s, ||k|| = {mp.nstr(k_norm_val, 6)}")

    # Cluster-Projektion
    cluster_proj_sq = mp.mpf(0)
    k_cluster = mp.matrix(dim, 1)
    overlaps = []
    for k in cluster_idx:
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]
        ov = inner_product(k_normed, qk, dim)
        overlaps.append(float(abs(ov)))
        cluster_proj_sq += ov**2
        for i in range(dim):
            k_cluster[i, 0] += ov * qk[i, 0]

    cluster_proj = mp.sqrt(cluster_proj_sq)
    leakage = 1 - float(cluster_proj_sq)

    # Nicht-Cluster-Projektion (Orthogonalkomplement)
    k_perp = mp.matrix(dim, 1)
    for i in range(dim):
        k_perp[i, 0] = k_normed[i, 0] - k_cluster[i, 0]
    perp_norm = norm(k_perp, dim)

    print(f"\n  CLUSTER-ANALYSE:")
    print(f"    ||P_cluster k||     = {mp.nstr(cluster_proj, 12)}")
    print(f"    ||P_cluster k||^2   = {mp.nstr(cluster_proj_sq, 12)}")
    print(f"    Leckage eps         = {leakage:.4e}")
    print(f"    ||k_perp||          = {mp.nstr(perp_norm, 8)}")

    # Nicht-Cluster-Overlap-Verteilung
    print(f"\n  Nicht-Cluster-Overlaps (erste 5):")
    non_cluster_ovs = []
    for k in range(cluster_size, min(cluster_size + 5, dim)):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]
        ov = abs(inner_product(k_normed, qk, dim))
        non_cluster_ovs.append(float(ov))
        print(f"    k={k}: |ov|={mp.nstr(ov, 6)}, w={mp.nstr(w_sorted[k], 6)}")

    # F-Auswertung: roh vs. Cluster-Projektion vs. Perp-Komponente
    print(f"\n  F(gamma) ZERLEGUNG (roh = cluster + perp):")
    print(f"  {'gamma':>12s}  {'|F_raw|':>12s}  {'|F_cluster|':>12s}  {'|F_perp|':>12s}  {'F_perp/F_raw':>12s}")

    k_cl_coeffs = mp.matrix(dim, 1)
    k_perp_coeffs = mp.matrix(dim, 1)
    for i in range(dim):
        k_cl_coeffs[i, 0] = k_cluster[i, 0] * k_norm_val
        k_perp_coeffs[i, 0] = k_perp[i, 0] * k_norm_val

    f_raw_list = []
    f_cl_list = []
    f_perp_list = []

    for gamma in RIEMANN_ZEROS:
        F_raw = evaluate_F_even(k_coeffs, N_val, L_mp, gamma)
        F_cl = evaluate_F_even(k_cl_coeffs, N_val, L_mp, gamma)
        F_perp = evaluate_F_even(k_perp_coeffs, N_val, L_mp, gamma)

        f_raw_abs = abs(float(F_raw))
        f_cl_abs = abs(float(F_cl))
        f_perp_abs = abs(float(F_perp))

        ratio = f_perp_abs / f_raw_abs if f_raw_abs > 0 else 0
        f_raw_list.append(f_raw_abs)
        f_cl_list.append(f_cl_abs)
        f_perp_list.append(f_perp_abs)

        print(f"  {gamma:12.6f}  {f_raw_abs:12.3e}  {f_cl_abs:12.3e}  {f_perp_abs:12.3e}  {ratio:12.4f}")

    # F an Nicht-Nullstellen
    non_zeros = [10.0, 17.5, 23.0]
    print(f"\n  F(z) an Nicht-Nullstellen:")
    print(f"  {'z':>12s}  {'|F_raw|':>12s}  {'|F_cluster|':>12s}  {'|F_perp|':>12s}")
    for z in non_zeros:
        F_raw = evaluate_F_even(k_coeffs, N_val, L_mp, z)
        F_cl = evaluate_F_even(k_cl_coeffs, N_val, L_mp, z)
        F_perp = evaluate_F_even(k_perp_coeffs, N_val, L_mp, z)
        print(f"  {z:12.6f}  {abs(float(F_raw)):12.3e}  {abs(float(F_cl)):12.3e}  {abs(float(F_perp)):12.3e}")

    return {
        "lambda": lam,
        "cluster_size": cluster_size,
        "dim": dim,
        "gap": float(gap),
        "leakage": leakage,
        "perp_norm": float(perp_norm),
        "F_raw_g1": f_raw_list[0] if f_raw_list else 0,
        "F_cl_g1": f_cl_list[0] if f_cl_list else 0,
        "F_perp_g1": f_perp_list[0] if f_perp_list else 0,
    }


def main():
    mp.mp.dps = DPS
    print("=" * 80)
    print("C2/MS2: LAMBDA-SKALIERUNG DER OUT-OF-CLUSTER-LECKAGE")
    print("=" * 80)
    print(f"  N = {N}, DPS = {DPS}")
    print(f"  Lambdas: {LAMBDAS}")

    results = []
    for lam in LAMBDAS:
        r = run_leakage_analysis(lam, N, DPS)
        if r:
            results.append(r)

    print(f"\n{'='*80}")
    print(f"  LECKAGE-SKALIERUNG ZUSAMMENFASSUNG")
    print(f"{'='*80}")
    print(f"  {'lam':>6s}  {'Cl/Dim':>8s}  {'Gap':>10s}  {'Leckage':>12s}  {'||k_perp||':>12s}  {'|F_raw(g1)|':>12s}  {'|F_cl(g1)|':>12s}  {'|F_perp(g1)|':>12s}")
    for r in results:
        cl_dim = f"{r['cluster_size']}/{r['dim']}"
        print(f"  {r['lambda']:6.1f}  {cl_dim:>8s}  {r['gap']:10.2e}  {r['leakage']:12.4e}  {r['perp_norm']:12.4e}  {r['F_raw_g1']:12.3e}  {r['F_cl_g1']:12.3e}  {r['F_perp_g1']:12.3e}")

    # Skalierungs-Fit: leakage ~ lambda^alpha?
    if len(results) >= 3:
        lams = [r["lambda"] for r in results]
        leaks = [r["leakage"] for r in results]
        valid = [(l, e) for l, e in zip(lams, leaks) if e > 0]
        if len(valid) >= 2:
            log_lams = [np.log(l) for l, _ in valid]
            log_leaks = [np.log(e) for _, e in valid]
            # Linearer Fit: log(leak) = alpha * log(lam) + const
            n_pts = len(log_lams)
            mean_x = sum(log_lams) / n_pts
            mean_y = sum(log_leaks) / n_pts
            ss_xx = sum((x - mean_x)**2 for x in log_lams)
            ss_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(log_lams, log_leaks))
            alpha = ss_xy / ss_xx if ss_xx > 0 else 0
            print(f"\n  Skalierungs-Exponent: Leckage ~ lambda^alpha, alpha = {alpha:.2f}")
            print(f"  (alpha < 0 => Leckage verschwindet fuer lambda -> inf)")


if __name__ == "__main__":
    main()
