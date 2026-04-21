"""
c2_quasimode_scaling.py — Quasimode-Bound Skalierung ueber lambda

Messe r (Residual), g (Gap), r/g (Bound) und ||P_perp k|| fuer lambda=3,5
mit angemessenem N, um die Skalierung der drei Groessen zu bestimmen.

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

CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
]


def project_k_onto_even_basis(lam, N_val):
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


def run_quasimode(lam, N_val):
    mp.mp.dps = DPS
    dim = N_val + 1
    L_mp = 2 * mp.log(mp.mpf(lam))

    print(f"\n{'='*80}")
    print(f"  QUASIMODE: lambda={lam}, N={N_val}, dim={dim}")
    print(f"{'='*80}")

    t0 = time.time()
    M_full, _ = build_QW_mp(N_val, lam, chi_func=chi_trivial, q_mod=1,
                            include_W02=True, sign_WR=-1, verbose=False)
    M_even, _ = project_to_parity(M_full, N_val, parity="even")
    w_sorted, Q_sorted = diagonalize_mp(M_even, verbose=False)
    w_min = w_sorted[0]
    print(f"  QW + eigsy in {time.time()-t0:.1f}s, w_min={mp.nstr(w_min, 8)}")

    t0 = time.time()
    k_coeffs, _ = project_k_onto_even_basis(lam, N_val)
    k_norm_val = norm(k_coeffs, dim)
    k_normed = mp.matrix(dim, 1)
    for i in range(dim):
        k_normed[i, 0] = k_coeffs[i, 0] / k_norm_val
    print(f"  k_lambda in {time.time()-t0:.1f}s")

    # Overlaps
    overlaps = []
    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]
        overlaps.append(inner_product(k_normed, qk, dim))

    # QW * k
    QWk = M_even * k_normed
    rho_global = inner_product(k_normed, QWk, dim)

    # Effektiver Cluster mit verschiedenen delta-Werten
    deltas = [1e-10, 1e-6, 1e-4, 1e-2, 1e-1, 0.5, 1.0]

    print(f"\n  {'delta':>10s}  {'#Cl':>5s}  {'||P_perp||':>12s}  {'Residual r':>12s}  {'Gap g':>12s}  {'r/g':>12s}  {'Tightness':>10s}")

    results = []
    for delta in deltas:
        cl_idx = [k for k in range(dim) if abs(float(w_sorted[k] - w_min)) < delta]
        cl_size = len(cl_idx)
        if cl_size == 0 or cl_size == dim:
            continue

        cl_proj_sq = sum(float(overlaps[k])**2 for k in cl_idx)
        p_perp = np.sqrt(max(0, 1 - cl_proj_sq))

        rho_w = sum(float(overlaps[k])**2 * float(w_sorted[k]) for k in cl_idx)
        rho = rho_w / cl_proj_sq if cl_proj_sq > 0 else float(w_min)

        r_sq = sum(float(overlaps[k])**2 * (float(w_sorted[k]) - rho)**2 for k in range(dim))
        r = np.sqrt(r_sq)

        non_cl = [k for k in range(dim) if k not in cl_idx]
        g = min(abs(float(w_sorted[k]) - rho) for k in non_cl) if non_cl else float('inf')

        bound = r / g if g > 0 else float('inf')
        tight = p_perp / bound if bound > 0 and bound < float('inf') else 0

        print(f"  {delta:10.0e}  {cl_size:5d}  {p_perp:12.4e}  {r:12.4e}  {g:12.4e}  {bound:12.4e}  {tight:10.4f}")

        results.append({
            "delta": delta, "cl_size": cl_size, "p_perp": p_perp,
            "r": r, "g": g, "bound": bound, "tight": tight,
        })

    # Eigenspektrum-Ueberblick (erste 15)
    print(f"\n  Eigenwert-Spektrum (erste 15):")
    for k in range(min(15, dim)):
        print(f"    k={k:3d}: w={float(w_sorted[k]):12.6f}  gap={float(w_sorted[k]-w_min):12.4e}  |c|={abs(float(overlaps[k])):12.4e}")

    return results


def main():
    mp.mp.dps = DPS
    print("=" * 80)
    print("QUASIMODE-SKALIERUNG")
    print("=" * 80)

    all_results = {}
    for cfg in CONFIGS:
        res = run_quasimode(cfg["lam"], cfg["N"])
        all_results[cfg["lam"]] = res

    # Vergleich bei delta=0.01 (effektiver Cluster)
    print(f"\n{'='*80}")
    print(f"  SKALIERUNG BEI delta=0.01 (effektiver Cluster)")
    print(f"{'='*80}")
    print(f"  {'lambda':>8s}  {'#Cl':>5s}  {'||P_perp||':>12s}  {'Residual r':>12s}  {'Gap g':>12s}  {'r/g':>12s}  {'Tight':>8s}")
    for lam, res_list in all_results.items():
        for r in res_list:
            if r["delta"] == 1e-2:
                print(f"  {lam:8.1f}  {r['cl_size']:5d}  {r['p_perp']:12.4e}  {r['r']:12.4e}  {r['g']:12.4e}  {r['bound']:12.4e}  {r['tight']:8.4f}")


if __name__ == "__main__":
    main()
