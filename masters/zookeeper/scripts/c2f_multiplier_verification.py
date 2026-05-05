"""
c2f_multiplier_verification.py — Verifiziert den Resolvent-Multiplikator

Testet: r_j = -<e_j, g>/<e_j, f> = alpha/(t_j - w_min)

Wobei:
  f = P_0 H u_tilde,  g = P_0 k_full
  T = P_0 H P_0,  {e_j, t_j} = T-Eigenbasis
  alpha = <u_tilde, k_full>,  w_min = kleinster A-Eigenwert

Autor: LG (Opus 4.6, Session 26)
Datum: 2026-04-19
"""

import os, sys, time
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, diagonalize_mp, chi_trivial,
)
from c2_approximation_test import (
    k_lambda_value, inner_product, norm,
)
from c2_poisson_decomposition import project_to_fourier

DPS = int(os.environ.get("DPS", 50))
CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
    {"lam": 7.0, "N": 85},
]


def build_all(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    W02 = mp.matrix(dim, dim)
    H = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            W02[i, j] = Aq[i, j] - Ah[i, j]
            H[i, j] = Ah[i, j]

    col0 = mp.matrix(dim, 1)
    for i in range(dim):
        col0[i, 0] = W02[i, 0]
    cn = norm(col0, dim)
    ut = mp.matrix(dim, 1)
    for i in range(dim):
        ut[i, 0] = col0[i, 0] / cn

    return Aq, H, ut, dim, lam_mp, L_mp


def run(lam, N):
    print(f"\n{'='*60}")
    print(f"  lambda = {lam}, N = {N}")
    print(f"{'='*60}")

    t0 = time.time()
    Aq, H, ut, dim, lam_mp, L_mp = build_all(lam, N)
    print(f"  Build: {time.time()-t0:.0f}s", flush=True)

    # Diag A
    t0 = time.time()
    ws, Qs = diagonalize_mp(Aq, verbose=False)
    wmin = float(ws[0])
    print(f"  Diag A: {time.time()-t0:.0f}s, w_min = {wmin:.6f}", flush=True)

    # Cluster
    cl_end = 0
    for k in range(dim):
        if float(ws[k] - ws[0]) < 1e-10:
            cl_end = k

    # H u_tilde
    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))

    # k_full
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    t0 = time.time()
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
    print(f"  Proj k: {time.time()-t0:.0f}s", flush=True)

    alpha = float(inner_product(ut, kn, dim))
    print(f"  alpha = {alpha:.6f}", flush=True)

    # T = P_0 H P_0
    T = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            T[i, j] = H[i, j] - ut[i, 0] * Hut[j, 0] - Hut[i, 0] * ut[j, 0] + mp.mpf(uHu) * ut[i, 0] * ut[j, 0]

    t0 = time.time()
    ts, Es = diagonalize_mp(T, verbose=False)
    print(f"  Diag T: {time.time()-t0:.0f}s", flush=True)

    # f = P_0 H ut,  g = P_0 k
    f_lam = mp.matrix(dim, 1)
    g_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]
        g_lam[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]

    # Compute r_j = -<e_j,g>/<e_j,f> and model = alpha/(t_j - w_min)
    print(f"\n  Cluster-Ende: k={cl_end} (Cl = {cl_end+1}/{dim})")
    print(f"\n  {'j':>3}  {'t_j':>10}  {'r_j':>12}  {'model':>12}  {'ratio':>10}  {'|1-r|':>10}")
    print(f"  {'-'*65}")

    ratios = []
    for jj in range(dim):
        tj = float(ts[jj])
        ej = mp.matrix(dim, 1)
        for i in range(dim):
            ej[i, 0] = Es[i, jj]

        fj = float(inner_product(ej, f_lam, dim))
        gj = float(inner_product(ej, g_lam, dim))

        if abs(fj) < 1e-30 or abs(tj) < 1e-12:
            continue

        rj = -gj / fj
        denom = tj - wmin
        if abs(denom) < 1e-12:
            continue
        model = alpha / denom
        if abs(model) < 1e-30:
            continue

        ratio = rj / model
        err = abs(1 - ratio)

        is_cluster = (jj <= cl_end)
        marker = "  [CL]" if is_cluster else ""

        if not is_cluster or abs(rj) > 0.01:
            print(f"  {jj:3d}  {tj:10.4f}  {rj:12.4f}  {model:12.4f}  {ratio:10.7f}  {err:10.2e}{marker}")

        if not is_cluster and err < 0.1:
            ratios.append(err)

    if ratios:
        arr = np.array(ratios)
        print(f"\n  Off-Cluster Statistik ({len(arr)} Moden):")
        print(f"    Median |1-ratio| = {np.median(arr):.2e}")
        print(f"    Max    |1-ratio| = {np.max(arr):.2e}")
        print(f"    Mean   |1-ratio| = {np.mean(arr):.2e}")

    print(f"\n  => Multiplikator phi(t) = alpha/(t - w_min)")
    print(f"     alpha = {alpha:.6f}, w_min = {wmin:.6f}")


if __name__ == "__main__":
    for cfg in CONFIGS:
        run(cfg["lam"], cfg["N"])
    print("\n=== FERTIG ===")
