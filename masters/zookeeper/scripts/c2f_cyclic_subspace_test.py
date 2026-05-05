"""
c2f_cyclic_subspace_test.py — Zyklischer Unterraum-Test

Testet: Liegt g = P_0 k_full im zyklischen Unterraum C_f = span{f, Tf, T^2f, ...}?

Zwei Methoden:
  1. Krylov-Distanz: dist(g, span{f, Tf, ..., T^m f}) fuer m = 0,1,...
  2. Spektraler Support: Gibt es e_j mit <e_j,f> ≈ 0 aber <e_j,g> ≠ 0?

Falls g in C_f, existiert phi mit g = phi(T)f, und phi(t_j) = <e_j,g>/<e_j,f>.
Vorzeichenkohaerenz ist dann aequivalent zu phi(t_j) <= 0.

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


def mat_vec(M, v, dim):
    r = mp.matrix(dim, 1)
    for i in range(dim):
        r[i, 0] = sum(M[i, j] * v[j, 0] for j in range(dim))
    return r


def run(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    print(f"\n{'='*65}")
    print(f"  lambda = {lam}, N = {N}, dim = {dim}")
    print(f"{'='*65}")

    # Build operators
    t0 = time.time()
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

    print(f"  Build: {time.time()-t0:.0f}s", flush=True)

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

    # f = P_0 H ut,  g = P_0 k
    f_lam = mp.matrix(dim, 1)
    g_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]
        g_lam[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]

    g_norm = float(norm(g_lam, dim))
    f_norm = float(norm(f_lam, dim))

    # T = P_0 H P_0
    T = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            T[i, j] = H[i, j] - ut[i, 0] * Hut[j, 0] - Hut[i, 0] * ut[j, 0] + mp.mpf(uHu) * ut[i, 0] * ut[j, 0]

    # === Part 1: Krylov distance ===
    print(f"\n  === KRYLOV-DISTANZ dist(g, span{{f, Tf, ..., T^m f}}) ===")
    print(f"  {'m':>3}  {'dist':>14}  {'dist/||g||':>14}  {'log10(dist/g)':>14}")

    import math

    # Lanczos-style Krylov basis: q_0 = f/||f||, then q_{m+1} from T q_m
    basis_vecs = []  # orthonormal basis vectors

    max_m = min(dim - 1, 40)

    for m in range(max_m + 1):
        if m == 0:
            cur = mp.matrix(dim, 1)
            for i in range(dim):
                cur[i, 0] = f_lam[i, 0]
        else:
            cur = mat_vec(T, basis_vecs[-1], dim)

        # Gram-Schmidt against all existing basis vectors
        for bv in basis_vecs:
            ov = inner_product(bv, cur, dim)
            for i in range(dim):
                cur[i, 0] -= ov * bv[i, 0]

        cur_norm = float(norm(cur, dim))
        if cur_norm < 1e-40:
            print(f"  {m:3d}  Krylov space exhausted (norm = {cur_norm:.2e})")
            break

        normed = mp.matrix(dim, 1)
        for i in range(dim):
            normed[i, 0] = cur[i, 0] / mp.mpf(cur_norm)
        basis_vecs.append(normed)

        # Project g onto current Krylov space
        g_proj = mp.matrix(dim, 1)
        for bv in basis_vecs:
            ov = inner_product(bv, g_lam, dim)
            for i in range(dim):
                g_proj[i, 0] += ov * bv[i, 0]

        resid = mp.matrix(dim, 1)
        for i in range(dim):
            resid[i, 0] = g_lam[i, 0] - g_proj[i, 0]

        dist = float(norm(resid, dim))
        rel_dist = dist / g_norm

        log_dist = math.log10(rel_dist) if rel_dist > 1e-50 else -50
        print(f"  {m:3d}  {dist:14.6e}  {rel_dist:14.6e}  {log_dist:14.2f}")

        if rel_dist < 1e-40:
            print(f"  => g liegt EXAKT im Krylov-Unterraum (m={m})")
            break

    # === Part 2: Spectral support check ===
    print(f"\n  === SPEKTRALER SUPPORT-CHECK ===")
    t0 = time.time()
    ts, Es = diagonalize_mp(T, verbose=False)
    print(f"  Diag T: {time.time()-t0:.0f}s", flush=True)

    print(f"\n  {'j':>3}  {'t_j':>10}  {'<e,f>':>12}  {'<e,g>':>12}  {'phi_j':>12}  {'sign':>5}")

    n_f_zero = 0
    n_g_leak = 0
    phi_values = []

    for jj in range(dim):
        tj = float(ts[jj])
        ej = mp.matrix(dim, 1)
        for i in range(dim):
            ej[i, 0] = Es[i, jj]

        fj = float(inner_product(ej, f_lam, dim))
        gj = float(inner_product(ej, g_lam, dim))

        if abs(tj) < 1e-12:
            continue

        if abs(fj) < 1e-25:
            n_f_zero += 1
            if abs(gj) > 1e-20:
                n_g_leak += 1
                print(f"  {jj:3d}  {tj:10.4f}  {fj:12.4e}  {gj:12.4e}  {'LEAK!':>12}  {'!!!':>5}")
            continue

        phi_j = gj / fj
        sign = "+" if phi_j > 0 else "-"
        phi_values.append((jj, tj, phi_j, fj, gj))

        if abs(phi_j) > 0.01 or jj < 5 or jj > dim - 3:
            print(f"  {jj:3d}  {tj:10.4f}  {fj:12.4e}  {gj:12.4e}  {phi_j:12.6f}  {sign:>5}")

    print(f"\n  f-Nullmoden: {n_f_zero}")
    print(f"  g-Leckage (f=0, g≠0): {n_g_leak}")

    if phi_values:
        phis = [p[2] for p in phi_values]
        n_pos = sum(1 for p in phis if p > 0)
        n_neg = sum(1 for p in phis if p < 0)
        print(f"  phi > 0: {n_pos}")
        print(f"  phi < 0: {n_neg}")
        print(f"  phi = 0: {len(phis) - n_pos - n_neg}")

        if n_pos == 0:
            print(f"  => PERFEKTE NEGATIVITAET: phi(t_j) <= 0 fuer alle j")


if __name__ == "__main__":
    for cfg in CONFIGS:
        run(cfg["lam"], cfg["N"])
    print(f"\n=== FERTIG ===")
