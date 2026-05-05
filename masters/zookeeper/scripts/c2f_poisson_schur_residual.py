"""
c2f_poisson_schur_residual.py — Schur-Komplement-Defekte der Poisson-Konstruktion

Testet die NICHT-ZIRKULAERE Bruecke: Erfuellt k_full die Schur-Komplement-Gleichungen?

k_full = alpha * u_tilde + g,  g = P_0 k_full

Schur-Defekte:
  r_perp := (T - w* I) g + alpha f       [P_0-Raum-Residuum]
  r_sc   := (C_tilde + uHu - w*) alpha + <f, g>   [Secular-Residuum]

Falls beide klein: g approx -alpha (T - w* I)^{-1} f OHNE MS2 vorauszusetzen.

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


def run(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    print(f"\n{'='*65}")
    print(f"  lambda = {lam}, N = {N}, dim = {dim}")
    print(f"{'='*65}")

    # === Build A = C_tilde |ut><ut| + H ===
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

    # u_tilde (normalized first column of W02)
    col0 = mp.matrix(dim, 1)
    for i in range(dim):
        col0[i, 0] = W02[i, 0]
    cn = norm(col0, dim)
    ut = mp.matrix(dim, 1)
    for i in range(dim):
        ut[i, 0] = col0[i, 0] / cn

    Ct = float(sum(W02[i, i] for i in range(dim)))
    print(f"  Build: {time.time()-t0:.0f}s, C_tilde = {Ct:.4f}", flush=True)

    # === Diagonalize A to get w_min ===
    t0 = time.time()
    ws, Qs = diagonalize_mp(Aq, verbose=False)
    wmin = float(ws[0])
    print(f"  Diag A: {time.time()-t0:.0f}s, w_min = {wmin:.6f}", flush=True)

    # === H u_tilde ===
    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))

    # === f = P_0 H ut = Hut - uHu * ut ===
    f_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]

    f_norm = float(norm(f_lam, dim))

    # === k_full (Poisson construction) ===
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

    # === Decompose k_full = alpha * ut + g ===
    alpha = float(inner_product(ut, kn, dim))

    g = mp.matrix(dim, 1)
    for i in range(dim):
        g[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]

    g_norm = float(norm(g, dim))

    # === T = P_0 H P_0 ===
    T = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            T[i, j] = H[i, j] - ut[i, 0] * Hut[j, 0] - Hut[i, 0] * ut[j, 0] + mp.mpf(uHu) * ut[i, 0] * ut[j, 0]

    # === Schur Residual 1: r_perp = (T - w* I) g + alpha f ===
    t0 = time.time()
    Tg = mp.matrix(dim, 1)
    for i in range(dim):
        Tg[i, 0] = sum(T[i, j] * g[j, 0] for j in range(dim))

    wmin_mp = mp.mpf(wmin)
    r_perp = mp.matrix(dim, 1)
    for i in range(dim):
        r_perp[i, 0] = Tg[i, 0] - wmin_mp * g[i, 0] + mp.mpf(alpha) * f_lam[i, 0]

    r_perp_norm = float(norm(r_perp, dim))
    t_res = time.time() - t0

    # === Schur Residual 2: r_sc = (C_tilde + uHu - w*) alpha + <f, g> ===
    fg = float(inner_product(f_lam, g, dim))
    r_sc = (Ct + uHu - wmin) * alpha + fg

    # === Full residual: (A - w* I) k_full ===
    Akn = mp.matrix(dim, 1)
    for i in range(dim):
        Akn[i, 0] = sum(Aq[i, j] * kn[j, 0] for j in range(dim))

    full_res = mp.matrix(dim, 1)
    for i in range(dim):
        full_res[i, 0] = Akn[i, 0] - wmin_mp * kn[i, 0]

    full_res_norm = float(norm(full_res, dim))

    # === Rayleigh quotient ===
    kAk = float(inner_product(kn, Akn, dim))
    rayleigh = kAk  # since ||k|| = 1

    # === Report ===
    print(f"\n  === SCHUR-DEFEKTE ===")
    print(f"  alpha       = {alpha:.6f}")
    print(f"  ||g||       = {g_norm:.6f}")
    print(f"  ||f||       = {f_norm:.6f}")
    print(f"  C_tilde     = {Ct:.4f}")
    print(f"  <u,Hu>      = {uHu:.6f}")
    print(f"  w*=w_min    = {wmin:.6f}")
    print(f"  <f, g>      = {fg:.6e}")

    print(f"\n  --- Residuen ---")
    print(f"  ||r_perp||  = {r_perp_norm:.6e}   (Nullraum-Defekt)")
    print(f"  |r_sc|      = {abs(r_sc):.6e}   (Secular-Defekt)")
    print(f"  ||(A-w*I)k||= {full_res_norm:.6e}   (Voller Eigenvektor-Defekt)")
    print(f"  Rayleigh    = {rayleigh:.10f}   (vs w_min = {wmin:.10f})")
    print(f"  |Ray - w_min| = {abs(rayleigh - wmin):.6e}")

    print(f"\n  --- Relative Residuen ---")
    print(f"  ||r_perp||/||g||     = {r_perp_norm/g_norm:.6e}")
    print(f"  ||r_perp||/(|a|·||f||) = {r_perp_norm/(abs(alpha)*f_norm):.6e}")
    print(f"  |r_sc|/(|a|·|C+uHu-w|) = {abs(r_sc)/(abs(alpha)*abs(Ct+uHu-wmin)):.6e}")

    # === Decompose r_perp in A-eigenbasis ===
    print(f"\n  --- r_perp in A-Eigenbasis (Top-5) ---")
    overlaps = []
    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]
        ov = float(inner_product(qk, r_perp, dim))
        overlaps.append((k, float(ws[k]), ov))

    overlaps.sort(key=lambda x: -abs(x[2]))
    for k, wk, ov in overlaps[:8]:
        wdiff = wk - wmin
        cl = "[CL]" if abs(wdiff) < 1e-10 else ""
        print(f"    k={k:3d}  w={wk:10.4f}  <q_k,r_perp>={ov:12.4e}  gap={wdiff:.4e}  {cl}")

    return {
        "lam": lam, "N": N,
        "r_perp": r_perp_norm,
        "r_sc": abs(r_sc),
        "full_res": full_res_norm,
        "rayleigh_err": abs(rayleigh - wmin),
        "alpha": alpha,
        "g_norm": g_norm,
        "f_norm": f_norm,
        "r_perp_rel": r_perp_norm / g_norm,
    }


if __name__ == "__main__":
    results = []
    for cfg in CONFIGS:
        res = run(cfg["lam"], cfg["N"])
        results.append(res)

    # === Scaling summary ===
    print(f"\n{'='*65}")
    print(f"  SKALIERUNG DER SCHUR-DEFEKTE")
    print(f"{'='*65}")
    print(f"  {'lam':>5}  {'||r_perp||':>12}  {'|r_sc|':>12}  {'||(A-w)k||':>12}  {'||r_p||/||g||':>14}  {'|Ray-w|':>12}")
    for r in results:
        print(f"  {r['lam']:5.1f}  {r['r_perp']:12.4e}  {r['r_sc']:12.4e}  {r['full_res']:12.4e}  {r['r_perp_rel']:14.4e}  {r['rayleigh_err']:12.4e}")

    if len(results) >= 2:
        import math
        for name, key in [("||r_perp||", "r_perp"), ("|r_sc|", "r_sc"),
                          ("||(A-w)k||", "full_res"), ("||r_p||/||g||", "r_perp_rel")]:
            vals = [(r["lam"], r[key]) for r in results if r[key] > 0]
            if len(vals) >= 2:
                lams = [math.log(v[0]) for v in vals]
                logs = [math.log(v[1]) for v in vals]
                # Simple 2-point from first and last
                exp2 = (logs[-1] - logs[0]) / (lams[-1] - lams[0])
                print(f"  {name:20s}  Exponent (2pt) = {exp2:.2f}")

    print(f"\n=== FERTIG ===")
