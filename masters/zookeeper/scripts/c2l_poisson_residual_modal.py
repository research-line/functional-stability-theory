"""
c2l_poisson_residual_modal.py — Modale Zerlegung des Poisson-Residuals

Berechnet r_j = <e_j, (H - w_min I) k_full> fuer alle T-Eigenmoden.

Aequivalent (C2j.1): r_j = (t_j - w_min)*g_j + alpha*f_j

Ziel: Messe den "Residual-Dominanz-Ratio" kappa_j = |r_j| / (alpha * |f_j|)
auf J_reg. C2i.2 sagt: kappa_j < 1 => phi_j < 0 (Vorzeichen erzwungen).

Zerlegung:
  - g_j^cl = -alpha_cl * f_j / (t_j - w_min)  [Cluster-Anteil, C2l.2]
  - g_j^perp = g_j - g_j^cl                    [Off-Cluster-Anteil]
  - r_j = (t_j - w_min) * g_j^perp             [Residual = Off-Cluster * Abstand]

Lambda-Skalierung: Wie skaliert kappa_j, ||g_perp||, ||r_off|| mit lambda?

Autor: LG (Opus 4.6, Session 27)
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

    print(f"\n{'='*72}")
    print(f"  POISSON-RESIDUAL MODAL:  lambda={lam}, N={N}, dim={dim}")
    print(f"{'='*72}")

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
    C_tilde = float(cn * cn)

    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))
    print(f"  Build: {time.time()-t0:.0f}s", flush=True)

    # k_full
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
    alpha = float(inner_product(ut, kn, dim))

    # f = P_0 H ut, g = P_0 k
    f_lam = mp.matrix(dim, 1)
    g_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]
        g_lam[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]

    # T = P_0 H P_0
    T = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            T[i, j] = (H[i, j]
                        - ut[i, 0] * Hut[j, 0]
                        - Hut[i, 0] * ut[j, 0]
                        + mp.mpf(uHu) * ut[i, 0] * ut[j, 0])

    # Diagonalize T and A
    t0 = time.time()
    ts, Es = diagonalize_mp(T, verbose=False)
    ws, Qs = diagonalize_mp(Aq, verbose=False)
    print(f"  Diag T+A: {time.time()-t0:.0f}s", flush=True)

    w_min = float(ws[0])

    # Cluster identification (A-eigenvalues)
    cl_A = [k for k in range(dim) if abs(float(ws[k]) - w_min) < 1e-10]
    cl_size = len(cl_A)

    # Cluster projection of k_full: alpha_cl
    c_coeffs = []
    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]
        c_coeffs.append(float(inner_product(qk, kn, dim)))

    alpha_cl = 0.0
    for k in cl_A:
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]
        u_comp = float(inner_product(ut, qk, dim))
        alpha_cl += c_coeffs[k] * u_comp

    cl_sq = sum(c_coeffs[k]**2 for k in cl_A)
    off_sq = sum(c_coeffs[k]**2 for k in range(dim) if k not in cl_A)

    print(f"\n  alpha = {alpha:.6f}, alpha_cl = {alpha_cl:.6f}, ratio = {alpha_cl/alpha:.6f}")
    print(f"  ||P_Cl k||^2 = {cl_sq:.8f}, ||P_perp k||^2 = {off_sq:.2e}")
    print(f"  Cluster A: {cl_size}/{dim}, w_min = {w_min:.6f}")

    # T-eigenvalues and components
    t_vals = []
    t_evecs = []
    fj_arr = []
    gj_arr = []
    t_is_zero = []

    for j in range(dim):
        tj = float(ts[j])
        ev = mp.matrix(dim, 1)
        for i in range(dim):
            ev[i, 0] = Es[i, j]

        is_zero = abs(tj) < 1e-12
        t_is_zero.append(is_zero)
        t_vals.append(tj)
        t_evecs.append(ev)

        if not is_zero:
            fj_arr.append(float(inner_product(ev, f_lam, dim)))
            gj_arr.append(float(inner_product(ev, g_lam, dim)))
        else:
            fj_arr.append(0.0)
            gj_arr.append(0.0)

    # === MODAL RESIDUAL ===
    print(f"\n  === MODALE RESIDUAL-ZERLEGUNG ===")
    print(f"  r_j = (t_j - w_min)*g_j + alpha*f_j")
    print(f"  kappa_j = |r_j| / (alpha * |f_j|)   [< 1 => Vorzeichen erzwungen]")
    print(f"  g_j^perp = r_j / (t_j - w_min)       [Off-Cluster-Anteil]")

    print(f"\n  {'j':>4}  {'t_j':>10}  {'f_j':>10}  {'g_j':>10}  {'r_j':>12}  "
          f"{'kappa_j':>10}  {'|g_perp|':>10}  {'sign_ok':>8}")
    print(f"  {'-'*88}")

    kappa_vals = []
    r_vals = []
    g_perp_vals = []
    n_sign_ok = 0
    n_reg = 0

    for j in range(dim):
        if t_is_zero[j]:
            continue
        tj = t_vals[j]
        fj = fj_arr[j]
        gj = gj_arr[j]

        denom = tj - w_min
        if abs(denom) < 1e-12:
            continue

        rj = denom * gj + alpha * fj
        g_perp_j = rj / denom

        if abs(fj) < 1e-25:
            continue

        kappa_j = abs(rj) / (abs(alpha) * abs(fj))
        sign_ok = (fj * gj < 0)

        n_reg += 1
        if sign_ok:
            n_sign_ok += 1

        kappa_vals.append(kappa_j)
        r_vals.append(abs(rj))
        g_perp_vals.append(abs(g_perp_j))

        if j < 10 or j > dim - 5 or kappa_j > 0.5:
            ok_str = "OK" if sign_ok else "FAIL"
            print(f"  {j:4d}  {tj:10.4f}  {fj:10.4e}  {gj:10.4e}  {rj:12.4e}  "
                  f"{kappa_j:10.6f}  {abs(g_perp_j):10.4e}  {ok_str:>8}")
        elif j == 10:
            print(f"  {'...':>4}")

    kappa_arr = np.array(kappa_vals)
    r_arr = np.array(r_vals)
    gp_arr = np.array(g_perp_vals)

    print(f"\n  === STATISTIK (J_reg: {n_reg} Moden) ===")
    print(f"  Sign OK: {n_sign_ok}/{n_reg}")
    print(f"  kappa_j < 1 (=> sign forced): {np.sum(kappa_arr < 1)}/{len(kappa_arr)}")
    print(f"  kappa_j Statistik:")
    print(f"    Max    = {np.max(kappa_arr):.6e}")
    print(f"    Median = {np.median(kappa_arr):.6e}")
    print(f"    Mean   = {np.mean(kappa_arr):.6e}")
    print(f"    Min    = {np.min(kappa_arr):.6e}")

    print(f"\n  |r_j| Statistik:")
    print(f"    Max    = {np.max(r_arr):.6e}")
    print(f"    RMS    = {np.sqrt(np.mean(r_arr**2)):.6e}")
    print(f"    ||r||  = {np.sqrt(np.sum(r_arr**2)):.6e}")

    print(f"\n  |g_j^perp| Statistik:")
    print(f"    Max    = {np.max(gp_arr):.6e}")
    print(f"    RMS    = {np.sqrt(np.mean(gp_arr**2)):.6e}")
    print(f"    ||g_perp|| = {np.sqrt(np.sum(gp_arr**2)):.6e}")

    # === BOUNDARY ANALYSIS ===
    print(f"\n  === BOUNDARY-ANALYSE ===")
    print(f"  (H - w_min) k_full: volle Residualnorm")

    Hk = mp.matrix(dim, 1)
    for i in range(dim):
        Hk[i, 0] = sum(H[i, jj] * kn[jj, 0] for jj in range(dim))

    res_full = mp.matrix(dim, 1)
    for i in range(dim):
        res_full[i, 0] = Hk[i, 0] - mp.mpf(w_min) * kn[i, 0]

    res_norm = float(norm(res_full, dim))
    res_u = float(inner_product(ut, res_full, dim))

    Ak = mp.matrix(dim, 1)
    for i in range(dim):
        Ak[i, 0] = Hk[i, 0] + mp.mpf(C_tilde) * mp.mpf(alpha) * ut[i, 0]
    Ak_res = mp.matrix(dim, 1)
    for i in range(dim):
        Ak_res[i, 0] = Ak[i, 0] - mp.mpf(w_min) * kn[i, 0]
    Ak_res_norm = float(norm(Ak_res, dim))

    print(f"  ||(H - w_min)k||   = {res_norm:.6e}")
    print(f"  ||(A - w_min)k||   = {Ak_res_norm:.6e}")
    print(f"  <u, (H-w_min)k>    = {res_u:.6e}")
    print(f"  -C_tilde * alpha   = {-C_tilde * alpha:.6e}  (expected u-component)")
    print(f"  Ratio              = {res_u / (-C_tilde * alpha):.6f}")

    return {
        "lam": lam,
        "kappa_max": float(np.max(kappa_arr)),
        "kappa_median": float(np.median(kappa_arr)),
        "kappa_mean": float(np.mean(kappa_arr)),
        "r_rms": float(np.sqrt(np.mean(r_arr**2))),
        "r_norm": float(np.sqrt(np.sum(r_arr**2))),
        "gp_norm": float(np.sqrt(np.sum(gp_arr**2))),
        "off_cl_sq": off_sq,
        "Ak_res_norm": Ak_res_norm,
        "n_sign_ok": n_sign_ok,
        "n_reg": n_reg,
    }


if __name__ == "__main__":
    results = []
    for cfg in CONFIGS:
        r = run(cfg["lam"], cfg["N"])
        results.append(r)

    if len(results) >= 2:
        print(f"\n{'='*72}")
        print(f"  LAMBDA-SKALIERUNG")
        print(f"{'='*72}")
        print(f"  {'Metrik':>20}  ", end="")
        for r in results:
            print(f"{'lam='+str(r['lam']):>12}  ", end="")
        if len(results) >= 2:
            print(f"{'Exponent':>10}")
        else:
            print()

        metrics = [
            ("kappa_max", "kappa_max"),
            ("kappa_median", "kappa_median"),
            ("||r||", "r_norm"),
            ("||g_perp||", "gp_norm"),
            ("||P_perp k||^2", "off_cl_sq"),
            ("||(A-w*)k||", "Ak_res_norm"),
        ]

        import math
        for name, key in metrics:
            vals = [r[key] for r in results]
            print(f"  {name:>20}  ", end="")
            for v in vals:
                print(f"{v:12.4e}  ", end="")
            if len(vals) >= 2 and vals[0] > 0 and vals[-1] > 0:
                lams = [r["lam"] for r in results]
                exp = math.log(vals[-1]/vals[0]) / math.log(lams[-1]/lams[0])
                print(f"{exp:10.2f}")
            else:
                print()

        print(f"\n  Sign-Match:")
        for r in results:
            print(f"    lambda={r['lam']}: {r['n_sign_ok']}/{r['n_reg']}")

    print(f"\n=== FERTIG ===")
