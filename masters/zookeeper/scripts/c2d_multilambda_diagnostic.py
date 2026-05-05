"""
c2d_multilambda_diagnostic.py — Multi-Lambda Diagnostik fuer D1/D2

Berechnet fuer lambda = 3, 5, 7:
1. alpha = <u_tilde, k_full> (Konvergenz-Test)
2. ||u_perp||^2 = 1 - sum_cl |beta_k|^2 (D1 = u-Lokalisierung)
3. Spektral-Gap: w_{cl+1} - w_{cl} (Weg B fuer D2)
4. Zaehler/Nenner-Split der beta_k-Formel
5. Per-EV: |N_k|, |D_k|, |N_k/D_k| Skalierung

Autor: LG (Opus 4.6, Session 24)
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
    k_lambda_value, h_educated_guess, inner_product, norm,
)
from c2_poisson_decomposition import project_to_fourier

DPS = int(os.environ.get("DPS", 50))
CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
    {"lam": 7.0, "N": 85},
]


def analyze(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

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

    ws, Qs = diagonalize_mp(Aq, verbose=False)
    wmin = ws[0]
    Ct = float(sum(W02[i, i] for i in range(dim)))

    col0 = mp.matrix(dim, 1)
    for i in range(dim):
        col0[i, 0] = W02[i, 0]
    cn = norm(col0, dim)
    ut = mp.matrix(dim, 1)
    for i in range(dim):
        ut[i, 0] = col0[i, 0] / cn
    t_qw = time.time() - t0

    # <u_tilde, H u_tilde>
    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))

    # k_full Fourier-Projektion
    t0 = time.time()

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
    t_proj = time.time() - t0

    alpha = float(inner_product(ut, kn, dim))

    # Cluster-Bestimmung
    cl_end = 0
    for k in range(dim):
        if float(ws[k] - wmin) < 1e-10:
            cl_end = k

    # Spektral-Gap
    if cl_end + 1 < dim:
        gap_cluster = float(ws[cl_end + 1] - ws[cl_end])
    else:
        gap_cluster = 0.0

    print(f"  QW={t_qw:.0f}s Proj={t_proj:.0f}s")
    print(f"  L={L:.4f} C~={Ct:.4f} <u,Hu>={uHu:.5f}")
    print(f"  Cluster: k=0..{cl_end} ({cl_end+1}/{dim})")
    print(f"  Spektral-Gap (cl_end -> cl_end+1): {gap_cluster:.4e}")
    print(f"  alpha = <u_tilde, k_full> = {alpha:.8f}")
    print(f"  alpha^2 = {alpha**2:.8f}")

    # Per-EV Analyse
    sum_beta2_cl = 0.0
    sum_beta2_nc = 0.0
    sum_num2_nc = 0.0
    sum_denom2_inv_nc = 0.0
    first_nc_data = []

    print(f"\n  {'k':>4s}  {'beta_k':>12s}  {'|N_k|':>12s}  {'D_k':>12s}  "
          f"{'|N_k|/D_k':>12s}  {'gap':>10s}  {'c_k':>12s}")

    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]

        beta_num = float(inner_product(ut, qk, dim))
        wk = float(ws[k])
        gap = wk - float(wmin)

        # v_k = q_k - beta_k u_tilde
        vk = mp.matrix(dim, 1)
        for i in range(dim):
            vk[i, 0] = qk[i, 0] - mp.mpf(beta_num) * ut[i, 0]

        # Zaehler: <u_tilde, H v_k>
        Hvk = mp.matrix(dim, 1)
        for i in range(dim):
            Hvk[i, 0] = sum(H[i, j] * vk[j, 0] for j in range(dim))
        uHv = float(inner_product(ut, Hvk, dim))
        num_k = abs(uHv)

        # Nenner
        denom = Ct + uHu - wk

        # c_k
        ck = float(inner_product(kn, qk, dim))

        if k > cl_end:
            sum_beta2_nc += beta_num ** 2
            sum_num2_nc += uHv ** 2
            sum_denom2_inv_nc += 1.0 / denom ** 2
            if len(first_nc_data) < 6:
                first_nc_data.append(dict(
                    k=k, beta=beta_num, num=uHv, denom=denom,
                    gap=gap, ck=ck
                ))
        else:
            sum_beta2_cl += beta_num ** 2

        if k <= cl_end + 1 or (k > cl_end and k < cl_end + 6):
            marker = " <-- CL-RAND" if k == cl_end else ""
            print(f"  {k:4d}  {beta_num:12.4e}  {num_k:12.4e}  {denom:12.4f}  "
                  f"{num_k/abs(denom):12.4e}  {gap:10.2e}  {ck:12.4e}{marker}")

    u_perp2 = 1.0 - sum_beta2_cl
    mean_num2 = sum_num2_nc / max(dim - cl_end - 1, 1)

    print(f"\n  === AGGREGATE ===")
    print(f"  ||u_perp||^2 = {u_perp2:.8f}")
    print(f"  sum |beta|^2 (nc) = {sum_beta2_nc:.8f}")
    print(f"  sum |N_k|^2 (nc) = {sum_num2_nc:.6e}")
    print(f"  mean |N_k|^2 (nc) = {mean_num2:.6e}")

    return dict(
        lam=lam, L=L, Ct=Ct, uHu=uHu, alpha=alpha,
        cl=cl_end + 1, dim=dim, gap=gap_cluster,
        u_perp2=u_perp2, beta2_nc=sum_beta2_nc,
        num2_nc=sum_num2_nc, mean_num2=mean_num2,
        first_nc=first_nc_data
    )


def main():
    results = []
    for cfg in CONFIGS:
        lam, N = cfg["lam"], cfg["N"]
        print(f"\n{'='*80}")
        print(f"  LAMBDA = {lam}, N = {N}")
        print(f"{'='*80}")
        r = analyze(lam, N)
        results.append(r)

    if len(results) < 2:
        return

    print(f"\n{'='*80}")
    print(f"  MULTI-LAMBDA SKALIERUNG")
    print(f"{'='*80}")

    def exp_fit(vals, lams):
        if len(vals) < 2:
            return float("nan")
        log_v = [np.log(abs(v)) for v in vals if abs(v) > 0]
        log_l = [np.log(l) for l, v in zip(lams, vals) if abs(v) > 0]
        if len(log_v) < 2:
            return float("nan")
        # Least-squares linear fit
        n = len(log_v)
        sx = sum(log_l)
        sy = sum(log_v)
        sxx = sum(x ** 2 for x in log_l)
        sxy = sum(x * y for x, y in zip(log_l, log_v))
        return (n * sxy - sx * sy) / (n * sxx - sx ** 2)

    lams = [r["lam"] for r in results]

    rows = [
        ("alpha", [r["alpha"] for r in results]),
        ("alpha^2", [r["alpha"] ** 2 for r in results]),
        ("||u_perp||^2", [r["u_perp2"] for r in results]),
        ("sum |beta|^2 (nc)", [r["beta2_nc"] for r in results]),
        ("sum |N_k|^2 (nc)", [r["num2_nc"] for r in results]),
        ("mean |N_k|^2 (nc)", [r["mean_num2"] for r in results]),
        ("Spektral-Gap", [r["gap"] for r in results]),
        ("L^2 = C~/ratio", [r["L"] ** 2 for r in results]),
        ("L^{-4}", [r["L"] ** (-4) for r in results]),
    ]

    header = f"  {'Metrik':<22s}"
    for r in results:
        header += f"  {'lam='+str(r['lam']):>12s}"
    header += f"  {'Exp(fit)':>10s}"
    print(f"\n{header}")

    for name, vals in rows:
        line = f"  {name:<22s}"
        for v in vals:
            line += f"  {v:12.6e}"
        exp = exp_fit(vals, lams)
        line += f"  {exp:+10.4f}"
        print(line)

    # Alpha-Konvergenz
    print(f"\n  === ALPHA-KONVERGENZ ===")
    for r in results:
        print(f"  lam={r['lam']:5.1f}: alpha = {r['alpha']:.8f}")

    # D1-Zerlegung: Zaehler vs Nenner
    print(f"\n  === D1-ZERLEGUNG: ZAEHLER vs NENNER ===")
    print(f"  D1-Rate = Zaehler-Rate + Nenner-Rate")
    num_vals = [r["num2_nc"] for r in results]
    beta_vals = [r["beta2_nc"] for r in results]
    print(f"  sum|beta|^2 Exp:   {exp_fit(beta_vals, lams):+.4f}")
    print(f"  sum|N_k|^2 Exp:    {exp_fit(num_vals, lams):+.4f}")
    print(f"  L^{{-4}} Exp:        {exp_fit([r['L']**(-4) for r in results], lams):+.4f}")
    print(f"  Extra (ueber L^-4): {exp_fit(beta_vals, lams) - exp_fit([r['L']**(-4) for r in results], lams):+.4f}")

    # Per-EV Vergleich (1. non-cluster)
    print(f"\n  === ERSTER NICHT-CLUSTER-EV ===")
    nc_data = [r["first_nc"][0] for r in results if r["first_nc"]]
    if len(nc_data) >= 2:
        nc_betas = [abs(d["beta"]) for d in nc_data]
        nc_nums = [abs(d["num"]) for d in nc_data]
        nc_denoms = [abs(d["denom"]) for d in nc_data]
        nc_cks = [abs(d["ck"]) for d in nc_data]
        nc_gaps = [d["gap"] for d in nc_data]

        print(f"  |beta| Exp:  {exp_fit(nc_betas, lams):+.4f}")
        print(f"  |N_k| Exp:   {exp_fit(nc_nums, lams):+.4f}")
        print(f"  |D_k| Exp:   {exp_fit(nc_denoms, lams):+.4f}")
        print(f"  |c_k| Exp:   {exp_fit(nc_cks, lams):+.4f}")
        print(f"  gap Exp:     {exp_fit(nc_gaps, lams):+.4f}")
        print(f"\n  Per-Lambda:")
        for r, d in zip(results, nc_data):
            print(f"    lam={r['lam']:5.1f}: k={d['k']:3d}, |beta|={abs(d['beta']):.4f}, "
                  f"|N|={abs(d['num']):.4e}, D={abs(d['denom']):.4f}, "
                  f"|c|={abs(d['ck']):.4e}, gap={d['gap']:.4e}")


if __name__ == "__main__":
    main()
