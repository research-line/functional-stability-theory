"""
c2d_beta_verification.py — Verifiziere Proposition C2d.2 (beta_k-Formel)

Formel: beta_k = -<u_tilde, H v_k> / (C_tilde + <u_tilde, H u_tilde> - w_k)
wobei v_k = P_0 q_k = q_k - beta_k u_tilde

Teste:
1. Stimmt die Formel fuer alle k?
2. Pol-Decoupling D1: sum_{k not in Cl} |beta_k|^2 vs L^{-4}
3. Lambda-Skalierung von beta_k

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
from c2_approximation_test import inner_product, norm

DPS = int(os.environ.get("DPS", 50))

CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
]


def analyze(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L = float(2 * mp.log(lam_mp))

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

    # <u_tilde, H u_tilde>
    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))

    print(f"  QW+H in {time.time()-t0:.0f}s | L={L:.3f} C~={Ct:.3f} <u,Hu>={uHu:.4f}")

    cl_end = 0
    for k in range(dim):
        if float(ws[k] - wmin) < 1e-10:
            cl_end = k

    print(f"  Cluster: k=0..{cl_end} ({cl_end+1}/{dim})")
    print(f"\n  === C2d.2 VERIFIKATION: beta_k Formel ===")
    print(f"  {'k':>4s}  {'beta_num':>12s}  {'beta_C2d2':>12s}  {'rel_err':>10s}  "
          f"{'Nenner':>10s}  {'gap':>10s}")

    sum_beta2_nc = 0
    sum_beta2_cl = 0

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

        # <u_tilde, H v_k>
        Hvk = mp.matrix(dim, 1)
        for i in range(dim):
            Hvk[i, 0] = sum(H[i, j] * vk[j, 0] for j in range(dim))
        uHv = float(inner_product(ut, Hvk, dim))

        denom = Ct + uHu - wk
        beta_formula = -uHv / denom if abs(denom) > 1e-50 else float('inf')

        rel_err = abs(beta_formula - beta_num) / max(abs(beta_num), 1e-30)

        if k > cl_end:
            sum_beta2_nc += beta_num ** 2
        else:
            sum_beta2_cl += beta_num ** 2

        if k <= cl_end + 8 or k >= dim - 3:
            print(f"  {k:4d}  {beta_num:12.6e}  {beta_formula:12.6e}  {rel_err:10.2e}  "
                  f"{denom:10.3f}  {gap:10.2e}")

    print(f"\n  === D1: POL-DECOUPLING ===")
    print(f"  sum |beta_k|^2 (Cluster):      {sum_beta2_cl:.6f}")
    print(f"  sum |beta_k|^2 (Non-Cluster):   {sum_beta2_nc:.6f}")
    print(f"  L^4 = {L**4:.3f}")
    print(f"  sum_nc / L^{-4} = {sum_beta2_nc * L**4:.3f}")
    print(f"  ||P_perp u_tilde||^2 = {1 - sum_beta2_cl:.6f}")

    # u_tilde cluster content
    u_cl = sum_beta2_cl
    u_nc = sum_beta2_nc
    print(f"\n  u_tilde Lokalisierung:")
    print(f"    Cluster:      {u_cl*100:.1f}%")
    print(f"    Non-Cluster:  {u_nc*100:.1f}%")
    print(f"    Check:        {u_cl+u_nc:.6f} (soll = 1.000)")

    return dict(lam=lam, L=L, Ct=Ct, uHu=uHu, cl=cl_end+1,
                beta2_nc=sum_beta2_nc, beta2_cl=sum_beta2_cl)


def main():
    results = []
    for cfg in CONFIGS:
        lam, N = cfg["lam"], cfg["N"]
        print(f"\n{'='*80}")
        print(f"  LAMBDA = {lam}, N = {N}")
        print(f"{'='*80}")
        r = analyze(lam, N)
        results.append(r)

    if len(results) >= 2:
        r0, r1 = results[0], results[1]
        ll = np.log(r1["lam"] / r0["lam"])
        print(f"\n{'='*80}")
        print(f"  LAMBDA-SKALIERUNG")
        print(f"{'='*80}")

        rows = [
            ("sum beta_nc^2", r0["beta2_nc"], r1["beta2_nc"]),
            ("beta_nc^2 * L^4", r0["beta2_nc"] * r0["L"]**4,
             r1["beta2_nc"] * r1["L"]**4),
        ]
        print(f"  {'Metrik':<22s}  {'lam='+str(r0['lam']):>10s}  {'lam='+str(r1['lam']):>10s}  "
              f"{'Exp':>8s}")
        for name, v0, v1 in rows:
            exp = np.log(v1 / v0) / ll if v0 > 0 and v1 > 0 else float('nan')
            print(f"  {name:<22s}  {v0:10.4e}  {v1:10.4e}  {exp:+8.2f}")


if __name__ == "__main__":
    main()
