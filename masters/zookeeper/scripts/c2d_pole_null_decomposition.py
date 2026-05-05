"""
c2d_pole_null_decomposition.py — Zerlege c_k in Pol- und Nullraum-Beitrag

c_k = <q_k, k_full> = beta_k * alpha + <v_k, k_full>
                       ^--- Pol-Beitrag    ^--- Nullraum-Beitrag

wobei alpha = <u_tilde, k_full> und v_k = q_k - beta_k * u_tilde.

Frage: Kommt die Power-Law-Rate von c_k aus D1 (beta_k -> 0)
oder aus D2 (Pol-Null-Cancellation wird besser)?

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
    for i in range(dim):
        for j in range(dim):
            W02[i, j] = Aq[i, j] - Ah[i, j]

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

    cl_end = 0
    for k in range(dim):
        if float(ws[k] - wmin) < 1e-10:
            cl_end = k

    print(f"  QW={t_qw:.0f}s Proj={t_proj:.0f}s | Cl={cl_end+1}/{dim} alpha={alpha:.6f}")

    # Per-EV Zerlegung
    print(f"\n  {'k':>4s}  {'c_k':>12s}  {'pole=ba':>12s}  {'null':>12s}  "
          f"{'|pole/c|':>10s}  {'|c/b|':>10s}  {'|c/b^2|':>10s}  {'beta_k':>10s}")

    sum_pole2 = 0
    sum_null2 = 0
    sum_cross = 0
    sum_c2 = 0
    sum_beta2 = 0
    first_nc = []

    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]

        bk = float(inner_product(ut, qk, dim))
        ck = float(inner_product(kn, qk, dim))
        gap = float(ws[k] - wmin)

        pole = bk * alpha
        null_part = ck - pole

        if k > cl_end:
            sum_pole2 += pole ** 2
            sum_null2 += null_part ** 2
            sum_cross += pole * null_part
            sum_c2 += ck ** 2
            sum_beta2 += bk ** 2
            if len(first_nc) < 12:
                first_nc.append(dict(k=k, ck=ck, pole=pole, null=null_part,
                                     beta=bk, gap=gap))

        ratio_pc = abs(pole / ck) if abs(ck) > 1e-30 else 0
        ratio_cb = abs(ck / bk) if abs(bk) > 1e-30 else 0
        ratio_cb2 = abs(ck / bk ** 2) if abs(bk) > 1e-20 else 0

        if k <= cl_end + 1 or (k > cl_end and k < cl_end + 12):
            marker = " <-- CL-RAND" if k == cl_end else ""
            print(f"  {k:4d}  {ck:12.4e}  {pole:12.4e}  {null_part:12.4e}  "
                  f"{ratio_pc:10.0f}  {ratio_cb:10.4e}  {ratio_cb2:10.4e}  "
                  f"{bk:10.4e}{marker}")

    print(f"\n  === AGGREGAT (Nicht-Cluster) ===")
    print(f"  Sigma |c_k|^2     = {sum_c2:.4e}  (= Leckage)")
    print(f"  Sigma |pole|^2    = {sum_pole2:.4e}  (= alpha^2 * Sigma beta^2)")
    print(f"  Sigma |null|^2    = {sum_null2:.4e}")
    print(f"  2*Sigma pole*null = {2*sum_cross:.4e}  (Cross-Term)")
    print(f"  Check: pole+null+cross = {sum_pole2+sum_null2+2*sum_cross:.4e}  (soll = c)")
    print(f"  Cancellation: |cross|/(pole+null) = "
          f"{abs(2*sum_cross)/(sum_pole2+sum_null2):.6f}")
    print(f"  alpha = {alpha:.6f}, alpha^2 = {alpha**2:.6f}")
    print(f"  Sigma beta^2 (nc) = {sum_beta2:.6f}")

    return dict(lam=lam, L=L, alpha=alpha, cl=cl_end + 1, dim=dim,
                sum_c2=sum_c2, sum_pole2=sum_pole2, sum_null2=sum_null2,
                sum_cross=sum_cross, sum_beta2=sum_beta2, first_nc=first_nc)


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
    r0, r1 = results[0], results[1]
    ll = np.log(r1["lam"] / r0["lam"])

    print(f"\n{'='*80}")
    print(f"  LAMBDA-SKALIERUNG: D1 vs D2")
    print(f"{'='*80}")

    def exp_fit(v0, v1):
        if abs(v0) > 0 and abs(v1) > 0:
            return np.log(abs(v1) / abs(v0)) / ll
        return float("nan")

    rows = [
        ("Sigma |c|^2 (Leckage)", r0["sum_c2"], r1["sum_c2"]),
        ("Sigma |pole|^2", r0["sum_pole2"], r1["sum_pole2"]),
        ("Sigma |null|^2", r0["sum_null2"], r1["sum_null2"]),
        ("2*cross", abs(2 * r0["sum_cross"]), abs(2 * r1["sum_cross"])),
        ("Sigma |beta|^2 (nc)", r0["sum_beta2"], r1["sum_beta2"]),
        ("alpha^2", r0["alpha"] ** 2, r1["alpha"] ** 2),
    ]

    print(f"\n  {'Metrik':<25s}  {'lam='+str(r0['lam']):>12s}  {'lam='+str(r1['lam']):>12s}  "
          f"{'Exp':>8s}")
    for name, v0, v1 in rows:
        e = exp_fit(v0, v1)
        print(f"  {name:<25s}  {v0:12.4e}  {v1:12.4e}  {e:+8.2f}")

    # Cancellation-Zerlegung der Rate
    exp_c = exp_fit(r0["sum_c2"], r1["sum_c2"])
    exp_p = exp_fit(r0["sum_pole2"], r1["sum_pole2"])
    exp_n = exp_fit(r0["sum_null2"], r1["sum_null2"])
    exp_b = exp_fit(r0["sum_beta2"], r1["sum_beta2"])

    print(f"\n  Rate-Zerlegung:")
    print(f"    Leckage Sigma|c|^2:    {exp_c:+.2f}")
    print(f"    davon D1 (Sigma|b|^2): {exp_b:+.2f}  ({abs(exp_b/exp_c)*100:.0f}% der Rate)")
    print(f"    davon D2 (Cancell.):   {exp_c - exp_b:+.2f}  ({abs((exp_c-exp_b)/exp_c)*100:.0f}% der Rate)")

    # Per-EV Vergleich (1. non-cluster)
    print(f"\n  Cluster-Grenz-EV Vergleich:")
    nc0 = r0["first_nc"][0]
    nc1 = r1["first_nc"][0]
    print(f"    lam={r0['lam']}: k={nc0['k']}, |c|={abs(nc0['ck']):.2e}, "
          f"|pole|={abs(nc0['pole']):.2e}, |beta|={abs(nc0['beta']):.3f}")
    print(f"    lam={r1['lam']}: k={nc1['k']}, |c|={abs(nc1['ck']):.2e}, "
          f"|pole|={abs(nc1['pole']):.2e}, |beta|={abs(nc1['beta']):.3f}")
    print(f"    |c|-Exp:    {exp_fit(abs(nc0['ck']), abs(nc1['ck'])):+.2f}")
    print(f"    |pole|-Exp: {exp_fit(abs(nc0['pole']), abs(nc1['pole'])):+.2f}")
    print(f"    |beta|-Exp: {exp_fit(abs(nc0['beta']), abs(nc1['beta'])):+.2f}")
    print(f"    D2-Beitrag: {exp_fit(abs(nc0['ck']), abs(nc1['ck'])) - exp_fit(abs(nc0['pole']), abs(nc1['pole'])):+.2f}")


if __name__ == "__main__":
    main()
