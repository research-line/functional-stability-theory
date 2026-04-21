"""
c2_cancellation_lambda_scan.py — Lambda-Skalierung der Drei-Ebenen-Cancellation

Scannt b_k/a_k Ratio, u_tilde-Kopplung, Pol/Null-Split ueber mehrere lambda.
Ziel: Wie skaliert |1+b_k/a_k| mit lambda? Gap-Abhaengigkeit?

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
    k_lambda_value, h_educated_guess, inner_product, norm, evaluate_F_even,
)
from c2_poisson_decomposition import project_to_fourier

DPS = int(os.environ.get("DPS", 50))
GAMMA1 = 14.134725141734693

CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
]


def analyze(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

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

    t0 = time.time()

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    def kn1(x):
        return mp.sqrt(mp.exp(x)) * h_educated_guess(mp.exp(x))

    cf = project_to_fourier(kfull, L_mp, N)
    cn1 = project_to_fourier(kn1, L_mp, N)
    cn2 = mp.matrix(dim, 1)
    for i in range(dim):
        cn2[i, 0] = cf[i, 0] - cn1[i, 0]

    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    p1 = mp.matrix(dim, 1)
    p2 = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
        p1[i, 0] = cn1[i, 0] / nf
        p2[i, 0] = cn2[i, 0] / nf
    t_proj = time.time() - t0

    gamma = mp.mpf(GAMMA1)
    eg = mp.matrix(dim, 1)
    eg[0, 0] = 1 / gamma
    for m in range(1, dim):
        pole = 2 * mp.pi * m / L_mp
        eg[m, 0] = mp.sqrt(mp.mpf(2)) * gamma / (gamma**2 - pole**2)
    alpha_g = float(inner_product(ut, eg, dim))

    cl_end = 0
    for k in range(dim):
        if float(ws[k] - wmin) < 1e-10:
            cl_end = k

    ratios = []
    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]
        ak = float(inner_product(p1, qk, dim))
        bk = float(inner_product(p2, qk, dim))
        ck = float(inner_product(kn, qk, dim))
        bku = float(inner_product(ut, qk, dim))
        gap = float(ws[k] - wmin)
        dev = abs(ck / ak) if abs(ak) > 1e-20 else 0
        ratios.append(dict(k=k, a=ak, b=bk, c=ck, gap=gap, beta=bku, dev=dev))

    un1 = float(inner_product(ut, p1, dim))
    un2 = float(inner_product(ut, p2, dim))
    n1sq = float(inner_product(p1, p1, dim))
    n2sq = float(inner_product(p2, p2, dim))

    Fn1 = float(inner_product(p1, eg, dim))
    Fn2 = float(inner_product(p2, eg, dim))
    Ff = float(inner_product(kn, eg, dim))

    nc = [r for r in ratios if r["k"] > cl_end]
    max_dev = max(r["dev"] for r in nc) if nc else 0
    rms_dev = np.sqrt(np.mean([r["dev"] ** 2 for r in nc])) if nc else 0
    leak = sum(r["c"] ** 2 for r in nc)

    return dict(
        lam=lam, N=N, dim=dim, wmin=float(wmin), Ct=Ct,
        cl=cl_end + 1, t_qw=t_qw, t_proj=t_proj,
        ratios=ratios, max_dev=max_dev, rms_dev=rms_dev, leak=leak,
        un1=un1, un2=un2, n1sq=n1sq, n2sq=n2sq,
        pf1=un1 ** 2 / n1sq * 100 if n1sq > 0 else 0,
        pf2=un2 ** 2 / n2sq * 100 if n2sq > 0 else 0,
        Fn1=Fn1, Fn2=Fn2, Ff=Ff, alpha_g=alpha_g,
        Fp1=alpha_g * un1, Fp2=alpha_g * un2,
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
        print(f"  QW={r['t_qw']:.0f}s Proj={r['t_proj']:.0f}s "
              f"| w_min={r['wmin']:.4f} C~={r['Ct']:.2f} Cl={r['cl']}/{r['dim']}")
        print(f"  Max|1+b/a|={r['max_dev']:.2e} RMS={r['rms_dev']:.2e} Leak={r['leak']:.2e}")

        print(f"\n  Nicht-Cluster-Ratios (erste 10):")
        print(f"  {'k':>4s}  {'|1+b/a|':>10s}  {'gap':>10s}  {'|c_k|^2':>10s}  {'|beta|':>8s}")
        for rr in r["ratios"][r["cl"]:r["cl"] + 10]:
            print(f"  {rr['k']:4d}  {rr['dev']:10.2e}  {rr['gap']:10.2e}  "
                  f"{rr['c']**2:10.2e}  {abs(rr['beta']):8.2e}")

        print(f"\n  Pol-Kopplung: un1={r['un1']:.4f}({r['pf1']:.1f}%) "
              f"un2={r['un2']:.4f}({r['pf2']:.1f}%)")
        print(f"  F(g1): Fn1={r['Fn1']:.4e} Fn2={r['Fn2']:.4e} Ff={r['Ff']:.4e}")
        print(f"         Fp1={r['Fp1']:.4e} Fp2={r['Fp2']:.4e} "
              f"Fpol={r['Fp1']+r['Fp2']:.4e}")

    # === VERGLEICH ===
    if len(results) < 2:
        return
    r0, r1 = results[0], results[1]
    ll = np.log(r1["lam"] / r0["lam"])

    print(f"\n{'='*80}")
    print(f"  LAMBDA-SKALIERUNG (lam {r0['lam']} -> {r1['lam']})")
    print(f"{'='*80}")

    def exp_fit(v0, v1):
        if v0 > 0 and v1 > 0:
            return np.log(v1 / v0) / ll
        return float("nan")

    rows = [
        ("Max |1+b/a|", r0["max_dev"], r1["max_dev"]),
        ("RMS |1+b/a|", r0["rms_dev"], r1["rms_dev"]),
        ("Leakage ||c||^2", r0["leak"], r1["leak"]),
        ("|F_full(g1)|", abs(r0["Ff"]), abs(r1["Ff"])),
        ("|F_n1|+|F_n2|", abs(r0["Fn1"]) + abs(r0["Fn2"]),
         abs(r1["Fn1"]) + abs(r1["Fn2"])),
    ]
    print(f"\n  {'Metrik':<22s}  {'lam='+str(r0['lam']):>10s}  {'lam='+str(r1['lam']):>10s}  "
          f"{'Exponent':>10s}")
    for name, v0, v1 in rows:
        e = exp_fit(v0, v1)
        print(f"  {name:<22s}  {v0:10.2e}  {v1:10.2e}  {e:+10.2f}")

    # Gap-Korrelation
    print(f"\n  Gap-Korrelation (log|1+b/a| vs log gap, non-cluster, gap>1e-6):")
    for r in results:
        nc = [rr for rr in r["ratios"]
              if rr["k"] > r["cl"] and rr["gap"] > 1e-6 and rr["dev"] > 1e-20]
        if len(nc) < 3:
            print(f"    lam={r['lam']}: zu wenig Datenpunkte ({len(nc)})")
            continue
        lg = np.array([np.log10(rr["gap"]) for rr in nc])
        ld = np.array([np.log10(rr["dev"]) for rr in nc])
        slope, intercept = np.polyfit(lg, ld, 1)
        print(f"    lam={r['lam']}: slope={slope:+.3f}, intercept={intercept:.2f} "
              f"(n={len(nc)} pts, gap range [{nc[0]['gap']:.1e}, {nc[-1]['gap']:.1e}])")

    # Cluster-Grenze: Vergleich des ersten non-cluster EV
    print(f"\n  Cluster-Grenz-Vergleich:")
    for r in results:
        k0 = r["cl"]
        if k0 < r["dim"]:
            rr = r["ratios"][k0]
            print(f"    lam={r['lam']}: k={k0}, gap={rr['gap']:.2e}, "
                  f"|1+b/a|={rr['dev']:.2e}, |c_k|={abs(rr['c']):.2e}")


if __name__ == "__main__":
    main()
