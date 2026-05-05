"""
c2d_delta_resolvent_test.py — Test der D2-Faktorisierung c_k = beta_k * delta_k

Exakte Beziehung:
  c_k = beta_k * (alpha_lambda - m_lambda(w_k))
  m_lambda(w) := <(T_lambda - wI)^{-1} f_lambda, P_0 k_full>
  delta_k := alpha_lambda - m_lambda(w_k) = c_k / beta_k

Testet: Ist m_lambda(w_k) ~ alpha fuer alle Off-Cluster k?
        Traegt sup|delta_k| den D2-Exponent ~lambda^{-2.67}?

Autor: LG (Opus 4.6, Session 25)
Datum: 2026-04-19
"""

import os, sys, time, json
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

    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf

    alpha = float(inner_product(ut, kn, dim))

    cl_end = 0
    for k in range(dim):
        if float(ws[k] - wmin) < 1e-10:
            cl_end = k

    all_data = []
    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]

        beta_k = float(inner_product(ut, qk, dim))
        c_k = float(inner_product(kn, qk, dim))
        wk = float(ws[k])
        gap = wk - float(wmin)

        if abs(beta_k) > 1e-15:
            delta_k = c_k / beta_k
            m_k = alpha - delta_k
        else:
            delta_k = float('nan')
            m_k = float('nan')

        all_data.append(dict(
            k=k, beta=beta_k, c=c_k, delta=delta_k, m=m_k,
            w=wk, gap=gap, is_cl=(k <= cl_end)
        ))

    return dict(
        lam=lam, alpha=alpha, cl_end=cl_end, dim=dim,
        Ct=Ct, uHu=uHu, data=all_data
    )


def exp_fit(vals, lams):
    pairs = [(l, v) for l, v in zip(lams, vals) if abs(v) > 0]
    if len(pairs) < 2:
        return float('nan')
    log_l = [np.log(l) for l, _ in pairs]
    log_v = [np.log(abs(v)) for _, v in pairs]
    n = len(log_v)
    sx, sy = sum(log_l), sum(log_v)
    sxx = sum(x**2 for x in log_l)
    sxy = sum(x * y for x, y in zip(log_l, log_v))
    return (n * sxy - sx * sy) / (n * sxx - sx**2)


def main():
    results = []
    for cfg in CONFIGS:
        lam, N = cfg["lam"], cfg["N"]
        print(f"\n{'='*80}", flush=True)
        print(f"  LAMBDA = {lam}, N = {N}", flush=True)
        print(f"{'='*80}", flush=True)
        t0 = time.time()
        r = analyze(lam, N)
        dt = time.time() - t0

        print(f"  Time: {dt:.0f}s", flush=True)
        print(f"  alpha = {r['alpha']:.8f}", flush=True)
        print(f"  Cluster: k=0..{r['cl_end']} ({r['cl_end']+1}/{r['dim']})", flush=True)

        # --- Per-k Tabelle ---
        print(f"\n  {'k':>4s}  {'beta_k':>12s}  {'c_k':>12s}  {'delta_k':>12s}  "
              f"{'m_k':>12s}  {'|delta|':>10s}  {'gap':>10s}  {'zone':>8s}", flush=True)

        for d in r['data']:
            zone = "CLUSTER" if d['is_cl'] else "OFF"
            show = (d['k'] <= r['cl_end'] + 8
                    or d['k'] >= r['dim'] - 3
                    or d['k'] == r['cl_end'])
            if not show:
                continue
            if np.isnan(d['delta']):
                print(f"  {d['k']:4d}  {d['beta']:12.4e}  {d['c']:12.4e}  "
                      f"{'NaN':>12s}  {'NaN':>12s}  {'NaN':>10s}  "
                      f"{d['gap']:10.2e}  {zone:>8s}", flush=True)
            else:
                print(f"  {d['k']:4d}  {d['beta']:12.4e}  {d['c']:12.4e}  "
                      f"{d['delta']:12.4e}  {d['m']:12.8f}  {abs(d['delta']):10.4e}  "
                      f"{d['gap']:10.2e}  {zone:>8s}", flush=True)

        # --- Off-Cluster Statistik ---
        off = [d for d in r['data'] if not d['is_cl'] and not np.isnan(d['delta'])]
        off_deltas = [abs(d['delta']) for d in off]
        off_deltas2 = [d['delta']**2 for d in off]

        if off_deltas:
            sup_d = max(off_deltas)
            rms_d = np.sqrt(np.mean(off_deltas2))
            mean_d = np.mean(off_deltas)
            weighted_d2 = sum(d['beta']**2 * d['delta']**2 for d in off)
            sum_c2 = sum(d['c']**2 for d in off)
            sum_b2 = sum(d['beta']**2 for d in off)

            print(f"\n  === OFF-CLUSTER DELTA STATISTIK ===", flush=True)
            print(f"  N_off = {len(off_deltas)}", flush=True)
            print(f"  sup|delta_k|   = {sup_d:.6e}", flush=True)
            print(f"  RMS|delta_k|   = {rms_d:.6e}", flush=True)
            print(f"  mean|delta_k|  = {mean_d:.6e}", flush=True)
            print(f"  min|delta_k|   = {min(off_deltas):.6e}", flush=True)
            print(f"  sum|c_k|^2     = {sum_c2:.6e}  (= sum|beta|^2 * |delta|^2)", flush=True)
            print(f"  sum|beta_k|^2  = {sum_b2:.6e}", flush=True)
            print(f"  eff delta^2    = {sum_c2/sum_b2:.6e}  (=sum|c|^2/sum|beta|^2)", flush=True)

            r['sup_delta'] = sup_d
            r['rms_delta'] = rms_d
            r['mean_delta'] = mean_d
            r['eff_delta2'] = sum_c2 / sum_b2
            r['sum_c2'] = sum_c2
            r['sum_b2'] = sum_b2

        results.append(r)

    # --- Multi-Lambda Skalierung ---
    if len(results) < 2:
        return

    print(f"\n{'='*80}", flush=True)
    print(f"  MULTI-LAMBDA DELTA SKALIERUNG", flush=True)
    print(f"{'='*80}", flush=True)

    lams = [r['lam'] for r in results]

    metrics = [
        ("alpha",          [r['alpha'] for r in results]),
        ("sup|delta|",     [r.get('sup_delta', 0) for r in results]),
        ("sup|delta|^2",   [r.get('sup_delta', 0)**2 for r in results]),
        ("RMS|delta|",     [r.get('rms_delta', 0) for r in results]),
        ("RMS|delta|^2",   [r.get('rms_delta', 0)**2 for r in results]),
        ("eff delta^2",    [r.get('eff_delta2', 0) for r in results]),
        ("sum|c|^2",       [r.get('sum_c2', 0) for r in results]),
        ("sum|beta|^2",    [r.get('sum_b2', 0) for r in results]),
    ]

    header = f"\n  {'Metrik':<18s}"
    for r in results:
        header += f"  {'lam='+str(r['lam']):>12s}"
    header += f"  {'Exp(fit)':>10s}"
    print(header, flush=True)

    for name, vals in metrics:
        line = f"  {name:<18s}"
        for v in vals:
            line += f"  {v:12.6e}"
        exp = exp_fit(vals, lams)
        line += f"  {exp:+10.4f}"
        print(line, flush=True)

    print(f"\n  === D2-INTERPRETATION ===", flush=True)
    print(f"  Erwarteter D2-Exponent: delta^2 ~ lambda^{{-5.35}}", flush=True)
    print(f"  Also: delta ~ lambda^{{-2.675}}", flush=True)
    eff2_exp = exp_fit([r.get('eff_delta2', 0) for r in results], lams)
    sup2_exp = exp_fit([r.get('sup_delta', 0)**2 for r in results], lams)
    rms2_exp = exp_fit([r.get('rms_delta', 0)**2 for r in results], lams)
    c2_exp = exp_fit([r.get('sum_c2', 0) for r in results], lams)
    b2_exp = exp_fit([r.get('sum_b2', 0) for r in results], lams)
    print(f"  eff delta^2 Exp:  {eff2_exp:+.4f}  (sum|c|^2/sum|beta|^2)", flush=True)
    print(f"  sup|delta|^2 Exp: {sup2_exp:+.4f}", flush=True)
    print(f"  RMS|delta|^2 Exp: {rms2_exp:+.4f}", flush=True)
    print(f"  sum|c|^2 Exp:     {c2_exp:+.4f}", flush=True)
    print(f"  sum|beta|^2 Exp:  {b2_exp:+.4f}", flush=True)
    print(f"  D2 = sum|c|^2 - sum|beta|^2: {c2_exp:+.4f} - ({b2_exp:+.4f}) = {c2_exp - b2_exp:+.4f}", flush=True)

    # Herglotz-Profil: m_lambda(w_k) ueber k fuer jedes lambda
    print(f"\n  === HERGLOTZ-PROFIL: m_lambda(w_k) vs k ===", flush=True)
    for r in results:
        print(f"\n  lambda = {r['lam']}, alpha = {r['alpha']:.6f}", flush=True)
        print(f"  {'k':>4s}  {'m_k':>12s}  {'|m_k-alpha|':>12s}  {'zone':>8s}", flush=True)
        sel = [d for d in r['data']
               if d['k'] <= r['cl_end'] + 6
               or d['k'] in [r['cl_end'] - 1, r['cl_end'], r['cl_end'] + 1]]
        for d in sel:
            if np.isnan(d['m']):
                continue
            zone = "CLUSTER" if d['is_cl'] else "OFF"
            diff = abs(d['m'] - r['alpha'])
            marker = " <-- CL-RAND" if d['k'] == r['cl_end'] else ""
            print(f"  {d['k']:4d}  {d['m']:12.8f}  {diff:12.4e}  {zone:>8s}{marker}", flush=True)


if __name__ == "__main__":
    main()
