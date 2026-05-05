"""
c2d_cross_spectral_diagnostic.py — Kreuzspektral-Analyse fuer D2

Diagonalisiert T_lambda = P_0 H P_0 und berechnet:
  nu_{lambda,j} = <e_j, f_lambda> * <e_j, g_lambda>
  f_lambda = P_0 H u_tilde,  g_lambda = P_0 k_full

Testet:
  1. Vorzeichenkohaerenz der nu_j (off-cluster)
  2. Momente M_0, M_1, M_2 der Spektralmasse
  3. Cauchy-Summen-Verifikation: m_lambda(w_k) = sum nu_j/(t_j - w_k)
  4. Finite-Difference dm/dw fuer Flatness
  5. Vergleich mit delta_k = c_k/beta_k

Autor: LG (Opus 4.6, Session 25)
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


def build_operators(lam, N):
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

    Ct = float(sum(W02[i, i] for i in range(dim)))

    return Aq, H, ut, Ct, dim, lam_mp, L_mp


def analyze(lam, N):
    t0 = time.time()
    Aq, H, ut, Ct, dim, lam_mp, L_mp = build_operators(lam, N)
    t_build = time.time() - t0
    print(f"  Build: {t_build:.0f}s", flush=True)

    # --- Diagonalize A_lambda ---
    t0 = time.time()
    ws, Qs = diagonalize_mp(Aq, verbose=False)
    wmin = ws[0]
    t_diag = time.time() - t0
    print(f"  Diag A: {t_diag:.0f}s", flush=True)

    # Cluster
    cl_end = 0
    for k in range(dim):
        if float(ws[k] - wmin) < 1e-10:
            cl_end = k

    # --- H u_tilde ---
    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))

    # --- k_full ---
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    t0 = time.time()
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
    t_proj = time.time() - t0
    print(f"  Proj k: {t_proj:.0f}s", flush=True)

    alpha = float(inner_product(ut, kn, dim))

    # --- Build T_lambda = P_0 H P_0 ---
    # P_0 = I - |ut><ut|
    # T[i,j] = H[i,j] - ut[i]*(Hut)[j] - (Hut)[i]*ut[j] + uHu*ut[i]*ut[j]
    T = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            T[i, j] = H[i, j] - ut[i, 0] * Hut[j, 0] - Hut[i, 0] * ut[j, 0] + mp.mpf(uHu) * ut[i, 0] * ut[j, 0]

    # --- Diagonalize T_lambda ---
    t0 = time.time()
    ts, Es = diagonalize_mp(T, verbose=False)
    t_diag_t = time.time() - t0
    print(f"  Diag T: {t_diag_t:.0f}s", flush=True)

    # --- f_lambda = P_0 H ut = Hut - uHu * ut ---
    f_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]

    # --- g_lambda = P_0 k_full = kn - alpha * ut ---
    g_lam = mp.matrix(dim, 1)
    for i in range(dim):
        g_lam[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]

    # --- Kreuzspektralmasse nu_j ---
    nu = []
    t_vals = []
    for j in range(dim):
        ej = mp.matrix(dim, 1)
        for i in range(dim):
            ej[i, 0] = Es[i, j]

        fj = float(inner_product(ej, f_lam, dim))
        gj = float(inner_product(ej, g_lam, dim))
        tj = float(ts[j])
        nuj = fj * gj
        nu.append(dict(j=j, t=tj, f=fj, g=gj, nu=nuj))
        t_vals.append(tj)

    # Sort by t value
    nu.sort(key=lambda x: x['t'])

    # --- Identify T_lambda null eigenvalue (ut direction) ---
    null_idx = min(range(len(nu)), key=lambda i: abs(nu[i]['t']))
    print(f"  T null-EW: t={nu[null_idx]['t']:.2e} (idx {nu[null_idx]['j']})", flush=True)

    # --- Report nu values ---
    non_null = [n for n in nu if abs(n['t']) > 1e-12]

    # Sign analysis
    pos_nu = sum(1 for n in non_null if n['nu'] > 0)
    neg_nu = sum(1 for n in non_null if n['nu'] < 0)
    total_nu_pos = sum(n['nu'] for n in non_null if n['nu'] > 0)
    total_nu_neg = sum(n['nu'] for n in non_null if n['nu'] < 0)

    print(f"\n  === KREUZSPEKTRALMASS nu_j ===", flush=True)
    print(f"  N_total = {len(non_null)} (excl null-EV)", flush=True)
    print(f"  Positiv: {pos_nu} ({total_nu_pos:+.6e})", flush=True)
    print(f"  Negativ: {neg_nu} ({total_nu_neg:+.6e})", flush=True)
    print(f"  sum(nu) = M_0 = {total_nu_pos + total_nu_neg:.6e} = <f,g>", flush=True)

    # Verify M_0 = <f, g>
    fg_direct = float(inner_product(f_lam, g_lam, dim))
    print(f"  <f,g> direct = {fg_direct:.6e}", flush=True)

    # --- Moments ---
    t_center = np.median([n['t'] for n in non_null])
    M0 = sum(n['nu'] for n in non_null)
    M1 = sum((n['t'] - t_center) * n['nu'] for n in non_null)
    M2 = sum((n['t'] - t_center)**2 * n['nu'] for n in non_null)
    total_var = sum(abs(n['nu']) for n in non_null)

    print(f"\n  === MOMENTE (Zentrum t* = {t_center:.6f}) ===", flush=True)
    print(f"  M_0 = {M0:.6e}", flush=True)
    print(f"  M_1 = {M1:.6e}", flush=True)
    print(f"  M_2 = {M2:.6e}", flush=True)
    print(f"  TV  = {total_var:.6e} (Totalvariation)", flush=True)
    print(f"  |M_0|/TV = {abs(M0)/total_var:.6e}", flush=True)
    print(f"  |M_1|/TV = {abs(M1)/total_var:.6e}", flush=True)

    # --- Print first ~20 nu values ---
    print(f"\n  {'j':>4s}  {'t_j':>12s}  {'<e,f>':>12s}  {'<e,g>':>12s}  "
          f"{'nu_j':>12s}  {'sign':>5s}", flush=True)
    for idx, n in enumerate(non_null[:25]):
        sign_str = "+" if n['nu'] > 0 else "-"
        print(f"  {n['j']:4d}  {n['t']:12.6f}  {n['f']:12.4e}  {n['g']:12.4e}  "
              f"{n['nu']:12.4e}  {sign_str:>5s}", flush=True)
    if len(non_null) > 25:
        print(f"  ... ({len(non_null) - 25} more)", flush=True)

    # --- Cauchy-Summen-Verifikation ---
    print(f"\n  === CAUCHY-SUMME vs delta_k ===", flush=True)
    print(f"  {'k':>4s}  {'m_cauchy':>12s}  {'m_exact':>12s}  {'diff':>12s}  "
          f"{'delta_k':>12s}  {'zone':>8s}", flush=True)

    for k in range(dim):
        if k > cl_end + 6 and k < dim - 2:
            continue

        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]

        beta_k = float(inner_product(ut, qk, dim))
        c_k = float(inner_product(kn, qk, dim))
        wk = float(ws[k])

        if abs(beta_k) > 1e-15:
            delta_exact = c_k / beta_k
            m_exact = alpha - delta_exact
        else:
            delta_exact = float('nan')
            m_exact = float('nan')

        # Cauchy sum
        m_cauchy = sum(n['nu'] / (n['t'] - wk) for n in non_null
                       if abs(n['t'] - wk) > 1e-20)

        zone = "CLUSTER" if k <= cl_end else "OFF"

        if np.isnan(m_exact):
            print(f"  {k:4d}  {m_cauchy:12.6f}  {'NaN':>12s}  {'NaN':>12s}  "
                  f"{'NaN':>12s}  {zone:>8s}", flush=True)
        else:
            diff = abs(m_cauchy - m_exact)
            print(f"  {k:4d}  {m_cauchy:12.6f}  {m_exact:12.6f}  {diff:12.4e}  "
                  f"{delta_exact:12.4e}  {zone:>8s}", flush=True)

    # --- Finite Differences ---
    print(f"\n  === FINITE DIFFERENCES dm/dw (Off-Cluster) ===", flush=True)
    off_data = []
    for k in range(cl_end + 1, dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]

        beta_k = float(inner_product(ut, qk, dim))
        c_k = float(inner_product(kn, qk, dim))
        wk = float(ws[k])

        if abs(beta_k) > 1e-15:
            delta_k = c_k / beta_k
            m_k = alpha - delta_k
            off_data.append(dict(k=k, w=wk, m=m_k, delta=delta_k))

    print(f"  {'k':>4s}  {'w_k':>12s}  {'m_k':>12s}  {'dm/dw':>12s}  "
          f"{'delta_k':>12s}", flush=True)
    for idx in range(len(off_data)):
        d = off_data[idx]
        if idx > 0:
            dw = d['w'] - off_data[idx-1]['w']
            dm = d['m'] - off_data[idx-1]['m']
            dmdw = dm / dw if abs(dw) > 1e-20 else float('nan')
        else:
            dmdw = float('nan')
        if idx < 8 or idx >= len(off_data) - 3:
            print(f"  {d['k']:4d}  {d['w']:12.6f}  {d['m']:12.8f}  "
                  f"{dmdw:12.4e}  {d['delta']:12.4e}", flush=True)

    return dict(
        lam=lam, alpha=alpha, cl_end=cl_end, dim=dim,
        nu=non_null, M0=M0, M1=M1, M2=M2, TV=total_var,
        pos_count=pos_nu, neg_count=neg_nu,
    )


def main():
    results = []
    for cfg in CONFIGS:
        lam, N = cfg["lam"], cfg["N"]
        print(f"\n{'='*80}", flush=True)
        print(f"  LAMBDA = {lam}, N = {N}", flush=True)
        print(f"{'='*80}", flush=True)
        r = analyze(lam, N)
        results.append(r)

    if len(results) < 2:
        return

    # --- Multi-Lambda ---
    print(f"\n{'='*80}", flush=True)
    print(f"  MULTI-LAMBDA KREUZSPEKTRAL-SKALIERUNG", flush=True)
    print(f"{'='*80}", flush=True)

    lams = [r['lam'] for r in results]

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

    metrics = [
        ("M_0", [r['M0'] for r in results]),
        ("|M_0|", [abs(r['M0']) for r in results]),
        ("M_1", [r['M1'] for r in results]),
        ("|M_1|", [abs(r['M1']) for r in results]),
        ("M_2", [r['M2'] for r in results]),
        ("TV", [r['TV'] for r in results]),
        ("|M_0|/TV", [abs(r['M0'])/r['TV'] for r in results]),
        ("pos/neg ratio", [r['pos_count']/max(r['neg_count'],1) for r in results]),
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


if __name__ == "__main__":
    main()
