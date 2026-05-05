"""
c2bb_poisson_rayleigh_cancellation.py

Exakte algebraische Identitaet:

  mu/alpha = (w_bar - w_min) + <f_lam, k_perp>/alpha

wobei w_bar = C_tilde + <u_tilde, H u_tilde> (Rayleigh-Quotient von u_tilde fuer A).

Vier Schluesselgroessen pro lambda:
  1. (w_bar - w_min)           Rayleigh-Gap
  2. <f_lam, k_perp>/alpha     Poisson-Alignment
  3. mu/alpha                  Summe (= erster Off-Cluster-Moment)
  4. Cancellation-Ratio        |mu/alpha| / (|T1| + |T2|)

Secular-Identitaet: 1 + C_tilde * m_H(w_min) = 0
  mit m_H(w) = sum beta_i^2 / (h_i - w), H-Eigenwerte h_i, beta_i = <u_tilde, phi_i>.
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import mpmath as mp
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
    {"lam": 7.0, "N": 80},
]


def run_one(lam, N):
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
    H = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            W02[i, j] = Aq[i, j] - Ah[i, j]
            H[i, j] = Ah[i, j]

    # u_tilde
    col0 = mp.matrix(dim, 1)
    for i in range(dim):
        col0[i, 0] = W02[i, 0]
    cn = norm(col0, dim)
    ut = mp.matrix(dim, 1)
    for i in range(dim):
        ut[i, 0] = col0[i, 0] / cn

    # C_tilde, uHu
    Ctilde = float(sum(ut[i, 0] * sum(W02[i, j] * ut[j, 0] for j in range(dim))
                       for i in range(dim)))

    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))

    w_bar = Ctilde + uHu

    # k_lambda
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
    alpha = float(inner_product(ut, kn, dim))

    # k_perp = k - alpha u_tilde
    k_perp = mp.matrix(dim, 1)
    for i in range(dim):
        k_perp[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]
    k_perp_norm = float(norm(k_perp, dim))

    # f_lambda = P_0 H u_tilde = H u_tilde - uHu * u_tilde
    f_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]
    f_lam_norm = float(norm(f_lam, dim))

    # A-Diagonalisierung -> w_min
    ws_A, Qs = diagonalize_mp(Aq, verbose=False)
    w_min = float(ws_A[0])

    # H-Diagonalisierung -> Secular-Identitaet
    hs, Ps = diagonalize_mp(H, verbose=False)

    beta_sq = np.zeros(dim)
    for i_ev in range(dim):
        phi_i = mp.matrix(dim, 1)
        for r in range(dim):
            phi_i[r, 0] = Ps[r, i_ev]
        beta_sq[i_ev] = float(inner_product(ut, phi_i, dim))**2

    h_arr = np.array([float(hs[i]) for i in range(dim)])

    build_time = time.time() - t0

    # === Exakte Identitaet ===
    # mu = alpha * (w_bar - w_min) + <f_lam, k_perp>
    T1 = w_bar - w_min
    T2 = float(inner_product(f_lam, k_perp, dim)) / alpha
    mu_alpha = T1 + T2

    # Spektrale Kontrolle
    c_arr = np.zeros(dim)
    alpha_a = np.zeros(dim)
    w_arr = np.array([float(ws_A[a]) for a in range(dim)])
    for a in range(dim):
        qa = mp.matrix(dim, 1)
        for r in range(dim):
            qa[r, 0] = Qs[r, a]
        c_arr[a] = float(inner_product(qa, kn, dim))
        alpha_a[a] = float(inner_product(ut, qa, dim))

    cl_A = [a for a in range(dim) if abs(w_arr[a] - w_min) < 1e-10]
    noncl_A = [a for a in range(dim) if a not in cl_A]
    mu_spectral = sum((w_arr[a] - w_min) * alpha_a[a] * c_arr[a] for a in range(dim))
    mu_alpha_spectral = mu_spectral / alpha

    # Secular-Identitaet mit H-Eigenwerten
    m_H_wmin = sum(beta_sq[i] / (h_arr[i] - w_min) for i in range(dim)
                   if abs(h_arr[i] - w_min) > 1e-30)
    secular = 1 + Ctilde * m_H_wmin

    # Cancellation-Ratio
    denom = abs(T1) + abs(T2)
    cancel_ratio = abs(mu_alpha) / denom if denom > 1e-50 else float('inf')

    return {
        'lam': lam, 'N': N, 'dim': dim, 'build_time': build_time,
        'Ctilde': Ctilde, 'uHu': uHu, 'w_bar': w_bar, 'w_min': w_min,
        'alpha': alpha, 'k_perp_norm': k_perp_norm,
        'f_lam_norm': f_lam_norm,
        'T1': T1, 'T2': T2, 'mu_alpha': mu_alpha,
        'mu_alpha_spectral': mu_alpha_spectral,
        'cancel_ratio': cancel_ratio,
        'secular': secular,
        'n_cluster': len(cl_A),
    }


def main():
    print(f"C2bb: Poisson-Rayleigh Cancellation Diagnostic", flush=True)
    print(f"DPS = {DPS}", flush=True)

    results = []
    for cfg in CONFIGS:
        lam, N = cfg['lam'], cfg['N']
        print(f"\n{'='*100}", flush=True)
        print(f"  lambda = {lam}, N = {N}", flush=True)
        print(f"{'='*100}", flush=True)

        r = run_one(lam, N)
        results.append(r)

        print(f"  Build: {r['build_time']:.0f}s", flush=True)
        print(f"  dim={r['dim']}, Cluster={r['n_cluster']}", flush=True)
        print(f"", flush=True)
        print(f"  C_tilde          = {r['Ctilde']:.6f}", flush=True)
        print(f"  <u,Hu>           = {r['uHu']:.6f}", flush=True)
        print(f"  w_bar            = {r['w_bar']:.6f}", flush=True)
        print(f"  w_min            = {r['w_min']:.6f}", flush=True)
        print(f"  alpha            = {r['alpha']:.6f}", flush=True)
        print(f"  ||k_perp||       = {r['k_perp_norm']:.6e}", flush=True)
        print(f"  ||f_lam||        = {r['f_lam_norm']:.6e}", flush=True)
        print(f"", flush=True)
        print(f"  --- Exakte Identitaet ---", flush=True)
        print(f"  T1: (w_bar-w_min)       = {r['T1']:+.10e}", flush=True)
        print(f"  T2: <f,k_perp>/alpha    = {r['T2']:+.10e}", flush=True)
        print(f"  mu/alpha (Summe)        = {r['mu_alpha']:+.10e}", flush=True)
        print(f"  mu/alpha (spektral)     = {r['mu_alpha_spectral']:+.10e}", flush=True)
        print(f"  Konsistenz              = {abs(r['mu_alpha'] - r['mu_alpha_spectral']):.2e}", flush=True)
        print(f"  Cancellation-Ratio      = {r['cancel_ratio']:.6e}", flush=True)
        print(f"  Cancellation-Faktor     = {1/r['cancel_ratio']:.0f}x", flush=True)
        print(f"", flush=True)
        print(f"  Secular: 1+C*m_H(w_min) = {r['secular']:.6e}  (soll=0)", flush=True)

    # === Cross-lambda Tabelle ===
    print(f"\n{'='*100}", flush=True)
    print(f"CROSS-LAMBDA ZUSAMMENFASSUNG", flush=True)
    print(f"{'='*100}", flush=True)
    print(f"", flush=True)
    print(f"  {'lam':>5}  {'T1=w_bar-w_min':>16}  {'T2=<f,k_p>/a':>16}  "
          f"{'mu/alpha':>14}  {'Cancel':>12}  {'Factor':>8}  {'Secular':>12}", flush=True)
    print(f"  {'-'*95}", flush=True)
    for r in results:
        print(f"  {r['lam']:5.1f}  {r['T1']:+16.10e}  {r['T2']:+16.10e}  "
              f"{r['mu_alpha']:+14.6e}  {r['cancel_ratio']:12.4e}  "
              f"{1/r['cancel_ratio']:8.0f}x  {r['secular']:+12.4e}", flush=True)

    print(f"\n  {'lam':>5}  {'alpha':>8}  {'||k_perp||':>12}  {'||f_lam||':>12}  "
          f"{'Cluster':>7}  {'dim':>5}", flush=True)
    print(f"  {'-'*60}", flush=True)
    for r in results:
        print(f"  {r['lam']:5.1f}  {r['alpha']:8.4f}  {r['k_perp_norm']:12.4e}  "
              f"{r['f_lam_norm']:12.4e}  {r['n_cluster']:7d}  {r['dim']:5d}", flush=True)

    # Skalierung
    if len(results) >= 2:
        print(f"\n  --- Skalierung ---", flush=True)
        for i in range(1, len(results)):
            r0, r1 = results[i-1], results[i]
            lam_ratio = r1['lam'] / r0['lam']
            t1_ratio = abs(r1['T1']) / abs(r0['T1'])
            cancel_ratio = r1['cancel_ratio'] / r0['cancel_ratio']
            mu_ratio = abs(r1['mu_alpha']) / abs(r0['mu_alpha'])
            print(f"  lam {r0['lam']}->{r1['lam']}: "
                  f"T1 x{t1_ratio:.3f}  "
                  f"Cancel x{cancel_ratio:.3e}  "
                  f"|mu/a| x{mu_ratio:.3e}", flush=True)

    print(f"\n  BEWERTUNG:", flush=True)
    all_cancel = all(r['cancel_ratio'] < 0.01 for r in results)
    if all_cancel:
        print(f"  STABIL: Cancellation < 1% bei allen lambda.", flush=True)
        print(f"  => C2bb Poisson-Rayleigh Cancellation Lemma ist numerisch gestuetzt.", flush=True)
    else:
        fails = [r for r in results if r['cancel_ratio'] >= 0.01]
        print(f"  INSTABIL bei lambda = {[r['lam'] for r in fails]}", flush=True)


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
