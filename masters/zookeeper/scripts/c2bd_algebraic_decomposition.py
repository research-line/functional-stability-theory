"""
c2bd_algebraic_decomposition.py

Algebraische Zerlegung des Derivative Freezing.

Saubere Form:
  mu/alpha = <u|(A - w0)|k> / alpha = sigma - w0 + Ctilde

wobei sigma = <u|H|k> / alpha.

Daher:
  d/dL(mu/alpha) = sigma' - w0' + Ctilde'

Dieses Script berechnet alle drei Terme via finite Differenzen
und identifiziert den Cancellation-Mechanismus.
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import mpmath as mp
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, chi_trivial,
)
from c2_approximation_test import (
    h_educated_guess, k_lambda_value, norm,
)
from c2_poisson_decomposition import (
    h_hat_analytical, project_to_fourier,
)

DPS = int(os.environ.get("DPS", 25))
N_DEFAULT = 20


def to_np(M, dim):
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = float(M[i, j])
    return A


def compute_algebraic_data(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    A_np = to_np(Aq, dim)
    H_np = to_np(Ah, dim)
    W_np = A_np - H_np

    col0 = W_np[:, 0]
    cn = np.linalg.norm(col0)
    ut = col0 / cn

    Ct = float(ut @ W_np @ ut)
    uHu = float(ut @ H_np @ ut)

    ws_A, vs_A = np.linalg.eigh(A_np)
    w0 = ws_A[0]
    g0 = vs_A[:, 0]
    if g0[0] < 0:
        g0 = -g0

    hs, vs_H = np.linalg.eigh(H_np)
    betas = vs_H.T @ ut

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
    alpha = float(ut @ kn)

    # sigma = <u|H|k> / alpha
    uHk = float(ut @ H_np @ kn)
    sigma = uHk / alpha

    # Verify: mu/alpha = sigma - w0 + Ct
    mu_alpha_formula = sigma - w0 + Ct
    mu_alpha_direct = float(ut @ (A_np - w0 * np.eye(dim)) @ kn) / alpha

    # T1 = w_bar - w0, T2 via f_lam
    w_bar = Ct + uHu
    T1 = w_bar - w0
    k_perp = kn - alpha * ut
    f_lam = H_np @ ut - uHu * ut
    T2 = float(f_lam @ k_perp) / alpha
    mu_alpha_T1T2 = T1 + T2

    # Alignment: <g0|k> and <g0|u>
    g0k = float(g0 @ kn)
    g0u = float(g0 @ ut)

    # Stieltjes
    m_H = sum(b**2 / (h - w0) for b, h in zip(betas, hs) if abs(h - w0) > 1e-30)
    m_H_prime = sum(b**2 / (h - w0)**2 for b, h in zip(betas, hs) if abs(h - w0) > 1e-30)
    secular_check = 1 + Ct * m_H

    return {
        'lam': lam, 'L': L,
        'sigma': sigma, 'w0': w0, 'Ct': Ct,
        'mu_alpha_formula': mu_alpha_formula,
        'mu_alpha_direct': mu_alpha_direct,
        'mu_alpha_T1T2': mu_alpha_T1T2,
        'alpha': alpha, 'uHu': uHu,
        'T1': T1, 'T2': T2,
        'g0k': g0k, 'g0u': g0u,
        'm_H': m_H, 'm_H_prime': m_H_prime,
        'secular_check': secular_check,
        'hs': list(hs), 'betas': list(betas),
    }


def run_decomposition(lam_base, N, dlam=0.01):
    L_base = 2 * np.log(lam_base)
    dL = 2 * dlam / lam_base

    print(f"\n{'='*90}")
    print(f"  ALGEBRAISCHE ZERLEGUNG: lam={lam_base}, N={N}, dlam={dlam}")
    print(f"  L={L_base:.6f}, dL={dL:.6e}")
    print(f"{'='*90}")

    t0 = time.time()
    d_m = compute_algebraic_data(lam_base - dlam, N)
    d_0 = compute_algebraic_data(lam_base, N)
    d_p = compute_algebraic_data(lam_base + dlam, N)
    print(f"  3 Punkte in {time.time()-t0:.0f}s", flush=True)

    # Verify clean form
    print(f"\n  --- VERIFIKATION: mu/alpha = sigma - w0 + Ct ---")
    print(f"  sigma       = {d_0['sigma']:+.10e}")
    print(f"  w0          = {d_0['w0']:+.10e}")
    print(f"  Ct          = {d_0['Ct']:+.10e}")
    print(f"  sigma-w0+Ct = {d_0['mu_alpha_formula']:+.10e}")
    print(f"  <u|(A-w0)|k>/a = {d_0['mu_alpha_direct']:+.10e}")
    print(f"  T1+T2       = {d_0['mu_alpha_T1T2']:+.10e}")
    print(f"  1+Ct*mH     = {d_0['secular_check']:+.6e}  (sollte ~0)")

    # Finite differences of the three components
    sigma_prime = (d_p['sigma'] - d_m['sigma']) / (2 * dL)
    w0_prime = (d_p['w0'] - d_m['w0']) / (2 * dL)
    Ct_prime = (d_p['Ct'] - d_m['Ct']) / (2 * dL)
    mu_alpha_prime = (d_p['mu_alpha_formula'] - d_m['mu_alpha_formula']) / (2 * dL)

    print(f"\n  --- L-ABLEITUNGEN (3-Term-Zerlegung) ---")
    print(f"  sigma'      = {sigma_prime:+.10e}")
    print(f"  -w0'        = {-w0_prime:+.10e}")
    print(f"  Ct'         = {Ct_prime:+.10e}")
    print(f"  ---")
    print(f"  sigma'-w0'+Ct' = {sigma_prime - w0_prime + Ct_prime:+.10e}")
    print(f"  d/dL(mu/a)     = {mu_alpha_prime:+.10e}")

    # How big are the individual terms vs. their sum?
    terms = [sigma_prime, -w0_prime, Ct_prime]
    sum_terms = sum(terms)
    max_term = max(abs(t) for t in terms)
    print(f"\n  max(|term|)    = {max_term:+.10e}")
    print(f"  |sum|          = {abs(sum_terms):.10e}")
    print(f"  Cancel-Grad    = {abs(sum_terms)/max_term:.6e}  ({abs(sum_terms)/max_term*100:.4f}%)")

    # Which pairs cancel?
    print(f"\n  --- PAARWEISE CANCELLATION ---")
    print(f"  sigma' + Ct'   = {sigma_prime + Ct_prime:+.10e}  (= w0' + d/dL(mu/a))")
    print(f"  sigma' - w0'   = {sigma_prime - w0_prime:+.10e}  (= -Ct' + d/dL(mu/a))")
    print(f"  -w0' + Ct'     = {-w0_prime + Ct_prime:+.10e}  (= -sigma' + d/dL(mu/a))")

    # Secular identity check
    dL_mH_explicit = 0.0
    u0_base = d_0['w0']
    for bm, hm, bp, hp in zip(d_m['betas'], d_m['hs'], d_p['betas'], d_p['hs']):
        dm = bm**2 / (hm - u0_base) if abs(hm - u0_base) > 1e-30 else 0
        dp = bp**2 / (hp - u0_base) if abs(hp - u0_base) > 1e-30 else 0
        dL_mH_explicit += (dp - dm) / (2 * dL)

    mHp = d_0['m_H_prime']
    secular_LHS = dL_mH_explicit + mHp * w0_prime
    secular_RHS = Ct_prime / d_0['Ct']**2

    print(f"\n  --- SECULAR-IDENTITAET ---")
    print(f"  dL_mH          = {dL_mH_explicit:+.10e}")
    print(f"  mH' * w0'      = {mHp * w0_prime:+.10e}")
    print(f"  LHS = dL_mH + mH'w0' = {secular_LHS:+.10e}")
    print(f"  RHS = Ct'/Ct^2        = {secular_RHS:+.10e}")
    print(f"  Rel. Fehler: {abs(secular_LHS - secular_RHS)/abs(secular_RHS):.4e}")

    # Alignment diagnostics
    alpha_prime = (d_p['alpha'] - d_m['alpha']) / (2 * dL)
    g0k_prime = (d_p['g0k'] - d_m['g0k']) / (2 * dL)
    g0u_prime = (d_p['g0u'] - d_m['g0u']) / (2 * dL)

    print(f"\n  --- ALIGNMENT-DYNAMIK ---")
    print(f"  alpha = <u|k>  = {d_0['alpha']:.8f},  alpha' = {alpha_prime:+.6e}")
    print(f"  <g0|k>         = {d_0['g0k']:.8f},  <g0|k>' = {g0k_prime:+.6e}")
    print(f"  <g0|u>         = {d_0['g0u']:.8f},  <g0|u>' = {g0u_prime:+.6e}")

    # Decomposition of sigma' into H-part, u-part, k-part, alpha-part
    # sigma = <u|H|k>/alpha
    # sigma' = [<u'|H|k> + <u|H'|k> + <u|H|k'>]/alpha - sigma * alpha'/alpha
    # We can't easily compute <u'|H|k> etc. but we can check:
    # sigma' = d/dL(<u|H|k>) / alpha - sigma * alpha'/alpha
    uHk_prime = (d_p['sigma'] * d_p['alpha'] - d_m['sigma'] * d_m['alpha']) / (2 * dL)
    sigma_from_parts = uHk_prime / d_0['alpha'] - d_0['sigma'] * alpha_prime / d_0['alpha']

    print(f"\n  --- SIGMA-ZERLEGUNG ---")
    print(f"  <u|H|k> = sigma*alpha = {d_0['sigma']*d_0['alpha']:+.10e}")
    print(f"  d/dL(<u|H|k>)         = {uHk_prime:+.10e}")
    print(f"  sigma * alpha'/alpha  = {d_0['sigma'] * alpha_prime / d_0['alpha']:+.10e}")
    print(f"  sigma' (direkt)       = {sigma_prime:+.10e}")
    print(f"  sigma' (aus Teilen)   = {sigma_from_parts:+.10e}")

    return {
        'lam': lam_base, 'L': L_base,
        'sigma_prime': sigma_prime, 'w0_prime': w0_prime, 'Ct_prime': Ct_prime,
        'mu_alpha_prime': mu_alpha_prime,
        'dL_mH': dL_mH_explicit, 'mHp': mHp,
        'sigma': d_0['sigma'], 'w0': d_0['w0'], 'Ct': d_0['Ct'],
        'alpha': d_0['alpha'], 'alpha_prime': alpha_prime,
    }


def main():
    print(f"C2bd: Algebraische Zerlegung des Derivative Freezing")
    print(f"DPS={DPS}, N={N_DEFAULT}")

    results = []
    for lam in [3.0, 5.0, 7.0]:
        r = run_decomposition(lam, N_DEFAULT, dlam=0.01)
        results.append(r)

    # Cross-lambda summary
    print(f"\n{'='*90}")
    print(f"  CROSS-LAMBDA ZUSAMMENFASSUNG")
    print(f"{'='*90}")
    print(f"\n  d/dL(mu/a) = sigma' - w0' + Ct'")
    print(f"  {'lam':>5} {'sigma_p':>14} {'-w0_p':>14} {'Ct_p':>14} {'Sum':>14} {'Cancel%':>10}")
    for r in results:
        s = r['sigma_prime'] - r['w0_prime'] + r['Ct_prime']
        mx = max(abs(r['sigma_prime']), abs(r['w0_prime']), abs(r['Ct_prime']))
        cp = (1 - abs(s)/mx) * 100 if mx > 0 else 0
        print(f"  {r['lam']:5.1f} {r['sigma_prime']:+14.6e} {-r['w0_prime']:+14.6e} "
              f"{r['Ct_prime']:+14.6e} {s:+14.6e} {cp:10.2f}%")

    print(f"\n  Welche ZWEI Terme canceln hauptsaechlich?")
    for r in results:
        sp, wp, cp = r['sigma_prime'], r['w0_prime'], r['Ct_prime']
        pairs = [
            ("sigma'+Ct'", sp + cp, "=> w0' dominiert"),
            ("sigma'-w0'", sp - wp, "=> -Ct' dominiert"),
            ("-w0'+Ct'", -wp + cp, "=> -sigma' dominiert"),
        ]
        best = min(pairs, key=lambda x: abs(x[1]))
        print(f"  lam={r['lam']:.0f}: Beste Paar-Cancellation: {best[0]} = {best[1]:+.6e} {best[2]}")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
