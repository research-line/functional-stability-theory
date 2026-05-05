"""
c2bd_secular_derivative.py

Numerische Verifikation der Secular-Ableitung (Route C aus C2BB):

  F(L, w) = 1 + C~(L) * m_H(L; w) = 0

Differenziert entlang w = u_0(L):

  dL_mH + mH' * u0' = C~' / C~^2

wobei:
  dL_mH = partial_L m_H(L; u_0)     -- Boundary-Beitrag
  mH'   = partial_w m_H(L; u_0)     -- Slope der Stieltjes-Trafo
  u0'   = d u_0 / dL                -- Eigenwertverschub
  C~'   = d C~ / dL                 -- Rang-1-Kopplungsaenderung

Ziel: Zeigen dass dL_mH ~ T2_boundary und mH'*u0' ~ R_bulk^lambda,
      und die Secular-Gleichung deren Cancellation erzwingt.
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
    h_educated_guess, k_lambda_value, inner_product, norm,
)
from c2_poisson_decomposition import project_to_fourier

DPS = int(os.environ.get("DPS", 30))


def to_np(M, dim):
    """mpmath matrix -> numpy array."""
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = float(M[i, j])
    return A


def compute_secular_data(lam, N):
    """Berechne alle Secular-Groessen bei gegebenem lambda.
    Nutzt numpy fuer Diagonalisierung (schnell), mpmath nur fuer Matrixbau."""
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    Aq_np = to_np(Aq, dim)
    Ah_np = to_np(Ah, dim)
    W02_np = Aq_np - Ah_np

    col0 = W02_np[:, 0]
    cn = np.linalg.norm(col0)
    ut_np = col0 / cn

    Ctilde = ut_np @ W02_np @ ut_np
    uHu = ut_np @ Ah_np @ ut_np
    w_bar = Ctilde + uHu

    ws_A, vs_A = np.linalg.eigh(Aq_np)
    u0 = ws_A[0]

    hs, vs_H = np.linalg.eigh(Ah_np)
    betas = vs_H.T @ ut_np

    def m_H_at(w):
        s = 0.0
        for b, h in zip(betas, hs):
            denom = h - w
            if abs(denom) < 1e-30:
                continue
            s += b**2 / denom
        return s

    def m_H_prime(w):
        s = 0.0
        for b, h in zip(betas, hs):
            denom = h - w
            if abs(denom) < 1e-30:
                continue
            s += b**2 / denom**2
        return s

    mH_u0 = m_H_at(u0)
    mHp_u0 = m_H_prime(u0)
    secular = 1 + Ctilde * mH_u0

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf_mp = norm(cf, dim)
    kn_np = np.array([float(cf[i, 0] / nf_mp) for i in range(dim)])
    alpha = float(ut_np @ kn_np)

    T1 = w_bar - u0

    k_perp_np = kn_np - alpha * ut_np
    Hut_np = Ah_np @ ut_np
    f_lam_np = Hut_np - uHu * ut_np

    T2 = float(f_lam_np @ k_perp_np) / alpha
    mu_alpha = T1 + T2

    return {
        'lam': lam, 'L': L, 'dim': dim,
        'Ctilde': Ctilde, 'uHu': uHu, 'u0': u0, 'w_bar': w_bar,
        'mH_u0': mH_u0, 'mHp_u0': mHp_u0, 'secular': secular,
        'alpha': alpha, 'T1': T1, 'T2': T2, 'mu_alpha': mu_alpha,
        'h_vals': list(hs), 'betas': list(betas),
    }


def run(lam_base, N, dlam=0.005):
    """Secular-Ableitung bei lambda_base mit Schrittweite dlam."""
    mp.mp.dps = DPS
    L_base = 2 * np.log(lam_base)
    dL = 2 * dlam / lam_base

    print(f"\n{'='*90}")
    print(f"  C2bd SECULAR DERIVATIVE: lam={lam_base}, N={N}, dlam={dlam}")
    print(f"  L_base = {L_base:.6f}, dL = {dL:.6e}")
    print(f"{'='*90}")

    t0 = time.time()

    print(f"\n  Computing at lam = {lam_base - dlam:.4f} ...", flush=True)
    d_minus = compute_secular_data(lam_base - dlam, N)
    print(f"  Computing at lam = {lam_base:.4f} ...", flush=True)
    d_base = compute_secular_data(lam_base, N)
    print(f"  Computing at lam = {lam_base + dlam:.4f} ...", flush=True)
    d_plus = compute_secular_data(lam_base + dlam, N)

    print(f"  All three points computed in {time.time()-t0:.0f}s", flush=True)

    # --- Secular verification ---
    print(f"\n  --- Secular-Gleichung ---")
    for label, d in [("lam-", d_minus), ("lam0", d_base), ("lam+", d_plus)]:
        print(f"  {label}: C~={d['Ctilde']:.6f}, mH={d['mH_u0']:.6e}, "
              f"1+C~mH={d['secular']:.2e}, u0={d['u0']:.8f}")

    # --- Finite differences ---
    dL_Ctilde = (d_plus['Ctilde'] - d_minus['Ctilde']) / (2 * dL)
    dL_u0 = (d_plus['u0'] - d_minus['u0']) / (2 * dL)

    # partial_L m_H at FIXED w = u0_base
    u0_base = d_base['u0']
    mH_minus_at_u0base = sum(
        b**2 / (h - u0_base) for b, h in zip(d_minus['betas'], d_minus['h_vals'])
        if abs(h - u0_base) > 1e-30)
    mH_plus_at_u0base = sum(
        b**2 / (h - u0_base) for b, h in zip(d_plus['betas'], d_plus['h_vals'])
        if abs(h - u0_base) > 1e-30)
    dL_mH = (mH_plus_at_u0base - mH_minus_at_u0base) / (2 * dL)

    mHp = d_base['mHp_u0']
    Ct = d_base['Ctilde']

    print(f"\n  --- Finite Differences (dL = {dL:.6e}) ---")
    print(f"  C~'(L)       = {dL_Ctilde:+.10e}")
    print(f"  u_0'(L)      = {dL_u0:+.10e}")
    print(f"  dL_mH        = {dL_mH:+.10e}  (partial_L m_H at fixed u0)")
    print(f"  mH'(u0)      = {mHp:+.10e}  (partial_w m_H)")
    print(f"  mH'*u0'      = {mHp * dL_u0:+.10e}")

    # --- Verify secular derivative identity ---
    lhs = dL_mH + mHp * dL_u0
    rhs = dL_Ctilde / Ct**2

    print(f"\n  --- Secular-Ableitung: dL_mH + mH'*u0' = C~'/C~^2 ---")
    print(f"  LHS = dL_mH + mH'*u0' = {lhs:+.10e}")
    print(f"  RHS = C~'/C~^2         = {rhs:+.10e}")
    print(f"  Differenz              = {abs(lhs - rhs):.2e}")
    print(f"  Rel. Fehler            = {abs(lhs - rhs) / abs(rhs):.2e}")

    # --- Decomposition in Bulk vs Boundary ---
    print(f"\n  --- Zerlegung in Boundary vs Rayleigh ---")
    print(f"  dL_mH    = {dL_mH:+.10e}  (BOUNDARY: wie H-Spektrum sich mit L aendert)")
    print(f"  mH'*u0'  = {mHp * dL_u0:+.10e}  (RAYLEIGH: wie u0 sich mit L verschiebt)")
    print(f"  Summe    = {lhs:+.10e}  (= C~'/C~^2, erzwungen durch Secular)")

    ratio_boundary = dL_mH / lhs if abs(lhs) > 1e-30 else float('nan')
    ratio_rayleigh = mHp * dL_u0 / lhs if abs(lhs) > 1e-30 else float('nan')
    print(f"  Anteil Boundary: {ratio_boundary:.4f}  ({ratio_boundary*100:.1f}%)")
    print(f"  Anteil Rayleigh: {ratio_rayleigh:.4f}  ({ratio_rayleigh*100:.1f}%)")

    # --- Connection to C2bb/C2bc ---
    print(f"\n  --- Verbindung zu C2bb/C2bc ---")
    print(f"  T1 = w_bar - u0     = {d_base['T1']:+.10e}")
    print(f"  T2                   = {d_base['T2']:+.10e}")
    print(f"  mu/alpha             = {d_base['mu_alpha']:+.10e}")
    print(f"  alpha                = {d_base['alpha']:.6f}")

    # Relate dL_mH to T2_boundary and mH'*u0' to R_bulk
    # mu/alpha = R_bulk + T2_boundary
    # Hypothesis: R_bulk ~ const * mH'*u0' and T2_boundary ~ const * dL_mH

    print(f"\n  --- Skalierungsvergleich ---")
    print(f"  dL_mH        = {dL_mH:+.10e}")
    print(f"  mH'*u0'      = {mHp * dL_u0:+.10e}")
    print(f"  Vorzeichen: dL_mH {'>' if dL_mH > 0 else '<'} 0, "
          f"mH'*u0' {'>' if mHp*dL_u0 > 0 else '<'} 0")

    if dL_mH * mHp * dL_u0 < 0:
        print(f"  *** GEGENLAEUIG — Cancellation durch Secular erzwungen! ***")
    else:
        print(f"  WARNUNG: Gleichsinnig — keine Cancellation")

    # Effective cancellation ratio
    if abs(dL_mH) + abs(mHp * dL_u0) > 0:
        cancel = abs(lhs) / (abs(dL_mH) + abs(mHp * dL_u0))
        print(f"  Cancel-Ratio = {cancel:.6e}  "
              f"(1 = keine Cancel., 0 = perfekte Cancel.)")

    return {
        'lam': lam_base, 'L': L_base, 'dL': dL,
        'Ctilde': Ct, 'dL_Ctilde': dL_Ctilde,
        'u0': u0_base, 'dL_u0': dL_u0,
        'mH': d_base['mH_u0'], 'mHp': mHp,
        'dL_mH': dL_mH, 'mHp_u0p': mHp * dL_u0,
        'lhs': lhs, 'rhs': rhs,
        'T1': d_base['T1'], 'T2': d_base['T2'],
        'mu_alpha': d_base['mu_alpha'],
    }


def main():
    print(f"C2bd: Secular Derivative Verification (Route C)", flush=True)
    print(f"DPS = {DPS}", flush=True)

    results = []
    for cfg in [
        {"lam": 3.0, "N": 30, "dlam": 0.005},
        {"lam": 5.0, "N": 55, "dlam": 0.005},
    ]:
        r = run(cfg['lam'], cfg['N'], cfg['dlam'])
        results.append(r)

    print(f"\n{'='*90}")
    print(f"CROSS-LAMBDA ZUSAMMENFASSUNG")
    print(f"{'='*90}")

    print(f"\n  {'Groesse':<25} {'lam=3':>15} {'lam=5':>15} {'Ratio':>10}")
    print(f"  {'-'*65}")
    for key in ['dL_mH', 'mHp_u0p', 'lhs', 'rhs', 'mu_alpha']:
        v3 = results[0][key]
        v5 = results[1][key]
        ratio = v5 / v3 if abs(v3) > 1e-30 else float('nan')
        print(f"  {key:<25} {v3:+15.6e} {v5:+15.6e} {ratio:10.3f}")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
