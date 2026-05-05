"""
c2bd_derivative_matching.py

Praeziser Test des Derivative-Freezing-Matchings:

  d/dL(mu/alpha) = d/dL(R_bulk) + d/dL(T2_boundary) ≈ 0

und Vergleich mit der Secular-Ableitung:

  partial_L m_H ↔ d/dL(T2_boundary)  (gleiche Richtung?)
  m_H' u_0'     ↔ d/dL(R_bulk)       (gleiche Richtung?)

Verwendet numpy-Diagonalisierung fuer Geschwindigkeit.
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import mpmath as mp
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, chi_trivial,
    mangoldt_mp, rho_even_mp,
)
from c2_approximation_test import (
    h_educated_guess, k_lambda_value, norm,
)
from c2_poisson_decomposition import (
    h_hat_analytical, project_to_fourier,
)

DPS = int(os.environ.get("DPS", 25))
NGRID = 15


def to_np(M, dim):
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = float(M[i, j])
    return A


def K_global(z):
    z_mp = mp.mpf(z)
    if z_mp >= 0:
        ez = mp.exp(z_mp)
        s = mp.mpf(0)
        for n in range(1, 500):
            val = h_educated_guess(mp.mpf(n) * ez)
            s += val
            if n > 3 and abs(val) < mp.mpf("1e-40"):
                break
        return mp.exp(z_mp / 2) * s
    else:
        emz = mp.exp(-z_mp)
        s = mp.mpf(0)
        for m in range(1, 500):
            val = h_hat_analytical(mp.mpf(m) * emz)
            s += val
            if m > 3 and abs(val) < mp.mpf("1e-40"):
                break
        return mp.exp(-z_mp / 2) * s


def compute_full_data(lam, N):
    """Berechne Secular + Bulk/Boundary-Daten."""
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

    # Build operators
    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    Aq_np = to_np(Aq, dim)
    Ah_np = to_np(Ah, dim)
    W02_np = Aq_np - Ah_np

    col0 = W02_np[:, 0]
    cn = np.linalg.norm(col0)
    ut = col0 / cn

    Ctilde = float(ut @ W02_np @ ut)
    uHu = float(ut @ Ah_np @ ut)
    w_bar = Ctilde + uHu

    ws_A, _ = np.linalg.eigh(Aq_np)
    u0 = ws_A[0]

    hs, vs_H = np.linalg.eigh(Ah_np)
    betas = vs_H.T @ ut

    def m_H_at(w):
        return sum(b**2 / (h - w) for b, h in zip(betas, hs) if abs(h - w) > 1e-30)

    def m_H_prime(w):
        return sum(b**2 / (h - w)**2 for b, h in zip(betas, hs) if abs(h - w) > 1e-30)

    mH = m_H_at(u0)
    mHp = m_H_prime(u0)

    # k, alpha, T1, T2
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
    alpha = float(ut @ kn)

    T1 = w_bar - u0
    k_perp = kn - alpha * ut
    Hut = Ah_np @ ut
    f_lam = Hut - uHu * ut
    T2 = float(f_lam @ k_perp) / alpha
    mu_alpha = T1 + T2

    # Residual
    Hk = Ah_np @ kn
    res = Hk - u0 * kn
    res_u = float(ut @ res)

    # B_L on small grid (Fourier projection)
    x_grid = np.linspace(0.02 * L, 0.98 * L, NGRID)
    B_values = np.zeros(NGRID)

    for i, x in enumerate(x_grid):
        x_mp = mp.mpf(x)
        # B'_L (prime boundary)
        ex = mp.exp(x_mp)
        emx = mp.exp(-x_mp)
        s1 = mp.mpf(0)
        for m in range(1, 2000):
            arg = mp.mpf(m) * ex
            hval = h_educated_guess(arg)
            if m > 3 and abs(hval) < mp.mpf("1e-40"):
                break
            threshold = float(mp.exp(L_mp - x_mp))
            ep = mp.mpf(0)
            for k in range(2, m + 1):
                if m % k == 0:
                    lk = mangoldt_mp(k)
                    if lk > 0 and k > threshold:
                        ep += lk
            if ep > 0:
                s1 += ep * hval
        s1 *= mp.exp(x_mp / 2)

        s2 = mp.mpf(0)
        for r in range(1, 2000):
            arg = mp.mpf(r) * emx
            hhval = h_hat_analytical(arg)
            if r > 3 and abs(hhval) < mp.mpf("1e-40"):
                break
            threshold2 = float(mp.exp(x_mp))
            em = mp.mpf(0)
            for k in range(2, r + 1):
                if r % k == 0:
                    lk = mangoldt_mp(k)
                    if lk > 0 and k > threshold2:
                        em += lk
            if em > 0:
                s2 += em * hhval
        s2 *= mp.exp(-x_mp / 2)
        bp = float(s1 + s2)

        # B_R_L (Archimedes boundary)
        Y_MAX = 12.0
        a1 = float(L_mp - x_mp)
        if a1 < 0.01:
            a1 = 0.01
        I1 = mp.mpf(0)
        if a1 < Y_MAX:
            try:
                I1 = mp.quad(lambda y: rho_even_mp(y) * K_global(float(x_mp + y)),
                             [a1, Y_MAX], maxdegree=5)
            except Exception:
                pass
        try:
            I2 = mp.quad(lambda y: rho_even_mp(y) * K_global(float(x_mp - y)),
                         [float(x_mp), Y_MAX], maxdegree=5)
        except Exception:
            I2 = mp.mpf(0)

        B_values[i] = bp + float(I1) + float(I2)

    # Fourier projection of B_L
    dx = x_grid[1] - x_grid[0]
    inv_sqrt_L = 1.0 / np.sqrt(L)
    sqrt_2_L = np.sqrt(2.0 / L)
    b_raw = np.zeros(dim)
    b_raw[0] = np.sum(B_values * inv_sqrt_L) * dx
    for n in range(1, dim):
        basis_n = sqrt_2_L * np.cos(2.0 * np.pi * n * x_grid / L)
        b_raw[n] = np.sum(B_values * basis_n) * dx
    nf_float = float(nf)
    b_norm = b_raw / nf_float

    boundary_u = float(ut @ b_norm)
    ebulk_u = res_u - boundary_u

    T2_boundary = boundary_u / alpha
    T2_bulk = ebulk_u / alpha - uHu + u0
    R_bulk = T1 + T2_bulk

    return {
        'lam': lam, 'L': L,
        'Ctilde': Ctilde, 'uHu': uHu, 'u0': u0, 'w_bar': w_bar,
        'mH': mH, 'mHp': mHp, 'alpha': alpha,
        'T1': T1, 'T2': T2, 'mu_alpha': mu_alpha,
        'T2_bulk': T2_bulk, 'T2_boundary': T2_boundary, 'R_bulk': R_bulk,
        'h_vals': list(hs), 'betas': list(betas),
    }


def run(lam_base, N, dlam=0.01):
    L_base = 2 * np.log(lam_base)
    dL = 2 * dlam / lam_base

    print(f"\n{'='*90}")
    print(f"  DERIVATIVE MATCHING: lam={lam_base}, N={N}, dlam={dlam}")
    print(f"  L={L_base:.6f}, dL={dL:.6e}")
    print(f"{'='*90}")

    t0 = time.time()
    print(f"  lam = {lam_base - dlam:.4f} ...", flush=True)
    d_m = compute_full_data(lam_base - dlam, N)
    print(f"    done ({time.time()-t0:.0f}s)", flush=True)

    t1 = time.time()
    print(f"  lam = {lam_base:.4f} ...", flush=True)
    d_0 = compute_full_data(lam_base, N)
    print(f"    done ({time.time()-t1:.0f}s)", flush=True)

    t2 = time.time()
    print(f"  lam = {lam_base + dlam:.4f} ...", flush=True)
    d_p = compute_full_data(lam_base + dlam, N)
    print(f"    done ({time.time()-t2:.0f}s)", flush=True)

    # Finite differences
    dL_Ct = (d_p['Ctilde'] - d_m['Ctilde']) / (2 * dL)
    dL_u0 = (d_p['u0'] - d_m['u0']) / (2 * dL)
    dL_T1 = (d_p['T1'] - d_m['T1']) / (2 * dL)
    dL_T2 = (d_p['T2'] - d_m['T2']) / (2 * dL)
    dL_mua = (d_p['mu_alpha'] - d_m['mu_alpha']) / (2 * dL)

    dL_T2bd = (d_p['T2_boundary'] - d_m['T2_boundary']) / (2 * dL)
    dL_T2bk = (d_p['T2_bulk'] - d_m['T2_bulk']) / (2 * dL)
    dL_Rbulk = (d_p['R_bulk'] - d_m['R_bulk']) / (2 * dL)

    # Secular at fixed u0_base
    u0_base = d_0['u0']
    mH_m = sum(b**2 / (h - u0_base) for b, h in zip(d_m['betas'], d_m['h_vals'])
               if abs(h - u0_base) > 1e-30)
    mH_p = sum(b**2 / (h - u0_base) for b, h in zip(d_p['betas'], d_p['h_vals'])
               if abs(h - u0_base) > 1e-30)
    dL_mH = (mH_p - mH_m) / (2 * dL)
    mHp = d_0['mHp']
    Ct = d_0['Ctilde']

    print(f"\n  --- Werte bei lam={lam_base} ---")
    print(f"  T1 = {d_0['T1']:+.6e}, T2 = {d_0['T2']:+.6e}")
    print(f"  mu/a = {d_0['mu_alpha']:+.6e}")
    print(f"  T2_bd = {d_0['T2_boundary']:+.6e}, R_bulk = {d_0['R_bulk']:+.6e}")

    print(f"\n  --- L-Ableitungen ---")
    print(f"  d/dL(T1)      = {dL_T1:+.6e}")
    print(f"  d/dL(T2)      = {dL_T2:+.6e}")
    print(f"  d/dL(mu/a)    = {dL_mua:+.6e}")
    print(f"  d/dL(T2_bd)   = {dL_T2bd:+.6e}")
    print(f"  d/dL(R_bulk)  = {dL_Rbulk:+.6e}")

    cancel = abs(dL_mua) / (abs(dL_T2bd) + abs(dL_Rbulk)) if (abs(dL_T2bd) + abs(dL_Rbulk)) > 0 else float('nan')
    print(f"\n  --- DERIVATIVE FREEZING ---")
    print(f"  d/dL(T2_bd) + d/dL(R_bulk) = {dL_T2bd + dL_Rbulk:+.6e}")
    print(f"  d/dL(mu/a)                 = {dL_mua:+.6e}")
    print(f"  Cancel-Ratio: {cancel:.4f}  ({cancel*100:.1f}%)")
    print(f"  => {(1-cancel)*100:.1f}% Cancellation in L-Ableitungen")

    print(f"\n  --- SECULAR-ABLEITUNG ---")
    print(f"  dL_mH    = {dL_mH:+.6e}  (Boundary)")
    print(f"  mHp*u0p  = {mHp * dL_u0:+.6e}  (Rayleigh)")
    lhs = dL_mH + mHp * dL_u0
    rhs = dL_Ct / Ct**2
    print(f"  LHS      = {lhs:+.6e}")
    print(f"  RHS      = {rhs:+.6e}")

    print(f"\n  --- MATCHING ---")
    if abs(dL_mH) > 1e-30:
        print(f"  d/dL(T2_bd)  / dL_mH   = {dL_T2bd / dL_mH:.6e}")
    if abs(mHp * dL_u0) > 1e-30:
        print(f"  d/dL(R_bulk) / mHp*u0p = {dL_Rbulk / (mHp * dL_u0):.6e}")
    if abs(rhs) > 1e-30:
        print(f"  d/dL(mu/a)   / Cp/C^2  = {dL_mua / rhs:.6e}")

    # Sign matching
    sign_bd = "+" if dL_T2bd > 0 else "-"
    sign_mH = "+" if dL_mH > 0 else "-"
    sign_Rb = "+" if dL_Rbulk > 0 else "-"
    sign_Rp = "+" if mHp * dL_u0 > 0 else "-"
    print(f"\n  --- VORZEICHEN ---")
    print(f"  d/dL(T2_bd) = {sign_bd},  dL_mH   = {sign_mH}  {'MATCH' if sign_bd == sign_mH else 'MISMATCH'}")
    print(f"  d/dL(R_bulk)= {sign_Rb},  mHp*u0p = {sign_Rp}  {'MATCH' if sign_Rb == sign_Rp else 'MISMATCH'}")

    return {
        'lam': lam_base, 'L': L_base, 'dL': dL,
        'dL_mua': dL_mua,
        'dL_T2bd': dL_T2bd, 'dL_Rbulk': dL_Rbulk,
        'dL_mH': dL_mH, 'mHp_u0p': mHp * dL_u0,
        'Ct_prime_Ct2': rhs,
        'Ct': Ct, 'mHp': mHp, 'dL_u0': dL_u0, 'dL_Ct': dL_Ct,
        'alpha': d_0['alpha'], 'mu_alpha': d_0['mu_alpha'],
        'T2_bd': d_0['T2_boundary'], 'R_bulk': d_0['R_bulk'],
    }


def multi_lambda_analysis(results):
    print(f"\n{'='*90}")
    print(f"  MULTI-LAMBDA KOEFFIZIENTEN-ANALYSE")
    print(f"{'='*90}")

    n = len(results)
    print(f"\n  {n} Datenpunkte: lam = {[r['lam'] for r in results]}")

    print(f"\n  --- Rohdaten ---")
    print(f"  {'lam':>5} {'d/dL(mu/a)':>14} {'dL_mH':>14} {'mHp*u0p':>14} {'Ct_p/Ct2':>14} {'alpha':>10}")
    for r in results:
        print(f"  {r['lam']:5.1f} {r['dL_mua']:+14.6e} {r['dL_mH']:+14.6e} "
              f"{r['mHp_u0p']:+14.6e} {r['Ct_prime_Ct2']:+14.6e} {r['alpha']:10.6f}")

    # Fit: d/dL(mu/a) = a * dL_mH + b * mHp_u0p + c * Ct'/Ct^2
    if n >= 2:
        y = np.array([r['dL_mua'] for r in results])
        X = np.array([[r['dL_mH'], r['mHp_u0p'], r['Ct_prime_Ct2']] for r in results])

        print(f"\n  --- Least-Squares Fit ---")
        print(f"  d/dL(mu/a) = a * dL_mH + b * mHp*u0p + c * Ct'/Ct^2")

        if n >= 3:
            coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
            a, b, c = coeffs
            fitted = X @ coeffs
            print(f"  a = {a:+.8e}")
            print(f"  b = {b:+.8e}")
            print(f"  c = {c:+.8e}")
            print(f"  Residuen: {[f'{v:+.2e}' for v in (y - fitted)]}")
        else:
            print(f"  (nur {n} Punkte, teste 2-Parameter-Fits)")

        # 2-param fits: alle 3 Kombinationen
        labels_2 = [
            ("a * dL_mH + b * mHp*u0p", [0, 1]),
            ("a * dL_mH + c * Ct'/Ct^2", [0, 2]),
            ("b * mHp*u0p + c * Ct'/Ct^2", [1, 2]),
        ]
        for label, idx in labels_2:
            X2 = X[:, idx]
            c2, res2, _, _ = np.linalg.lstsq(X2, y, rcond=None)
            fitted2 = X2 @ c2
            rel_err = np.max(np.abs(y - fitted2)) / np.max(np.abs(y)) if np.max(np.abs(y)) > 0 else float('nan')
            names = ["a", "b", "c"]
            coeff_str = ", ".join(f"{names[i]}={c2[j]:+.8e}" for j, i in enumerate(idx))
            print(f"\n  {label}")
            print(f"    {coeff_str}")
            print(f"    Max rel. Fehler: {rel_err:.2e}")
            print(f"    Residuen: {[f'{v:+.2e}' for v in (y - fitted2)]}")

    # Direkte Ratios
    print(f"\n  --- Direkte Ratios ---")
    print(f"  {'lam':>5} {'mu_a/dL_mH':>14} {'mu_a/mHp_u0p':>14} {'mu_a/Ct_p':>14} {'T2bd/dL_mH':>14} {'Rbulk/mHp_u0p':>14}")
    for r in results:
        r1 = r['dL_mua'] / r['dL_mH'] if abs(r['dL_mH']) > 1e-40 else float('nan')
        r2 = r['dL_mua'] / r['mHp_u0p'] if abs(r['mHp_u0p']) > 1e-40 else float('nan')
        r3 = r['dL_mua'] / r['Ct_prime_Ct2'] if abs(r['Ct_prime_Ct2']) > 1e-40 else float('nan')
        r4 = r['dL_T2bd'] / r['dL_mH'] if abs(r['dL_mH']) > 1e-40 else float('nan')
        r5 = r['dL_Rbulk'] / r['mHp_u0p'] if abs(r['mHp_u0p']) > 1e-40 else float('nan')
        print(f"  {r['lam']:5.1f} {r1:+14.6e} {r2:+14.6e} {r3:+14.6e} {r4:+14.6e} {r5:+14.6e}")

    # Skalierung mit alpha
    print(f"\n  --- Skalierung mit alpha ---")
    print(f"  {'lam':>5} {'alpha':>10} {'alpha^2':>12} {'mu_a/alpha':>14} {'mu_a/alpha^2':>14} {'dL_T2bd/alpha':>14}")
    for r in results:
        a = r['alpha']
        print(f"  {r['lam']:5.1f} {a:10.6f} {a**2:12.8f} "
              f"{r['dL_mua']/a:+14.6e} {r['dL_mua']/a**2:+14.6e} "
              f"{r['dL_T2bd']/a:+14.6e}")


def main():
    print(f"C2bd: Derivative Matching (Freezing Test) — Multi-Lambda")
    print(f"DPS={DPS}, NGRID={NGRID}")

    results = []
    for lam in [3.0, 4.0, 5.0, 7.0, 9.0]:
        N = 15
        dlam = 0.01
        r = run(lam, N, dlam=dlam)
        results.append(r)

    multi_lambda_analysis(results)


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
