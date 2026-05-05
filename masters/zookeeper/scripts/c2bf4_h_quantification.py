"""
c2bf4_h_quantification.py

Berechnet h = P_perp (A - u0 I) k und seine L-Ableitung
in der Notation von C2BI (Exact Poisson Branch Transfer).

P_perp = I - |u><u| (Projektor senkrecht zu u-tilde)
T = P_perp H P_perp (H eingeschraenkt auf Komplement von u)
R0 = (T - u0 I)^{-1} (Resolvent im Komplement)
f = P_perp H |u> = H|u> - <u|H|u> |u>
h = P_perp (A - u0 I) k = P_perp (H - u0 I) k  [da P_perp |u> = 0]

mu/alpha = (1/alpha) <f, R0 h>

Dieses Script quantifiziert ||h||, ||R0 h||, <f, R0 h> und deren L-Ableitungen.
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


def compute_h_data(lam, N):
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
    ut = col0 / cn  # u-tilde

    Ct = float(ut @ W_np @ ut)
    uHu = float(ut @ H_np @ ut)

    ws_A, vs_A = np.linalg.eigh(A_np)
    w0 = ws_A[0]

    hs, vs_H = np.linalg.eigh(H_np)

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
    alpha = float(ut @ kn)

    # P_perp = I - |ut><ut|
    P_perp = np.eye(dim) - np.outer(ut, ut)

    # f = P_perp H |ut> = H|ut> - <ut|H|ut> |ut>
    f_vec = H_np @ ut - uHu * ut

    # h = P_perp (H - w0 I) k  [= P_perp (A - w0 I) k since P_perp |u>=0]
    h_vec = P_perp @ (H_np - w0 * np.eye(dim)) @ kn

    # T = P_perp H P_perp (restricted to complement of ut)
    T_np = P_perp @ H_np @ P_perp

    # R0 = (T - w0 I)^{-1} in the complement space
    # Need pseudo-inverse since T has zero eigenvalue in ut-direction
    T_shifted = T_np - w0 * P_perp  # = P_perp (H - w0) P_perp
    # Regularize: add small term in ut-direction
    T_reg = T_shifted + 1e6 * np.outer(ut, ut)
    R0 = np.linalg.inv(T_reg)
    R0 = P_perp @ R0 @ P_perp  # Project back to complement

    R0h = R0 @ h_vec
    fR0h = float(f_vec @ R0h)
    mu_alpha_resolvent = fR0h / alpha

    # Direct check
    mu_alpha_direct = float(ut @ (A_np - w0 * np.eye(dim)) @ kn) / alpha

    # Norms
    h_norm = np.linalg.norm(h_vec)
    R0h_norm = np.linalg.norm(R0h)
    f_norm = np.linalg.norm(f_vec)

    # Spectral decomposition of h in H-eigenbasis
    h_coeffs = vs_H.T @ h_vec
    f_coeffs = vs_H.T @ f_vec

    # The pairing <f, R0 h> in H-eigenbasis
    # R0 in H-eigenbasis: diagonal with entries 1/(h_i - w0) for modes perp to ut
    # Actually T is not diagonal in H-eigenbasis, but the H-eigendecomposition gives insight

    return {
        'lam': lam, 'L': L, 'w0': w0, 'Ct': Ct, 'uHu': uHu, 'alpha': alpha,
        'h_norm': h_norm, 'R0h_norm': R0h_norm, 'f_norm': f_norm,
        'fR0h': fR0h, 'mu_alpha_resolvent': mu_alpha_resolvent,
        'mu_alpha_direct': mu_alpha_direct,
        'h_coeffs': h_coeffs, 'f_coeffs': f_coeffs,
        'hs': hs,
        'h_vec': h_vec, 'f_vec': f_vec, 'R0h': R0h,
        'ut': ut, 'kn': kn,
    }


def main():
    print(f"C2bf.4: Quantifizierung von h = P_perp(A-u0)k")
    print(f"DPS={DPS}, N={N_DEFAULT}")

    for lam in [3.0, 5.0, 7.0]:
        L = 2 * np.log(lam)
        dL = 2 * 0.01 / lam

        print(f"\n{'='*80}")
        print(f"  lam={lam}, L={L:.6f}")
        print(f"{'='*80}")

        t0 = time.time()
        d_m = compute_h_data(lam - 0.01, N_DEFAULT)
        d_0 = compute_h_data(lam, N_DEFAULT)
        d_p = compute_h_data(lam + 0.01, N_DEFAULT)
        print(f"  3 Punkte in {time.time()-t0:.0f}s", flush=True)

        print(f"\n  --- VERIFIKATION ---")
        print(f"  mu/a (resolvent) = {d_0['mu_alpha_resolvent']:+.10e}")
        print(f"  mu/a (direkt)    = {d_0['mu_alpha_direct']:+.10e}")

        print(f"\n  --- GROESSEN ---")
        print(f"  ||h||      = {d_0['h_norm']:.6e}")
        print(f"  ||R0 h||   = {d_0['R0h_norm']:.6e}")
        print(f"  ||f||      = {d_0['f_norm']:.6e}")
        print(f"  <f,R0 h>   = {d_0['fR0h']:+.6e}")
        print(f"  alpha      = {d_0['alpha']:.6f}")
        print(f"  mu/a       = <f,R0 h>/a = {d_0['fR0h']/d_0['alpha']:+.6e}")

        # Cauchy-Schwarz bound
        cs_bound = d_0['f_norm'] * d_0['R0h_norm']
        print(f"\n  --- CANCELLATION in <f, R0 h> ---")
        print(f"  |<f,R0 h>| = {abs(d_0['fR0h']):.6e}")
        print(f"  ||f||*||R0h|| = {cs_bound:.6e}")
        print(f"  Cancellation: {abs(d_0['fR0h'])/cs_bound:.6e}  ({abs(d_0['fR0h'])/cs_bound*100:.4f}%)")

        # L-derivatives of norms
        h_norm_prime = (d_p['h_norm'] - d_m['h_norm']) / (2 * dL)
        R0h_norm_prime = (d_p['R0h_norm'] - d_m['R0h_norm']) / (2 * dL)
        fR0h_prime = (d_p['fR0h'] - d_m['fR0h']) / (2 * dL)
        alpha_prime = (d_p['alpha'] - d_m['alpha']) / (2 * dL)

        print(f"\n  --- L-ABLEITUNGEN ---")
        print(f"  d/dL ||h||      = {h_norm_prime:+.6e}")
        print(f"  d/dL ||R0h||    = {R0h_norm_prime:+.6e}")
        print(f"  d/dL <f,R0h>    = {fR0h_prime:+.6e}")
        print(f"  d/dL alpha      = {alpha_prime:+.6e}")

        # Decompose d/dL(mu/a) = d/dL(<f,R0h>/alpha)
        # = [d/dL(<f,R0h>)]/alpha - (mu/a) * alpha'/alpha
        mu_a = d_0['fR0h'] / d_0['alpha']
        dmu_part1 = fR0h_prime / d_0['alpha']
        dmu_part2 = -mu_a * alpha_prime / d_0['alpha']
        dmu_total = (d_p['mu_alpha_direct'] - d_m['mu_alpha_direct']) / (2 * dL)

        print(f"\n  --- RESOLVENT-ZERLEGUNG von d/dL(mu/a) ---")
        print(f"  [d/dL <f,R0h>]/alpha  = {dmu_part1:+.6e}")
        print(f"  -(mu/a)(alpha'/alpha) = {dmu_part2:+.6e}")
        print(f"  Summe                 = {dmu_part1 + dmu_part2:+.6e}")
        print(f"  d/dL(mu/a) direkt     = {dmu_total:+.6e}")

        # How small is h in spectral terms?
        # h in H-eigenbasis: h_j = Σ_i <phi_i|h> delta_ij
        print(f"\n  --- h im H-Eigensystem (groesste 5 Komponenten) ---")
        h_abs = np.abs(d_0['h_coeffs'])
        top5 = np.argsort(h_abs)[-5:][::-1]
        for idx in top5:
            gap = d_0['hs'][idx] - d_0['w0']
            print(f"    Mode {idx:3d}: |h_j| = {h_abs[idx]:.4e}, "
                  f"h_j - w0 = {gap:+.4e}, |h_j|/|gap| = {h_abs[idx]/abs(gap):.4e}")

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
