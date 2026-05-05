"""
c2bs_cancellation_identity.py

Tests the projected frontier cancellation identity:
  <R0f, h_bulk> ≈ -<R0f, h_bd>

Defines:
  h_bulk = P_perp (H-w0) k_E   (bulk contribution to residual)
  h_bd   = P_perp (H-w0) k_B   (boundary contribution to residual)
  R0f    = sum_j f_j/(t_j-w0) v_j  (resolvent of T at w0, applied to f)
  f      = P_perp H ut           (coupling vector)

Output:
  S_bulk, S_bd, rho_L = S_bulk + S_bd
  mu/alpha (direct)
  cancel_ratio = |rho| / (|S_bulk| + |S_bd|)
  sign_ratio = S_bulk / S_bd
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
NGRID = 40


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


def compute_boundary_fourier(lam, N, L, L_mp, dim, nf_float):
    """Compute B_L Fourier coefficients (normalized by ||Pi K||)."""
    lam_mp = mp.mpf(lam)
    L_val = float(L_mp)
    x_grid = np.linspace(0.02 * L_val, 0.98 * L_val, NGRID)
    B_values = np.zeros(NGRID)

    for i, x in enumerate(x_grid):
        x_mp = mp.mpf(x)
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

    dx = x_grid[1] - x_grid[0]
    inv_sqrt_L = 1.0 / np.sqrt(L_val)
    sqrt_2_L = np.sqrt(2.0 / L_val)
    b_raw = np.zeros(dim)
    b_raw[0] = np.sum(B_values * inv_sqrt_L) * dx
    for n in range(1, dim):
        basis_n = sqrt_2_L * np.cos(2.0 * np.pi * n * x_grid / L_val)
        b_raw[n] = np.sum(B_values * basis_n) * dx

    return b_raw / nf_float


def compute_cancellation(lam, N):
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

    A_np = to_np(Aq, dim)
    H_np = to_np(Ah, dim)
    W_np = A_np - H_np

    # ut (rank-1 direction)
    col0 = W_np[:, 0]
    cn = np.linalg.norm(col0)
    ut = col0 / cn

    # Eigenvalues of A
    ws_A, _ = np.linalg.eigh(A_np)
    w0 = ws_A[0]

    P_perp = np.eye(dim) - np.outer(ut, ut)

    # T = P_perp H P_perp and its eigensystem
    T_np = P_perp @ H_np @ P_perp
    ts_all, vs_T_all = np.linalg.eigh(T_np)

    # Non-trivial T-eigenvalues
    t_mask = np.abs(ts_all) > 1e-10
    t_idx = np.where(t_mask)[0]
    ts = ts_all[t_idx]
    vs_T = vs_T_all[:, t_idx]
    t_sort = np.argsort(ts)
    ts = ts[t_sort]
    vs_T = vs_T[:, t_sort]
    n_T = len(ts)

    # f = P_perp H ut (coupling vector)
    f_vec = P_perp @ (H_np @ ut)

    # f_j = <v_j, f>
    f_j = np.array([float(vs_T[:, j] @ f_vec) for j in range(n_T)])

    # R0f = sum_j f_j/(t_j - w0) v_j  (resolvent in T-eigenbasis)
    R0f = np.zeros(dim)
    for j in range(n_T):
        if abs(ts[j] - w0) > 1e-14:
            R0f += (f_j[j] / (ts[j] - w0)) * vs_T[:, j]

    # k vector (Poisson)
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    nf_float = float(nf)
    kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
    alpha = float(ut @ kn)

    # B_L Fourier coefficients (normalized)
    k_B = compute_boundary_fourier(lam, N, L, L_mp, dim, nf_float)

    # Bulk = K - B
    k_E = kn - k_B

    # Residuals from each part: (H - w0) applied
    res_E = (H_np - w0 * np.eye(dim)) @ k_E
    res_B = (H_np - w0 * np.eye(dim)) @ k_B
    res_full = (H_np - w0 * np.eye(dim)) @ kn

    # Project to ut-perp
    h_bulk = P_perp @ res_E
    h_bd = P_perp @ res_B
    h_full = P_perp @ res_full

    # mu/alpha direct
    mu_alpha = float(ut @ (A_np - w0 * np.eye(dim)) @ kn) / alpha

    # R0f pairings
    pairing_bulk = float(R0f @ h_bulk)
    pairing_bd = float(R0f @ h_bd)
    pairing_full = float(R0f @ h_full)

    # Normalized: S = (1/alpha) * pairing
    S_bulk = pairing_bulk / alpha
    S_bd = pairing_bd / alpha
    rho_L = S_bulk + S_bd

    # Also compute mu-components along ut
    mu_E = float(ut @ res_E)
    mu_B = float(ut @ res_B)
    mu_full = float(ut @ res_full)

    # C_tilde contribution
    C_tilde = cn
    mu_from_A = mu_full + C_tilde * alpha

    return {
        'lam': lam, 'N': N, 'dim': dim, 'L': L,
        'w0': w0, 'alpha': alpha, 'C_tilde': C_tilde,
        'mu_alpha': mu_alpha,
        # R0f pairings
        'S_bulk': S_bulk,
        'S_bd': S_bd,
        'rho_L': rho_L,
        'pairing_full': pairing_full / alpha,
        # mu-channel
        'mu_E_over_alpha': mu_E / alpha,
        'mu_B_over_alpha': mu_B / alpha,
        'mu_sum_over_alpha': (mu_E + mu_B) / alpha,
        # Vector norms
        'h_bulk_norm': np.linalg.norm(h_bulk),
        'h_bd_norm': np.linalg.norm(h_bd),
        'h_full_norm': np.linalg.norm(h_full),
        'k_E_norm': np.linalg.norm(k_E),
        'k_B_norm': np.linalg.norm(k_B),
        # Cancellation metrics
        'cancel_ratio_R0f': abs(rho_L) / (abs(S_bulk) + abs(S_bd)) if (abs(S_bulk) + abs(S_bd)) > 1e-30 else 0,
        'cancel_ratio_mu': abs(mu_E + mu_B) / (abs(mu_E) + abs(mu_B)) if (abs(mu_E) + abs(mu_B)) > 1e-30 else 0,
        'cancel_ratio_h': np.linalg.norm(h_full) / (np.linalg.norm(h_bulk) + np.linalg.norm(h_bd)),
        'sign_ratio': S_bulk / S_bd if abs(S_bd) > 1e-30 else float('inf'),
        'sign_ratio_mu': mu_E / mu_B if abs(mu_B) > 1e-30 else float('inf'),
        # R0f norm for reference
        'R0f_norm': np.linalg.norm(R0f),
        'n_T': n_T,
    }


def main():
    print("C2bs: Cancellation Identity Test")
    print(f"DPS={DPS}, NGRID={NGRID}")
    print()
    print("Teste: <R0f, h_bulk> ≈ -<R0f, h_bd>")
    print("       mu_E/alpha ≈ -mu_B/alpha")
    print()

    results = []
    for lam in [3.0, 5.0, 7.0, 9.0]:
        for N in [15, 20]:
            t0 = time.time()
            d = compute_cancellation(lam, N)
            elapsed = time.time() - t0
            results.append(d)

            print(f"{'='*80}")
            print(f"  lam={lam}, N={N} ({elapsed:.0f}s)")
            print(f"{'='*80}")
            print(f"  w0={d['w0']:.8f}, alpha={d['alpha']:.6f}, C_tilde={d['C_tilde']:.6f}")
            print(f"  ||R0f|| = {d['R0f_norm']:.6e}, n_T = {d['n_T']}")
            print()

            print(f"  --- R0f-PAIRING (Resolvent-Kanal) ---")
            print(f"  S_bulk  = {d['S_bulk']:+.10e}")
            print(f"  S_bd    = {d['S_bd']:+.10e}")
            print(f"  rho_L   = S_bulk + S_bd = {d['rho_L']:+.10e}")
            print(f"  sign(S_bulk/S_bd) = {d['sign_ratio']:+.6f}")
            cr = d['cancel_ratio_R0f']
            print(f"  cancel_ratio = |rho|/(|S_bulk|+|S_bd|) = {cr:.6e}")
            print(f"  => Cancellation: {1/cr:.1f}x  ({'STARK' if cr < 0.01 else 'SCHWACH' if cr > 0.1 else 'MITTEL'})")
            print()

            print(f"  --- mu-KANAL (ut-Projektion) ---")
            print(f"  mu_E/alpha = {d['mu_E_over_alpha']:+.10e}")
            print(f"  mu_B/alpha = {d['mu_B_over_alpha']:+.10e}")
            print(f"  Summe      = {d['mu_sum_over_alpha']:+.10e}")
            print(f"  mu/alpha (direkt, via A) = {d['mu_alpha']:+.10e}")
            print(f"  sign(mu_E/mu_B) = {d['sign_ratio_mu']:+.6f}")
            cr_mu = d['cancel_ratio_mu']
            print(f"  cancel_ratio_mu = {cr_mu:.6e}")
            print()

            print(f"  --- VEKTOR-NORMEN ---")
            print(f"  ||h_bulk|| = {d['h_bulk_norm']:.6e}")
            print(f"  ||h_bd||   = {d['h_bd_norm']:.6e}")
            print(f"  ||h||      = {d['h_full_norm']:.6e}")
            print(f"  cancel_ratio_h = ||h||/(||h_bulk||+||h_bd||) = {d['cancel_ratio_h']:.6e}")
            print(f"  ||k_E|| = {d['k_E_norm']:.6f}, ||k_B|| = {d['k_B_norm']:.6f}")
            print()

    # Summary
    print("=" * 90)
    print("ZUSAMMENFASSUNG: CANCELLATION IDENTITY")
    print("=" * 90)
    print(f"{'lam':>5} {'N':>4} {'S_bulk':>14} {'S_bd':>14} {'rho_L':>14} {'S_b/S_d':>8} {'cancel':>10}")
    for d in results:
        print(f"{d['lam']:5.1f} {d['N']:4d} {d['S_bulk']:+14.6e} {d['S_bd']:+14.6e} "
              f"{d['rho_L']:+14.6e} {d['sign_ratio']:+8.4f} {d['cancel_ratio_R0f']:10.2e}")

    print()
    print(f"{'lam':>5} {'N':>4} {'mu_E/a':>14} {'mu_B/a':>14} {'mu/a(dir)':>14} {'sign':>8} {'cancel':>10}")
    for d in results:
        print(f"{d['lam']:5.1f} {d['N']:4d} {d['mu_E_over_alpha']:+14.6e} {d['mu_B_over_alpha']:+14.6e} "
              f"{d['mu_alpha']:+14.6e} {d['sign_ratio_mu']:+8.4f} {d['cancel_ratio_mu']:10.2e}")

    print()
    all_cancel = all(d['cancel_ratio_R0f'] < 0.01 for d in results)
    if all_cancel:
        print("ERGEBNIS: STARKE CANCELLATION in allen Faellen!")
        print("  S_bulk ≈ -S_bd mit Residual < 1% der Einzelterme.")
        print("  => Die projected frontier cancellation identity ist numerisch bestaetigt.")
    else:
        strong = sum(1 for d in results if d['cancel_ratio_R0f'] < 0.01)
        print(f"ERGEBNIS: {strong}/{len(results)} Faelle mit starker Cancellation (<1%).")
        for d in results:
            if d['cancel_ratio_R0f'] >= 0.01:
                print(f"  Schwach: lam={d['lam']}, N={d['N']}, ratio={d['cancel_ratio_R0f']:.4e}")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
