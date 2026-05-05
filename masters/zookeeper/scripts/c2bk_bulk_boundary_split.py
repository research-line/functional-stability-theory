"""
c2bk_bulk_boundary_split.py

Zerlegt h = P_0 (H - w0) k in Bulk- und Boundary-Anteil:

  h = h_bulk + h_boundary

wobei:
  h_boundary = P_0 * Fourier(B_L) / n_{L,N}
  h_bulk     = h - h_boundary

B_L ist der Boundary-Defekt aus C2M (H_L K_L = (H_inf K)|_[0,L] + B_L).
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
NGRID = 30


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


def compute_split(lam, N):
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

    ws_A, _ = np.linalg.eigh(A_np)
    w0 = ws_A[0]

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
    alpha = float(ut @ kn)

    P_perp = np.eye(dim) - np.outer(ut, ut)

    # Full residual: res = (H - w0) k
    res = H_np @ kn - w0 * kn

    # h = P_perp * res
    h_full = P_perp @ res

    # B_L computation on grid (from c2bd_derivative_matching)
    x_grid = np.linspace(0.02 * L, 0.98 * L, NGRID)
    B_values = np.zeros(NGRID)

    for i, x in enumerate(x_grid):
        x_mp = mp.mpf(x)
        ex = mp.exp(x_mp)
        emx = mp.exp(-x_mp)

        # Prime boundary
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

        # Archimedes boundary
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
    b_norm = b_raw / nf_float  # boundary vector in k-normalized Fourier space

    # Split
    h_boundary = P_perp @ b_norm
    h_bulk = h_full - h_boundary

    return {
        'lam': lam, 'L': L, 'w0': w0, 'alpha': alpha,
        'h_full': h_full, 'h_bulk': h_bulk, 'h_boundary': h_boundary,
        'h_norm': np.linalg.norm(h_full),
        'h_bulk_norm': np.linalg.norm(h_bulk),
        'h_boundary_norm': np.linalg.norm(h_boundary),
        'res': res, 'b_norm': b_norm, 'ut': ut, 'kn': kn,
        'mu_alpha': float(ut @ (A_np - w0 * np.eye(dim)) @ kn) / alpha,
    }


def main():
    print(f"C2bk: Bulk/Boundary-Split von h")
    print(f"DPS={DPS}, NGRID={NGRID}")

    results = []
    for lam in [3.0, 5.0, 7.0]:
        print(f"\n{'='*80}")
        print(f"  lam={lam}")
        print(f"{'='*80}")

        t0 = time.time()
        d = compute_split(lam, N=20)
        print(f"  Berechnung: {time.time()-t0:.0f}s", flush=True)
        results.append(d)

        print(f"\n  --- NORMEN ---")
        print(f"  ||h||          = {d['h_norm']:.6e}")
        print(f"  ||h_bulk||     = {d['h_bulk_norm']:.6e}")
        print(f"  ||h_boundary|| = {d['h_boundary_norm']:.6e}")
        print(f"  ||h_bulk|| + ||h_boundary|| = {d['h_bulk_norm'] + d['h_boundary_norm']:.6e}")
        ratio_b = d['h_boundary_norm'] / d['h_norm'] if d['h_norm'] > 0 else 0
        ratio_k = d['h_bulk_norm'] / d['h_norm'] if d['h_norm'] > 0 else 0
        print(f"  Boundary-Anteil: {ratio_b*100:.1f}%")
        print(f"  Bulk-Anteil:     {ratio_k*100:.1f}%")

        # ut-Projektionen (Skalarwerte)
        h_u = float(d['ut'] @ d['h_full'])
        hbulk_u = float(d['ut'] @ d['h_bulk'])
        hbnd_u = float(d['ut'] @ d['h_boundary'])
        print(f"\n  --- ut-PROJEKTIONEN ---")
        print(f"  <ut|h>          = {h_u:+.6e}  (sollte ~0 wegen P_perp)")
        print(f"  <ut|h_bulk>     = {hbulk_u:+.6e}")
        print(f"  <ut|h_boundary> = {hbnd_u:+.6e}")

        # Cosine angle between h_bulk and h_boundary
        if d['h_bulk_norm'] > 1e-40 and d['h_boundary_norm'] > 1e-40:
            cos_angle = float(d['h_bulk'] @ d['h_boundary']) / (d['h_bulk_norm'] * d['h_boundary_norm'])
            print(f"\n  --- GEOMETRIE ---")
            print(f"  cos(h_bulk, h_boundary) = {cos_angle:+.6f}")
            print(f"  h_bulk . h_boundary     = {float(d['h_bulk'] @ d['h_boundary']):+.6e}")

        # Verification
        print(f"\n  --- VERIFIKATION ---")
        print(f"  mu/alpha (direkt) = {d['mu_alpha']:+.6e}")

    # Cross-lambda summary
    print(f"\n{'='*80}")
    print(f"  ZUSAMMENFASSUNG")
    print(f"{'='*80}")
    print(f"  {'lam':>5} {'||h||':>12} {'||h_bulk||':>12} {'||h_bnd||':>12} {'%Bulk':>8} {'%Bnd':>8}")
    for d in results:
        pb = d['h_boundary_norm'] / d['h_norm'] * 100
        pk = d['h_bulk_norm'] / d['h_norm'] * 100
        print(f"  {d['lam']:5.1f} {d['h_norm']:12.4e} {d['h_bulk_norm']:12.4e} "
              f"{d['h_boundary_norm']:12.4e} {pk:8.1f} {pb:8.1f}")

    print(f"\n  Skalierung mit lambda:")
    for i in range(len(results) - 1):
        r1, r2 = results[i], results[i+1]
        if r1['h_bulk_norm'] > 1e-40 and r2['h_bulk_norm'] > 1e-40:
            p_bulk = np.log(r2['h_bulk_norm'] / r1['h_bulk_norm']) / np.log(r2['lam'] / r1['lam'])
            p_bnd = np.log(r2['h_boundary_norm'] / r1['h_boundary_norm']) / np.log(r2['lam'] / r1['lam'])
            print(f"    lam={r1['lam']:.0f}->{r2['lam']:.0f}: "
                  f"||h_bulk|| ~ lam^{p_bulk:.2f}, ||h_bnd|| ~ lam^{p_bnd:.2f}")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
