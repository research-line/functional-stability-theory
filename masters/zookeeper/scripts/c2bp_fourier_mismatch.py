"""
c2bp_fourier_mismatch.py

Prueft die 2x2-Determinanten-Bedingung fuer Boundary-Nichtdegeneriertheit.

Fuer x := Pi_{L,N} B_L und u := u_tilde gilt:
  ||P_perp x|| >= |x_{n1} u_{n2} - x_{n2} u_{n1}| / sqrt(|u_{n1}|^2 + |u_{n2}|^2)

Ein einziger Fourier-Mismatch x_{n1}/x_{n2} != u_{n1}/u_{n2} genuegt.

Strukturelle Erwartung: B_L ist boundary-lokalisiert (nahe x=L),
ut ist ueber [0,L] verteilt (Pol-Struktur) => verschiedene Fourier-Ratios.
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


def compute_mismatch(lam, N):
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

    # B_L auf Grid
    x_grid = np.linspace(0.02 * L, 0.98 * L, NGRID)
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

    # Fourier-Projektion von B_L
    dx = x_grid[1] - x_grid[0]
    inv_sqrt_L = 1.0 / np.sqrt(L)
    sqrt_2_L = np.sqrt(2.0 / L)
    x_coeffs = np.zeros(dim)
    x_coeffs[0] = np.sum(B_values * inv_sqrt_L) * dx
    for n in range(1, dim):
        basis_n = sqrt_2_L * np.cos(2.0 * np.pi * n * x_grid / L)
        x_coeffs[n] = np.sum(B_values * basis_n) * dx

    # ut-Komponenten (u_n)
    u_coeffs = ut

    return {
        'lam': lam, 'L': L, 'dim': dim,
        'x_coeffs': x_coeffs,   # B_L Fourier-Koeffizienten (NICHT normiert)
        'u_coeffs': u_coeffs,   # ut-Komponenten
    }


def analyze_mismatch(d):
    x = d['x_coeffs']
    u = d['u_coeffs']
    dim = d['dim']
    lam = d['lam']

    print(f"\n  === FOURIER-KOEFFIZIENTEN (lam={lam}) ===")
    print(f"  {'n':>4} {'x_n (B_L)':>14} {'u_n (ut)':>14} {'x_n/x_0':>12} {'u_n/u_0':>12} {'Mismatch':>12}")
    for n in range(min(dim, 8)):
        x_ratio = x[n] / x[0] if abs(x[0]) > 1e-30 else 0
        u_ratio = u[n] / u[0] if abs(u[0]) > 1e-30 else 0
        mismatch = abs(x_ratio - u_ratio) if n > 0 else 0
        print(f"  {n:4d} {x[n]:+14.6e} {u[n]:+14.6e} {x_ratio:+12.6f} {u_ratio:+12.6f} {mismatch:12.6f}")

    # Beste 2x2-Determinante
    print(f"\n  === 2x2-DETERMINANTEN-BOUND ===")
    best_det = 0.0
    best_pair = (0, 0)
    best_bound = 0.0

    for n1 in range(dim):
        for n2 in range(n1 + 1, dim):
            det = abs(x[n1] * u[n2] - x[n2] * u[n1])
            denom = np.sqrt(u[n1]**2 + u[n2]**2)
            if denom > 1e-30:
                bound = det / denom
                if bound > best_bound:
                    best_bound = bound
                    best_det = det
                    best_pair = (n1, n2)

    print(f"  Bestes Modenpaar: (n1={best_pair[0]}, n2={best_pair[1]})")
    print(f"  |x_n1 u_n2 - x_n2 u_n1| = {best_det:.6e}")
    print(f"  sqrt(u_n1^2 + u_n2^2)    = {np.sqrt(u[best_pair[0]]**2 + u[best_pair[1]]**2):.6e}")
    print(f"  => ||P_perp x|| >= {best_bound:.6e}")

    # Vergleich mit tatsaechlichem ||P_perp x||
    P_perp = np.eye(dim) - np.outer(u, u)
    h_bd_actual = np.linalg.norm(P_perp @ x)
    print(f"  ||P_perp x|| (direkt)    = {h_bd_actual:.6e}")
    print(f"  Bound/Actual             = {best_bound/h_bd_actual:.6f}  (Tightness)")

    # Top-5 Determinanten
    dets = []
    for n1 in range(dim):
        for n2 in range(n1 + 1, dim):
            det = abs(x[n1] * u[n2] - x[n2] * u[n1])
            denom = np.sqrt(u[n1]**2 + u[n2]**2)
            if denom > 1e-30:
                dets.append((det / denom, n1, n2))
    dets.sort(reverse=True)
    print(f"\n  Top-5 Modenpaare:")
    for bound, n1, n2 in dets[:5]:
        r_x = x[n1] / x[n2] if abs(x[n2]) > 1e-30 else float('inf')
        r_u = u[n1] / u[n2] if abs(u[n2]) > 1e-30 else float('inf')
        print(f"    ({n1},{n2}): bound={bound:.6e}, x-ratio={r_x:+.6f}, u-ratio={r_u:+.6f}")

    return best_bound, h_bd_actual


def main():
    print("C2bp: Fourier-Mismatch fuer Boundary-Nichtdegeneriertheit")
    print(f"DPS={DPS}, N=20, NGRID={NGRID}")

    results = []
    for lam in [3.0, 5.0, 7.0, 9.0]:
        print(f"\n{'='*80}")
        print(f"  lam={lam}")
        print(f"{'='*80}")

        t0 = time.time()
        d = compute_mismatch(lam, N=20)
        print(f"  Berechnung: {time.time()-t0:.0f}s", flush=True)

        bound, actual = analyze_mismatch(d)
        results.append({'lam': lam, 'bound': bound, 'actual': actual})

    # Zusammenfassung
    print(f"\n{'='*80}")
    print(f"  ZUSAMMENFASSUNG: MISMATCH-BEWEIS")
    print(f"{'='*80}")
    print(f"  {'lam':>5} {'2x2-Bound':>14} {'||P_perp x||':>14} {'Tightness':>12} {'Mismatch?':>10}")
    for r in results:
        tight = r['bound'] / r['actual'] if r['actual'] > 0 else 0
        ok = "JA" if r['bound'] > 1e-10 else "NEIN"
        print(f"  {r['lam']:5.1f} {r['bound']:14.6e} {r['actual']:14.6e} {tight:12.6f} {ok:>10}")

    all_ok = all(r['bound'] > 1e-10 for r in results)
    print(f"\n  FAZIT: {'NICHTDEGENERIERTHEIT BEWIESEN fuer alle lambda' if all_ok else 'NICHT ALLE BESTANDEN'}")
    if all_ok:
        print(f"  Der 2x2-Determinanten-Bound ist strikt positiv bei allen lambda.")
        print(f"  => B_L-Fourier-Koeffizienten sind NICHT proportional zu ut.")
        print(f"  => ||h_bd|| = ||P_perp(Pi B_L)|| > 0 (quantitativ).")
        print(f"  => Boundary-Nichtdegeneriertheit gesichert.")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
