"""
c2bc_bulk_boundary_decomposition.py

Bulk/Boundary-Zerlegung von T2 = <f_lam, k_perp>/alpha.

Aus C2bb: mu/alpha = T1 + T2, T1 = w_bar - w_min, T2 = <f_lam, k_perp>/alpha.
Aus C2m.3: (H_L - w*) K_L = (H_inf - w*) K|[0,L] + B_L mit B_L >= 0.

Zerlegung des Residuals res = (H - w_min)k in Fourier-Raum:
  res = E_bulk + boundary
  <u_tilde, res> = <u_tilde, E_bulk> + <u_tilde, boundary>
  mu/alpha = C_tilde + <u_tilde, res>/alpha

Conjecture (C2bb.3):
  Der Bulk-Term E_bulk traegt -(C_tilde + epsilon), so dass
  mu/alpha = epsilon_bulk + <u_tilde, boundary>/alpha.
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import mpmath as mp
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, diagonalize_mp, chi_trivial,
    rho_even_mp, mangoldt_mp,
)
from c2_approximation_test import (
    h_educated_guess, k_lambda_value, inner_product, norm,
)
from c2_poisson_decomposition import (
    h_hat_analytical, project_to_fourier,
)

DPS = int(os.environ.get("DPS", 50))
NGRID = 100


def K_direct(z_mp):
    ez = mp.exp(z_mp)
    s = mp.mpf(0)
    for n in range(1, 500):
        val = h_educated_guess(mp.mpf(n) * ez)
        s += val
        if n > 3 and abs(val) < mp.mpf("1e-45"):
            break
    return mp.exp(z_mp / 2) * s


def K_poisson(z_mp):
    emz = mp.exp(-z_mp)
    s = mp.mpf(0)
    for m in range(1, 500):
        val = h_hat_analytical(mp.mpf(m) * emz)
        s += val
        if m > 3 and abs(val) < mp.mpf("1e-45"):
            break
    return mp.exp(-z_mp / 2) * s


def K_global(z):
    z_mp = mp.mpf(z)
    return K_direct(z_mp) if z_mp >= 0 else K_poisson(z_mp)


def E_L_plus(m, x_mp, L_mp):
    threshold = float(mp.exp(L_mp - x_mp))
    s = mp.mpf(0)
    for k in range(2, m + 1):
        if m % k == 0:
            lk = mangoldt_mp(k)
            if lk > 0 and k > threshold:
                s += lk
    return s


def E_L_minus(r, x_mp):
    threshold = float(mp.exp(x_mp))
    s = mp.mpf(0)
    for k in range(2, r + 1):
        if r % k == 0:
            lk = mangoldt_mp(k)
            if lk > 0 and k > threshold:
                s += lk
    return s


def B_prime_L(x_mp, L_mp):
    ex = mp.exp(x_mp)
    emx = mp.exp(-x_mp)
    s1 = mp.mpf(0)
    for m in range(1, 5000):
        arg = mp.mpf(m) * ex
        hval = h_educated_guess(arg)
        if m > 3 and abs(hval) < mp.mpf("1e-45"):
            break
        ep = E_L_plus(m, x_mp, L_mp)
        if ep > 0:
            s1 += ep * hval
    s1 *= mp.exp(x_mp / 2)
    s2 = mp.mpf(0)
    for r in range(1, 5000):
        arg = mp.mpf(r) * emx
        hhval = h_hat_analytical(arg)
        if r > 3 and abs(hhval) < mp.mpf("1e-45"):
            break
        em = E_L_minus(r, x_mp)
        if em > 0:
            s2 += em * hhval
    s2 *= mp.exp(-x_mp / 2)
    return s1 + s2


Y_MAX = 12.0

def B_R_L(x_mp, L_mp):
    y_max_mp = mp.mpf(Y_MAX)
    a1 = L_mp - x_mp
    if a1 < mp.mpf("0.01"):
        a1 = mp.mpf("0.01")
    def int1(y):
        return rho_even_mp(y) * K_global(float(x_mp + y))
    I1 = mp.mpf(0)
    if a1 < y_max_mp:
        try:
            I1 = mp.quad(int1, [a1, y_max_mp], maxdegree=6)
        except Exception:
            I1 = mp.mpf(0)
    def int2(y):
        return rho_even_mp(y) * K_global(float(x_mp - y))
    try:
        I2 = mp.quad(int2, [x_mp, y_max_mp], maxdegree=6)
    except Exception:
        I2 = mp.mpf(0)
    return I1, I2


def run(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

    print(f"\n{'='*90}")
    print(f"  C2bc BULK/BOUNDARY DECOMPOSITION OF T2: lam={lam}, N={N}")
    print(f"{'='*90}")

    t0 = time.time()

    # === Phase 1: Operatoren, C2bb-Groessen ===
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

    Ctilde = float(sum(ut[i, 0] * sum(W02[i, j] * ut[j, 0] for j in range(dim))
                       for i in range(dim)))

    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))
    w_bar = Ctilde + uHu

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
    alpha = float(inner_product(ut, kn, dim))

    k_perp = mp.matrix(dim, 1)
    for i in range(dim):
        k_perp[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]

    f_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]

    ws, _ = diagonalize_mp(Aq, verbose=False)
    w_min = float(ws[0])

    # Residual res = (H - w_min)k
    Hk = mp.matrix(dim, 1)
    for i in range(dim):
        Hk[i, 0] = sum(H[i, j] * kn[j, 0] for j in range(dim))

    res = mp.matrix(dim, 1)
    for i in range(dim):
        res[i, 0] = Hk[i, 0] - mp.mpf(w_min) * kn[i, 0]

    res_u = float(inner_product(ut, res, dim))

    # C2bb identity
    T1 = w_bar - w_min
    T2 = float(inner_product(f_lam, k_perp, dim)) / alpha
    mu_alpha = T1 + T2

    print(f"  Phase 1 in {time.time()-t0:.0f}s", flush=True)
    print(f"  C_tilde      = {Ctilde:.6f}", flush=True)
    print(f"  uHu          = {uHu:.6f}", flush=True)
    print(f"  w_bar        = {w_bar:.6f}", flush=True)
    print(f"  w_min        = {w_min:.6f}", flush=True)
    print(f"  alpha        = {alpha:.6f}", flush=True)
    print(f"  T1           = {T1:+.10e}", flush=True)
    print(f"  T2           = {T2:+.10e}", flush=True)
    print(f"  mu/alpha     = {mu_alpha:+.10e}", flush=True)
    print(f"  <ut, res>    = {res_u:+.10e}", flush=True)
    print(f"  C+<ut,res>/a = {Ctilde + res_u/alpha:+.10e}  (soll = mu/alpha)", flush=True)

    # === Phase 2: B_L auf Grid berechnen ===
    print(f"\n  --- B_L Grid ({NGRID} Punkte) ---", flush=True)
    t1 = time.time()
    x_grid = np.linspace(0.02 * L, 0.98 * L, NGRID)
    B_values = np.zeros(NGRID)

    for i, x in enumerate(x_grid):
        x_mp = mp.mpf(x)
        bp = B_prime_L(x_mp, L_mp)
        br1, br2 = B_R_L(x_mp, L_mp)
        B_values[i] = float(bp + br1 + br2)
        if (i + 1) % 25 == 0:
            print(f"    {i+1}/{NGRID} ({time.time()-t1:.0f}s)", flush=True)

    print(f"  Grid in {time.time()-t1:.0f}s", flush=True)
    print(f"  max|B_L|     = {np.max(np.abs(B_values)):.6e}", flush=True)
    print(f"  B_L > 0      = {np.all(B_values >= -1e-20)}", flush=True)

    # === Phase 3: Fourier-Projektion von B_L ===
    dx = x_grid[1] - x_grid[0]
    inv_sqrt_L = 1.0 / np.sqrt(L)
    sqrt_2_over_L = np.sqrt(2.0 / L)

    b_raw = np.zeros(dim)
    b_raw[0] = np.sum(B_values * inv_sqrt_L) * dx
    for n in range(1, dim):
        basis_n = sqrt_2_over_L * np.cos(2.0 * np.pi * n * x_grid / L)
        b_raw[n] = np.sum(B_values * basis_n) * dx

    nf_float = float(nf)
    b_norm = b_raw / nf_float

    b_mp = mp.matrix(dim, 1)
    for i in range(dim):
        b_mp[i, 0] = mp.mpf(b_norm[i])

    boundary_u = float(inner_product(ut, b_mp, dim))
    ebulk_u = res_u - boundary_u

    # === Phase 4: Zerlegung ===
    print(f"\n  --- u_tilde-Zerlegung ---", flush=True)
    print(f"  <ut, res>        = {res_u:+.10e}  (total)", flush=True)
    print(f"  <ut, E_bulk>     = {ebulk_u:+.10e}  (Bulk)", flush=True)
    print(f"  <ut, boundary>   = {boundary_u:+.10e}  (Boundary B_L)", flush=True)
    print(f"  Summe            = {ebulk_u + boundary_u:+.10e}  (soll = <ut,res>)", flush=True)

    print(f"\n  --- mu/alpha Zerlegung ---", flush=True)
    term_C = Ctilde
    term_bulk = ebulk_u / alpha
    term_boundary = boundary_u / alpha

    print(f"  C_tilde          = {term_C:+.10e}", flush=True)
    print(f"  <ut,E_bulk>/a    = {term_bulk:+.10e}", flush=True)
    print(f"  <ut,boundary>/a  = {term_boundary:+.10e}", flush=True)
    print(f"  Summe            = {term_C + term_bulk + term_boundary:+.10e}", flush=True)
    print(f"  mu/alpha (C2bb)  = {mu_alpha:+.10e}", flush=True)
    print(f"  Konsistenz       = {abs(term_C + term_bulk + term_boundary - mu_alpha):.2e}", flush=True)

    # T2-Zerlegung
    T2_bulk = term_bulk + Ctilde - (Ctilde + uHu - w_min)
    T2_boundary = term_boundary

    # Alternative: T2 = <ut, res>/alpha + C - T1 = (ebulk_u + boundary_u)/alpha + C - T1
    # T2_via_bulk = ebulk_u/alpha + C - T1 = ebulk_u/alpha + C - (C + uHu - w_min)
    #             = ebulk_u/alpha - uHu + w_min
    # T2_via_boundary = boundary_u/alpha

    T2_from_bulk = ebulk_u / alpha - uHu + w_min
    T2_from_boundary = boundary_u / alpha

    print(f"\n  --- T2 Zerlegung (Vermutung C2bb.3) ---", flush=True)
    print(f"  T2 (total)       = {T2:+.10e}", flush=True)
    print(f"  T2_bulk          = {T2_from_bulk:+.10e}", flush=True)
    print(f"  T2_boundary      = {T2_from_boundary:+.10e}", flush=True)
    print(f"  T2_bulk + T2_bd  = {T2_from_bulk + T2_from_boundary:+.10e}", flush=True)

    print(f"\n  --- Vermutung C2bb.3: T2_bulk = -T1? ---", flush=True)
    print(f"  T1               = {T1:+.10e}", flush=True)
    print(f"  T2_bulk          = {T2_from_bulk:+.10e}", flush=True)
    print(f"  T2_bulk + T1     = {T2_from_bulk + T1:+.10e}  (soll klein)", flush=True)
    print(f"  Ratio            = {abs(T2_from_bulk + T1) / abs(T1):.6e}", flush=True)

    print(f"\n  --- Vermutung: mu/alpha = T2_boundary? ---", flush=True)
    print(f"  mu/alpha         = {mu_alpha:+.10e}", flush=True)
    print(f"  T2_boundary      = {T2_from_boundary:+.10e}", flush=True)
    print(f"  Diff             = {abs(mu_alpha - T2_from_boundary):.6e}", flush=True)

    # Anteil
    total_abs = abs(T2_from_bulk) + abs(T2_from_boundary)
    print(f"\n  --- Anteile ---", flush=True)
    print(f"  |T2_bulk|/total  = {abs(T2_from_bulk)/total_abs:.6f}  ({abs(T2_from_bulk)/total_abs*100:.1f}%)", flush=True)
    print(f"  |T2_bd|/total    = {abs(T2_from_boundary)/total_abs:.6f}  ({abs(T2_from_boundary)/total_abs*100:.1f}%)", flush=True)

    return {
        'lam': lam, 'T1': T1, 'T2': T2, 'mu_alpha': mu_alpha,
        'Ctilde': Ctilde, 'uHu': uHu, 'w_bar': w_bar, 'w_min': w_min,
        'alpha': alpha, 'res_u': res_u,
        'ebulk_u': ebulk_u, 'boundary_u': boundary_u,
        'T2_bulk': T2_from_bulk, 'T2_boundary': T2_from_boundary,
    }


def main():
    print(f"C2bc: Bulk/Boundary Decomposition of T2", flush=True)
    print(f"DPS = {DPS}, NGRID = {NGRID}", flush=True)

    results = []
    for cfg in [{"lam": 3.0, "N": 30}, {"lam": 5.0, "N": 55}]:
        r = run(cfg['lam'], cfg['N'])
        results.append(r)

    print(f"\n{'='*90}")
    print(f"ZUSAMMENFASSUNG", flush=True)
    print(f"{'='*90}")
    for r in results:
        print(f"\n  lam = {r['lam']}", flush=True)
        print(f"  T1 = {r['T1']:+.6e}, T2 = {r['T2']:+.6e}, mu/a = {r['mu_alpha']:+.6e}", flush=True)
        print(f"  T2_bulk = {r['T2_bulk']:+.6e} ({abs(r['T2_bulk'])/abs(r['T2'])*100:.1f}% von T2)", flush=True)
        print(f"  T2_bd   = {r['T2_boundary']:+.6e} ({abs(r['T2_boundary'])/abs(r['T2'])*100:.1f}% von T2)", flush=True)
        print(f"  T2_bulk + T1 = {r['T2_bulk']+r['T1']:+.6e} (Cancellation-Rest)", flush=True)
        print(f"  mu/a - T2_bd = {r['mu_alpha']-r['T2_boundary']:+.6e}", flush=True)


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
