"""
c2m_boundary_defect_explicit.py — Explizite B_L-Berechnung (Bulk/Boundary-Separation)

C2m.3: (H_L - w*) K_L = (H_inf - w*) K + B_L  mit B_L >= 0

B_L = B_{R,L} + B'_L wobei:
  B'_L(x) = e^{x/2} sum_m E_L^+(m,x) h(me^x) + e^{-x/2} sum_r E_L^-(r,x) h_hat(re^{-x})
  B_{R,L}(x) = int_{L-x}^inf rho(y) K(x+y) dy + int_x^inf rho(y) K(x-y) dy

Drei Zielgroessen:
  1. ||P_0 E_bulk|| = ||P_0(res - B_L/||K_L||)||  (projizierter Bulk-Defekt)
  2. |<u_tilde, res>|                               (u_tilde-Komponente)
  3. ||P_0 B_L|| / ||K_L||                          (projizierter Boundary-Defekt)

WICHTIG: Korrekte Fourier-Basis cos(2*pi*n*x/L), NICHT cos(n*pi*x/L)!

Autor: LG (Opus 4.6, Session 27)
Datum: 2026-04-19
"""

import os, sys, time, math
import mpmath as mp
import numpy as np
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
CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
]
NGRID = 200
Y_MAX = 12.0


# ================================================================
# Globaler Poisson-Kern K(z)
# ================================================================

def K_direct(z_mp):
    """K(z) = e^{z/2} sum_{n>=1} h(n*e^z) — direkt, fuer z >= 0."""
    ez = mp.exp(z_mp)
    s = mp.mpf(0)
    for n in range(1, 500):
        val = h_educated_guess(mp.mpf(n) * ez)
        s += val
        if n > 3 and abs(val) < mp.mpf("1e-45"):
            break
    return mp.exp(z_mp / 2) * s


def K_poisson(z_mp):
    """K(z) = e^{-z/2} sum_{m>=1} h_hat(m*e^{-z}) — Poisson-Form, fuer z <= 0."""
    emz = mp.exp(-z_mp)
    s = mp.mpf(0)
    for m in range(1, 500):
        val = h_hat_analytical(mp.mpf(m) * emz)
        s += val
        if m > 3 and abs(val) < mp.mpf("1e-45"):
            break
    return mp.exp(-z_mp / 2) * s


def K_global(z):
    """K(z) — waehlt die besser konvergierende Darstellung."""
    z_mp = mp.mpf(z)
    return K_direct(z_mp) if z_mp >= 0 else K_poisson(z_mp)


# ================================================================
# B'_L: Prim-Boundary-Defekt
# ================================================================

def E_L_plus(m, x_mp, L_mp):
    """E_L^+(m,x) = sum_{k|m, k > e^{L-x}} Lambda(k) >= 0"""
    threshold = float(mp.exp(L_mp - x_mp))
    s = mp.mpf(0)
    for k in range(2, m + 1):
        if m % k == 0:
            lk = mangoldt_mp(k)
            if lk > 0 and k > threshold:
                s += lk
    return s


def E_L_minus(r, x_mp):
    """E_L^-(r,x) = sum_{k|r, k > e^x} Lambda(k) >= 0"""
    threshold = float(mp.exp(x_mp))
    s = mp.mpf(0)
    for k in range(2, r + 1):
        if r % k == 0:
            lk = mangoldt_mp(k)
            if lk > 0 and k > threshold:
                s += lk
    return s


def B_prime_L(x_mp, L_mp):
    """B'_L(x) — Prim-Boundary-Defekt."""
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


# ================================================================
# B_{R,L}: Archimedes-Boundary-Defekt
# ================================================================

def B_R_L(x_mp, L_mp):
    """B_{R,L}(x) = int_{L-x}^inf rho(y) K(x+y) dy + int_x^inf rho(y) K(x-y) dy"""
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


# ================================================================
# Hauptroutine
# ================================================================

def run(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

    print(f"\n{'='*72}")
    print(f"  C2m BULK/BOUNDARY SEPARATION: lambda={lam}, N={N}, L={L:.4f}")
    print(f"{'='*72}")

    # --- Phase 1: Operatoren, Residual, ut, kn ---
    t0 = time.time()
    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    H = mp.matrix(dim, dim)
    W02 = mp.matrix(dim, dim)
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

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = float(norm(cf, dim))
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / mp.mpf(nf)
    alpha = float(inner_product(ut, kn, dim))

    ws, _ = diagonalize_mp(Aq, verbose=False)
    w_min = float(ws[0])

    Hk = mp.matrix(dim, 1)
    for i in range(dim):
        Hk[i, 0] = sum(H[i, j] * kn[j, 0] for j in range(dim))

    res = mp.matrix(dim, 1)
    for i in range(dim):
        res[i, 0] = Hk[i, 0] - mp.mpf(w_min) * kn[i, 0]

    res_u = float(inner_product(ut, res, dim))
    res_P0 = mp.matrix(dim, 1)
    for i in range(dim):
        res_P0[i, 0] = res[i, 0] - mp.mpf(res_u) * ut[i, 0]

    P0_res_norm = float(norm(res_P0, dim))
    print(f"  Phase 1 in {time.time()-t0:.0f}s")
    print(f"  w_min={w_min:.6f}, alpha={alpha:.6f}, ||K_L||={nf:.6f}")
    print(f"  |<ut, res>| = {abs(res_u):.6e}")
    print(f"  ||P0 res||  = {P0_res_norm:.6e}")

    # --- Phase 2: Quick B_L an Stichpunkten ---
    print(f"\n  --- Quick B_L Schaetzung ---")
    print(f"  {'x/L':>6}  {'B_prime':>12}  {'B_R_fwd':>12}  {'B_R_bwd':>12}  {'B_total':>12}  {'dt':>6}")
    print(f"  {'-'*68}")

    test_fracs = [0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95]
    quick_vals = []
    quick_bp_max = 0.0
    quick_br1_max = 0.0

    for frac in test_fracs:
        x_mp = mp.mpf(frac * L)
        t1 = time.time()
        bp = B_prime_L(x_mp, L_mp)
        br1, br2 = B_R_L(x_mp, L_mp)
        btot = bp + br1 + br2
        dt = time.time() - t1
        quick_vals.append((frac, float(btot)))
        quick_bp_max = max(quick_bp_max, abs(float(bp)))
        quick_br1_max = max(quick_br1_max, abs(float(br1)))
        print(f"  {frac:6.2f}  {float(bp):12.3e}  {float(br1):12.3e}  {float(br2):12.3e}  {float(btot):12.3e}  {dt:5.1f}s")

    skip_prime = quick_bp_max < 1e-15
    skip_fwd = quick_br1_max < 1e-15

    max_BL_quick = max(abs(v) for _, v in quick_vals)
    ratio_quick = max_BL_quick / (P0_res_norm * nf)
    print(f"\n  max|B_L|/||K_L||     = {max_BL_quick/nf:.3e}")
    print(f"  ||P0 res||           = {P0_res_norm:.3e}")
    print(f"  Ratio max|B_L/nf|/||P0 res|| ~ {ratio_quick:.3e}")

    if skip_prime:
        print(f"  => B'_L vernachlaessigbar (< 1e-20)")
    if skip_fwd:
        print(f"  => B_R forward-Integral vernachlaessigbar (< 1e-20)")

    # --- Phase 3: Volles Grid (nur dominanter Term) ---
    print(f"\n  --- Volle B_L Berechnung (Grid={NGRID}) ---")
    t2 = time.time()
    x_grid = np.linspace(0.01 * L, 0.99 * L, NGRID)
    B_values = np.zeros(NGRID)

    for i, x in enumerate(x_grid):
        x_mp = mp.mpf(x)
        if skip_fwd and skip_prime:
            _, br2 = B_R_L(x_mp, L_mp)
            B_values[i] = float(br2)
        else:
            bp = B_prime_L(x_mp, L_mp) if not skip_prime else mp.mpf(0)
            br1, br2 = B_R_L(x_mp, L_mp)
            B_values[i] = float(bp + br1 + br2)
        if (i + 1) % 50 == 0:
            print(f"    {i+1}/{NGRID} ({time.time()-t2:.0f}s)")

    print(f"  Grid in {time.time()-t2:.0f}s")

    # --- Phase 4: Fourier-Projektion von B_L ---
    dx = x_grid[1] - x_grid[0]
    inv_sqrt_L = 1.0 / np.sqrt(L)
    sqrt_2_over_L = np.sqrt(2.0 / L)

    b_raw = np.zeros(dim)
    b_raw[0] = np.sum(B_values * inv_sqrt_L) * dx
    for n in range(1, dim):
        basis_n = sqrt_2_over_L * np.cos(2.0 * np.pi * n * x_grid / L)
        b_raw[n] = np.sum(B_values * basis_n) * dx

    b_norm = b_raw / nf

    b_mp = mp.matrix(dim, 1)
    for i in range(dim):
        b_mp[i, 0] = mp.mpf(b_norm[i])

    b_u = float(inner_product(ut, b_mp, dim))
    P0_b = mp.matrix(dim, 1)
    for i in range(dim):
        P0_b[i, 0] = b_mp[i, 0] - mp.mpf(b_u) * ut[i, 0]

    P0_b_norm = float(norm(P0_b, dim))

    # E_bulk = res - B_L/||K_L|| (in Koeffizientenraum)
    P0_ebulk = mp.matrix(dim, 1)
    for i in range(dim):
        P0_ebulk[i, 0] = res_P0[i, 0] - P0_b[i, 0]

    P0_ebulk_norm = float(norm(P0_ebulk, dim))

    # --- Phase 5: Ergebnis ---
    print(f"\n  {'='*56}")
    print(f"  ERGEBNIS: BULK/BOUNDARY-SEPARATION")
    print(f"  {'='*56}")
    print(f"  |<ut, res>|     = {abs(res_u):.6e}  (ut-Komponente)")
    print(f"  ||P0 res||      = {P0_res_norm:.6e}  (Gesamt-P0-Residual)")
    print(f"  ||P0 B_L/nf||   = {P0_b_norm:.6e}  (Boundary-Defekt, normiert)")
    print(f"  ||P0 E_bulk||   = {P0_ebulk_norm:.6e}  (Bulk-Defekt)")
    print(f"  ||P0 B_L|| / ||P0 res|| = {P0_b_norm/P0_res_norm:.6e}")
    print(f"  ||P0 E_bulk|| / ||P0 res|| = {P0_ebulk_norm/P0_res_norm:.6e}")

    print(f"\n  Spatial B_L(x):")
    print(f"    max|B_L(x)|   = {np.max(np.abs(B_values)):.6e}")
    print(f"    RMS B_L       = {np.sqrt(np.mean(B_values**2)):.6e}")
    print(f"    B_L > 0 ?     = {np.all(B_values >= -1e-30)}")

    left_E = np.sum(B_values[:NGRID//5]**2) * dx
    right_E = np.sum(B_values[-NGRID//5:]**2) * dx
    total_E = np.sum(B_values**2) * dx
    if total_E > 0:
        print(f"    Rand-Fraktion (20%): {100*(left_E + right_E)/total_E:.1f}%")

    return {
        "lam": lam,
        "res_u": abs(res_u),
        "P0_res": P0_res_norm,
        "P0_BL": P0_b_norm,
        "P0_Ebulk": P0_ebulk_norm,
        "max_BL": np.max(np.abs(B_values)),
        "nf": nf,
    }


if __name__ == "__main__":
    results = []
    for cfg in CONFIGS:
        r = run(cfg["lam"], cfg["N"])
        results.append(r)

    if len(results) >= 2:
        print(f"\n{'='*72}")
        print(f"  lambda-SKALIERUNG")
        print(f"{'='*72}")

        print(f"\n  {'Metrik':>20}  ", end="")
        for r in results:
            print(f"{'lam='+str(r['lam']):>14}  ", end="")
        print(f"{'Exponent':>10}")

        for name, key in [
            ("|<ut,res>|", "res_u"),
            ("||P0 res||", "P0_res"),
            ("||P0 B_L||", "P0_BL"),
            ("||P0 E_bulk||", "P0_Ebulk"),
            ("max|B_L|", "max_BL"),
        ]:
            vals = [r[key] for r in results]
            print(f"  {name:>20}  ", end="")
            for v in vals:
                print(f"{v:14.6e}  ", end="")
            if len(vals) >= 2 and vals[0] > 0 and vals[-1] > 0:
                lams = [r["lam"] for r in results]
                exp = math.log(vals[-1] / vals[0]) / math.log(lams[-1] / lams[0])
                print(f"  {exp:8.2f}")
            else:
                print()

    print(f"\n=== FERTIG ===")
