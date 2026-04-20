"""
c2_approximation_test.py — Connes C2/MS2: Approximationsguete k_lambda vs. xi_lambda

MS2 (Missing Step 2): Ist der "Educated Guess" k_lambda eine gute Naeherung
fuer den QW-Grundzustands-Eigenvektor xi_lambda?

CCM 2025 (arXiv:2511.22755), §7:
    h(t) = (sqrt(3)/2^{11/4}) h_4(t) - (3/2^{17/4}) h_0(t)
    h_0(t) = pi^{-1/4} exp(-t^2/2)
    h_4(t) = (384*sqrt(pi))^{-1/2} (16t^4 - 48t^2 + 12) exp(-t^2/2)
    k_lambda(u) = E(h_lambda)(u) = u^{1/2} sum_{n=1}^infty h(n*u/lambda^{-1}... )

    Genauer: h_lambda(t) = lambda^{1/2} h(lambda*t), und
    k_lambda(u) = u^{1/2} sum_{n=1}^infty h_lambda(n*u)
                = u^{1/2} lambda^{1/2} sum_{n=1}^infty h(lambda*n*u)

Fourier-Projektion auf Even-Sektor der V_n-Basis:
    c_m = integral_0^L V_m^*(e^x / lambda) * k_lambda(e^x / lambda) dx
    mit x = log(lambda*u), u = lambda^{-1} e^x, du/u = dx

Messung:
    1. Overlap |<k_lambda, xi_lambda>| / (||k_lambda|| ||xi_lambda||)
    2. Eigengleichungs-Residuum ||QW*k_proj - epsilon*k_proj|| / ||k_proj||
    3. F-Funktion Nullstellen von k_proj vs. xi_lambda

Autor: LG (Opus 4.7 [1M], Session 22)
Datum: 2026-04-19
"""

import os
import sys
import time
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp,
    project_to_parity,
    diagonalize_mp,
    chi_trivial,
    RIEMANN_ZEROS,
)

DPS = int(os.environ.get("DPS", 50))
N = int(os.environ.get("N", 30))

LAMBDAS = [3.0, 5.0, 7.0, 10.0]

# =================================================================
# 1. Hermite-Funktionen h_0, h_4 und der "Educated Guess" h(t)
# =================================================================

def hermite_h0(t):
    """h_0(t) = pi^{-1/4} exp(-t^2/2)"""
    return mp.power(mp.pi, mp.mpf("-0.25")) * mp.exp(-t**2 / 2)

def hermite_h4(t):
    """h_4(t) = (384*sqrt(pi))^{-1/2} (16t^4 - 48t^2 + 12) exp(-t^2/2)"""
    norm = 1 / mp.sqrt(384 * mp.sqrt(mp.pi))
    poly = 16 * t**4 - 48 * t**2 + 12
    return norm * poly * mp.exp(-t**2 / 2)

def h_educated_guess(t):
    """h(t) = (sqrt(3)/2^{11/4}) h_4(t) - (3/2^{17/4}) h_0(t)
    CCM §7: spezifische Kombination mit verschwindendem Integral."""
    c4 = mp.sqrt(3) / mp.power(2, mp.mpf("11") / 4)
    c0 = mp.mpf(3) / mp.power(2, mp.mpf("17") / 4)
    return c4 * hermite_h4(t) - c0 * hermite_h0(t)


# =================================================================
# 2. Epstein-Summe: k_lambda(u) = u^{1/2} sum_{n>=1} h_lambda(n*u)
# =================================================================

def k_lambda_value(u, lam, n_terms=500):
    """k_lambda(u) = u^{1/2} * lambda^{1/2} * sum_{n=1}^{n_terms} h(lambda * n * u).
    h_lambda(t) = lambda^{1/2} h(lambda*t), also
    k_lambda(u) = u^{1/2} sum_n h_lambda(n*u) = u^{1/2} lam^{1/2} sum_n h(lam*n*u).
    """
    lam_mp = mp.mpf(lam)
    u_mp = mp.mpf(u)
    s = mp.mpf(0)
    for n in range(1, n_terms + 1):
        arg = lam_mp * n * u_mp
        val = h_educated_guess(arg)
        s += val
        if n > 5 and abs(val) < mp.mpf("1e-40"):
            break
    return mp.sqrt(u_mp) * mp.sqrt(lam_mp) * s


# =================================================================
# 3. Fourier-Projektion auf Even-Sektor
# =================================================================

def project_k_onto_even_basis(lam, N_val, dps_val):
    """Projiziere k_lambda auf die Even-Fourier-Basis {e_0, e_1, ..., e_N}.

    Even-Basis in u-Raum (L^2([lam^{-1}, lam], du/u)):
      V_0(u) = L^{-1/2}
      V_n(u) = L^{-1/2} exp(2*pi*i*n*log(lam*u)/L)  -->  cos-Kombination fuer even

    Even-Koeffizienten:
      c_0 = int_{lam^{-1}}^{lam} k_lambda(u) * L^{-1/2} du/u
      c_m = sqrt(2) * int_{lam^{-1}}^{lam} k_lambda(u) * L^{-1/2} cos(2*pi*m*log(lam*u)/L) du/u

    Substitution x = log(lam*u), x in [0, L], du/u = dx:
      c_0 = L^{-1/2} int_0^L k_lambda(lam^{-1} e^x) dx
      c_m = sqrt(2) * L^{-1/2} int_0^L k_lambda(lam^{-1} e^x) cos(2*pi*m*x/L) dx
    """
    mp.mp.dps = dps_val
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    inv_sqrt_L = 1 / mp.sqrt(L_mp)
    inv_sqrt2 = 1 / mp.sqrt(mp.mpf(2))
    dim = N_val + 1

    print(f"  Fourier-Projektion: lam={lam}, N={N_val}, L={float(L_mp):.4f}")

    def k_of_x(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(u, lam)

    coeffs = mp.matrix(dim, 1)

    # c_0
    def integrand_0(x):
        return k_of_x(x) * inv_sqrt_L

    coeffs[0, 0] = mp.quad(integrand_0, [0, L_mp], maxdegree=8)
    print(f"    c_0 = {mp.nstr(coeffs[0, 0], 8)}")

    # c_m fuer m = 1..N
    for m in range(1, dim):
        def integrand_m(x, m_val=m):
            return k_of_x(x) * mp.sqrt(mp.mpf(2)) * inv_sqrt_L * mp.cos(2 * mp.pi * m_val * x / L_mp)
        coeffs[m, 0] = mp.quad(integrand_m, [0, L_mp], maxdegree=8)
        if m <= 5 or m == dim - 1:
            print(f"    c_{m} = {mp.nstr(coeffs[m, 0], 8)}")

    return coeffs, L_mp


# =================================================================
# 4. Overlap und Residuum
# =================================================================

def inner_product(a, b, dim):
    """<a, b> im Even-Sektor (reell)."""
    return sum(a[i, 0] * b[i, 0] for i in range(dim))

def norm(a, dim):
    return mp.sqrt(inner_product(a, a, dim))

def evaluate_F_even(coeffs, N_val, L_mp, gamma):
    """F(z) = sum_j xi_j / (z - 2*pi*j/L) in der even-Basis.
    Even-Basis: e_0 -> V_0, e_k -> (V_k + V_{-k})/sqrt(2).
    F(z) = c_0/(z) + sum_{k=1}^N c_k/sqrt(2) * [1/(z-2*pi*k/L) + 1/(z+2*pi*k/L)]
         = c_0/z + sum_{k=1}^N c_k * sqrt(2) * z / (z^2 - (2*pi*k/L)^2)
    """
    g = mp.mpf(gamma)
    F = coeffs[0, 0] / g
    for k in range(1, N_val + 1):
        pole = 2 * mp.pi * k / L_mp
        F += coeffs[k, 0] * mp.sqrt(mp.mpf(2)) * g / (g**2 - pole**2)
    return F


# =================================================================
# 5. Hauptroutine
# =================================================================

def run_for_lambda(lam, N_val, dps_val):
    mp.mp.dps = dps_val
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    dim = N_val + 1

    print(f"\n{'='*80}")
    print(f"  LAMBDA = {lam}, N = {N_val}, L = {float(L_mp):.4f}, DPS = {dps_val}")
    print(f"{'='*80}")

    # --- Phase A: QW Grundzustand xi_lambda ---
    print("\n--- Phase A: QW-Matrix + Grundzustand xi_lambda ---")
    t0 = time.time()
    M_full, _ = build_QW_mp(N_val, lam, chi_func=chi_trivial, q_mod=1,
                            include_W02=True, sign_WR=-1, verbose=False)
    M_even, U_even = project_to_parity(M_full, N_val, parity="even")
    w_sorted, Q_sorted = diagonalize_mp(M_even, verbose=False)
    dt_qw = time.time() - t0

    w_min = w_sorted[0]
    xi = mp.matrix(dim, 1)
    for i in range(dim):
        xi[i, 0] = Q_sorted[i, 0]

    # Cluster identifizieren
    cluster_tol = mp.mpf("1e-10")
    cluster_idx = [k for k in range(dim) if abs(w_sorted[k] - w_min) < cluster_tol]
    cluster_size = len(cluster_idx)
    cluster_edge = cluster_idx[-1] if cluster_size > 1 else 0
    xi_edge = mp.matrix(dim, 1)
    for i in range(dim):
        xi_edge[i, 0] = Q_sorted[i, cluster_edge]

    print(f"  QW aufgebaut + diagonalisiert in {dt_qw:.1f}s")
    print(f"  w_min = {mp.nstr(w_min, 10)}")
    print(f"  Cluster-Groesse: {cluster_size} (Edge-Index k={cluster_edge})")
    print(f"  Niedrigste 6 EWs: {[mp.nstr(w_sorted[k], 6) for k in range(min(6, dim))]}")

    # --- Phase B: Educated Guess k_lambda ---
    print("\n--- Phase B: Educated Guess k_lambda (Hermite + Epstein) ---")
    t0 = time.time()
    k_coeffs, _ = project_k_onto_even_basis(lam, N_val, dps_val)
    dt_k = time.time() - t0
    print(f"  k_lambda projiziert in {dt_k:.1f}s")

    k_norm = norm(k_coeffs, dim)
    print(f"  ||k_lambda|| = {mp.nstr(k_norm, 8)}")

    if k_norm < mp.mpf("1e-30"):
        print("  WARNUNG: k_lambda ist nahezu Null — Epstein-Summe konvergiert nicht bei diesem lambda")
        return None

    # Normiere k_coeffs
    k_normed = mp.matrix(dim, 1)
    for i in range(dim):
        k_normed[i, 0] = k_coeffs[i, 0] / k_norm

    # --- Phase C: Overlap-Analyse ---
    print("\n--- Phase C: Overlap <k_lambda, xi_lambda> ---")

    # Overlap mit Grundzustand (k=0)
    ov_ground = abs(inner_product(k_normed, xi, dim))
    print(f"  |<k_normed, xi_ground>| = {mp.nstr(ov_ground, 8)}")

    # Overlap mit Cluster-Edge
    if cluster_size > 1:
        ov_edge = abs(inner_product(k_normed, xi_edge, dim))
        print(f"  |<k_normed, xi_edge>|   = {mp.nstr(ov_edge, 8)}")
    else:
        ov_edge = ov_ground

    # Overlap mit jedem Cluster-Eigenvektor
    print(f"\n  Overlap mit allen Cluster-EVs:")
    best_k = 0
    best_ov = mp.mpf(0)
    for k in cluster_idx:
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]
        ov = abs(inner_product(k_normed, qk, dim))
        marker = ""
        if ov > best_ov:
            best_ov = ov
            best_k = k
            marker = " <-- BEST"
        print(f"    k={k}: |overlap| = {mp.nstr(ov, 8)}{marker}")

    # Auch einige Nicht-Cluster-EVs testen
    print(f"\n  Overlap mit Nicht-Cluster-EVs:")
    for k in range(cluster_size, min(cluster_size + 5, dim)):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]
        ov = abs(inner_product(k_normed, qk, dim))
        print(f"    k={k}: |overlap| = {mp.nstr(ov, 8)}, w={mp.nstr(w_sorted[k], 6)}")

    # --- Phase D: Eigengleichungs-Residuum ---
    print("\n--- Phase D: Eigengleichungs-Residuum ---")

    # QW * k_normed
    QWk = M_even * k_normed

    # Rayleigh-Quotient: epsilon = <k, QW*k> / <k, k> = <k_normed, QW*k_normed>
    eps_rayleigh = inner_product(k_normed, QWk, dim)
    print(f"  Rayleigh-Quotient eps = <k, QW*k> = {mp.nstr(eps_rayleigh, 10)}")
    print(f"  w_min (Grundzustand)              = {mp.nstr(w_min, 10)}")
    print(f"  Differenz eps - w_min             = {mp.nstr(eps_rayleigh - w_min, 6)}")

    # Residuum ||QW*k - eps*k|| / ||k||
    res = mp.matrix(dim, 1)
    for i in range(dim):
        res[i, 0] = QWk[i, 0] - eps_rayleigh * k_normed[i, 0]
    res_norm = norm(res, dim)
    print(f"  ||QW*k - eps*k|| / ||k|| = {mp.nstr(res_norm, 6)}")

    # --- Phase E: F-Funktion Nullstellen ---
    print("\n--- Phase E: F-Funktion Nullstellen ---")

    print(f"\n  F(gamma) an Riemann-Nullstellen (k_lambda vs xi_edge):")
    print(f"  {'gamma':>12s}  {'F_k(gamma)':>15s}  {'F_xi(gamma)':>15s}  {'Ratio':>12s}")
    for gamma in RIEMANN_ZEROS:
        F_k = evaluate_F_even(k_coeffs, N_val, L_mp, gamma)
        F_xi = evaluate_F_even(xi_edge, N_val, L_mp, gamma)
        ratio = abs(F_k / F_xi) if abs(F_xi) > mp.mpf("1e-30") else mp.inf
        print(f"  {gamma:12.6f}  {float(F_k):15.6e}  {float(F_xi):15.6e}  {float(ratio):12.4e}")

    # F an Nicht-Nullstellen zum Vergleich
    non_zeros = [10.0, 17.5, 23.0, 28.0]
    print(f"\n  F(z) an Nicht-Nullstellen:")
    print(f"  {'z':>12s}  {'F_k(z)':>15s}  {'F_xi(z)':>15s}")
    for z in non_zeros:
        F_k = evaluate_F_even(k_coeffs, N_val, L_mp, z)
        F_xi = evaluate_F_even(xi_edge, N_val, L_mp, z)
        print(f"  {z:12.6f}  {float(F_k):15.6e}  {float(F_xi):15.6e}")

    # --- Phase F: Koeffizienten-Vergleich ---
    print("\n--- Phase F: Koeffizienten-Vergleich (erste 10) ---")
    xi_use = xi_edge
    xi_norm_val = norm(xi_use, dim)
    print(f"  {'m':>4s}  {'k_normed[m]':>15s}  {'xi_normed[m]':>15s}  {'Diff':>12s}")
    for m in range(min(10, dim)):
        km = k_normed[m, 0]
        xm = xi_use[m, 0] / xi_norm_val
        print(f"  {m:4d}  {float(km):15.8e}  {float(xm):15.8e}  {float(abs(km - xm)):12.4e}")

    return {
        "lambda": lam,
        "w_min": float(w_min),
        "cluster_size": cluster_size,
        "best_overlap_k": best_k,
        "best_overlap": float(best_ov),
        "rayleigh": float(eps_rayleigh),
        "residual": float(res_norm),
        "k_norm": float(k_norm),
    }


def main():
    mp.mp.dps = DPS
    print("=" * 80)
    print("C2/MS2: APPROXIMATIONSGUETE k_lambda vs. xi_lambda")
    print("=" * 80)
    print(f"  Hermite-Kombination: h = (sqrt(3)/2^(11/4)) h_4 - (3/2^(17/4)) h_0")
    print(f"  Epstein-Summe: k_lambda(u) = u^(1/2) lam^(1/2) sum_n h(lam*n*u)")
    print(f"  N = {N}, DPS = {DPS}")
    print(f"  Lambdas: {LAMBDAS}")

    # Vorab-Check: h(t) hat verschwindendes Integral
    int_h = mp.quad(h_educated_guess, [-20, 20])
    print(f"\n  Integral-Check: int h(t) dt = {mp.nstr(int_h, 6)} (soll ~0)")

    results = []
    for lam in LAMBDAS:
        r = run_for_lambda(lam, N, DPS)
        if r:
            results.append(r)

    # Zusammenfassung
    print(f"\n{'='*80}")
    print(f"  ZUSAMMENFASSUNG")
    print(f"{'='*80}")
    print(f"  {'lambda':>8s}  {'w_min':>12s}  {'Cluster':>8s}  {'best_k':>7s}  {'Overlap':>10s}  {'Rayleigh':>12s}  {'Residual':>12s}")
    for r in results:
        print(f"  {r['lambda']:8.1f}  {r['w_min']:12.6f}  {r['cluster_size']:8d}  {r['best_overlap_k']:7d}  {r['best_overlap']:10.6f}  {r['rayleigh']:12.6f}  {r['residual']:12.4e}")


if __name__ == "__main__":
    main()
