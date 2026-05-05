#!/usr/bin/env python3
"""
c2_server_leakage.py — Server-Lauf: Lambda-Skalierung der C2/MS2-Leckage mit N proportional lambda^2

Kernfrage: Verschwindet ||P_perp k_lambda|| fuer lambda -> inf bei angemessenem N?

Parameter-Grid (cluster/dim ~ 25-30%):
  lambda=3,  N=30   (dim=31,  cluster~5,   ~16%)
  lambda=5,  N=55   (dim=56,  cluster~16,  ~29%)
  lambda=7,  N=85   (dim=86,  cluster~25,  ~29%)
  lambda=10, N=120  (dim=121, cluster~35,  ~29%)
  lambda=15, N=180  (dim=181, cluster~55,  ~30%)

Geschaetzte Laufzeit: ~3 Stunden (ellmos-services, CCX13, 2 vCPU)

Autor: LG (Opus 4.7 [1M], Session 22)
Datum: 2026-04-19
"""

import os
import sys
import time
import json
import mpmath as mp
import numpy as np
from datetime import datetime

DPS = 50

CONFIGS = [
    {"lam": 3.0,  "N": 30},
    {"lam": 5.0,  "N": 55},
    {"lam": 7.0,  "N": 85},
    {"lam": 10.0, "N": 120},
    {"lam": 15.0, "N": 180},
]

RIEMANN_ZEROS = [
    14.134725141734693,
    21.022039638771554,
    25.010857580145688,
    30.424876125859513,
    32.935061587739189,
]

# =================================================================
# Matrix-Bausteine (aus dirichlet_ccm_fourier_mp.py, selbststaendig)
# =================================================================

def rho_even_mp(x):
    return mp.exp(x / 2) / (mp.exp(x) - mp.exp(-x))

def alpha_L_mp(n, L):
    if n == 0:
        return mp.mpf(0)
    def integrand(x):
        return mp.sin(2 * mp.pi * n * x / L) * rho_even_mp(x)
    return mp.quad(integrand, [0, L]) / mp.pi

def beta_L_mp(n, L):
    def integrand(x):
        return x * mp.cos(2 * mp.pi * n * x / L) * rho_even_mp(x)
    return mp.quad(integrand, [0, L]) / L

def gamma_L_mp(n, L, c_plus_w):
    def integrand(x):
        return (mp.cos(2 * mp.pi * n * x / L) - mp.exp(-x / 2)) * rho_even_mp(x)
    return mp.quad(integrand, [0, L]) + c_plus_w

def compute_c_plus_w_mp(L):
    eL2 = mp.exp(L / 2)
    term1 = mp.mpf("0.5") * mp.log((eL2 - 1) / (eL2 + 1))
    term2 = mp.atan(eL2)
    term3 = -mp.pi / 4
    term4 = mp.euler / 2
    term5 = mp.mpf("0.5") * mp.log(8 * mp.pi)
    return term1 + term2 + term3 + term4 + term5

def build_WR_matrix_mp(N, L_mp):
    c_plus_w = compute_c_plus_w_mp(L_mp)
    alpha_vals = {}
    beta_vals = {}
    gamma_vals = {}
    for n in range(-N, N + 1):
        alpha_vals[n] = alpha_L_mp(n, L_mp)
        beta_vals[n] = beta_L_mp(n, L_mp)
        gamma_vals[n] = gamma_L_mp(n, L_mp, c_plus_w)
    size = 2 * N + 1
    WR = mp.matrix(size, size)
    for i in range(size):
        n = i - N
        for j in range(size):
            m = j - N
            if n == m:
                WR[i, j] = 2 * gamma_vals[n] - 2 * beta_vals[n]
            else:
                WR[i, j] = (alpha_vals[m] - alpha_vals[n]) / (n - m)
    return WR

def build_W02_matrix_mp(N, L_mp):
    size = 2 * N + 1
    W02 = mp.matrix(size, size)
    pre = 32 * L_mp * mp.sinh(L_mp / 4) ** 2
    L2 = L_mp * L_mp
    for i in range(size):
        n = i - N
        for j in range(size):
            m = j - N
            num = L2 - 16 * mp.pi ** 2 * m * n
            den = (L2 + 16 * mp.pi ** 2 * m ** 2) * (L2 + 16 * mp.pi ** 2 * n ** 2)
            W02[i, j] = pre * num / den
    return W02

def is_prime_power(k):
    if k < 2:
        return None, 0
    for p in range(2, int(np.sqrt(k)) + 1):
        if k % p == 0:
            m = 0
            r = k
            while r % p == 0:
                r //= p
                m += 1
            if r == 1:
                return p, m
            else:
                return None, 0
    return k, 1

def mangoldt_mp(k):
    p, m = is_prime_power(k)
    return mp.log(mp.mpf(p)) if p else mp.mpf(0)

def q_kernel_mp(m, n, y, L_mp):
    if m == n:
        return 2 * (1 - y / L_mp) * mp.cos(2 * mp.pi * n * y / L_mp)
    else:
        return (mp.sin(2 * mp.pi * m * y / L_mp) - mp.sin(2 * mp.pi * n * y / L_mp)) / (mp.pi * (n - m))

def build_Wprime_matrix_mp(N, L_mp):
    size = 2 * N + 1
    Wp = mp.matrix(size, size)
    k_max = int(mp.exp(L_mp))
    for k in range(2, k_max + 1):
        lam_k = mangoldt_mp(k)
        if lam_k == 0:
            continue
        y = mp.log(mp.mpf(k))
        sqrt_k = mp.sqrt(mp.mpf(k))
        pre = lam_k / sqrt_k
        for i in range(size):
            n = i - N
            for j in range(size):
                m = j - N
                Wp[i, j] += pre * q_kernel_mp(m, n, y, L_mp)
    return Wp

def build_QW_mp(N, lam):
    L_mp = 2 * mp.log(mp.mpf(lam))
    print(f"  Baue QW: lambda={lam}, N={N}, L={float(L_mp):.4f}", flush=True)
    t0 = time.time()
    print(f"    W_R...", flush=True)
    WR = build_WR_matrix_mp(N, L_mp)
    print(f"    W_02...", flush=True)
    W02 = build_W02_matrix_mp(N, L_mp)
    print(f"    W_prime (k_max={int(mp.exp(L_mp))})...", flush=True)
    Wp = build_Wprime_matrix_mp(N, L_mp)
    size = 2 * N + 1
    M = mp.matrix(size, size)
    for i in range(size):
        for j in range(size):
            M[i, j] = -WR[i, j] - Wp[i, j] + W02[i, j]
    for i in range(size):
        for j in range(i + 1, size):
            avg = (M[i, j] + M[j, i]) / 2
            M[i, j] = avg
            M[j, i] = avg
    dt = time.time() - t0
    print(f"    QW fertig in {dt:.1f}s", flush=True)
    return M, L_mp

def project_to_even(M, N):
    size = 2 * N + 1
    dim = N + 1
    U = mp.matrix(size, dim)
    U[N, 0] = mp.mpf(1)
    inv_sqrt2 = 1 / mp.sqrt(mp.mpf(2))
    for k in range(1, N + 1):
        U[N + k, k] = inv_sqrt2
        U[N - k, k] = inv_sqrt2
    M_sub = U.T * M * U
    return M_sub, dim

def diagonalize_mp(M):
    n = M.rows
    print(f"    Diagonalisiere {n}x{n}...", flush=True)
    t0 = time.time()
    E, Q = mp.eigsy(M)
    dt = time.time() - t0
    print(f"    eigsy fertig in {dt:.1f}s", flush=True)
    w = [E[i] for i in range(n)]
    order = sorted(range(n), key=lambda i: float(w[i]))
    w_sorted = [w[i] for i in order]
    Q_sorted = mp.matrix(n, n)
    for new_idx, old_idx in enumerate(order):
        for row in range(n):
            Q_sorted[row, new_idx] = Q[row, old_idx]
    return w_sorted, Q_sorted

# =================================================================
# Hermite + Epstein (k_lambda)
# =================================================================

def hermite_h0(t):
    return mp.power(mp.pi, mp.mpf("-0.25")) * mp.exp(-t**2 / 2)

def hermite_h4(t):
    norm = 1 / mp.sqrt(384 * mp.sqrt(mp.pi))
    poly = 16 * t**4 - 48 * t**2 + 12
    return norm * poly * mp.exp(-t**2 / 2)

def h_educated_guess(t):
    c4 = mp.sqrt(3) / mp.power(2, mp.mpf("11") / 4)
    c0 = mp.mpf(3) / mp.power(2, mp.mpf("17") / 4)
    return c4 * hermite_h4(t) - c0 * hermite_h0(t)

def k_lambda_value(u, lam, n_terms=500):
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

def project_k_onto_even_basis(lam, N_val):
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    inv_sqrt_L = 1 / mp.sqrt(L_mp)
    dim = N_val + 1

    def k_of_x(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(u, lam)

    coeffs = mp.matrix(dim, 1)
    print(f"    Fourier-Projektion (dim={dim})...", flush=True)
    t0 = time.time()

    def integrand_0(x):
        return k_of_x(x) * inv_sqrt_L
    coeffs[0, 0] = mp.quad(integrand_0, [0, L_mp], maxdegree=8)

    for m in range(1, dim):
        def integrand_m(x, m_val=m):
            return k_of_x(x) * mp.sqrt(mp.mpf(2)) * inv_sqrt_L * mp.cos(2 * mp.pi * m_val * x / L_mp)
        coeffs[m, 0] = mp.quad(integrand_m, [0, L_mp], maxdegree=8)
        if m % 20 == 0:
            elapsed = time.time() - t0
            print(f"      m={m}/{dim-1} ({elapsed:.0f}s)", flush=True)

    dt = time.time() - t0
    print(f"    Fourier-Projektion fertig in {dt:.1f}s", flush=True)
    return coeffs, L_mp

# =================================================================
# Analyse-Hilfsfunktionen
# =================================================================

def inner_product(a, b, dim):
    return sum(a[i, 0] * b[i, 0] for i in range(dim))

def vec_norm(a, dim):
    return mp.sqrt(inner_product(a, a, dim))

def evaluate_F_even(coeffs, N_val, L_mp, gamma):
    g = mp.mpf(gamma)
    F = coeffs[0, 0] / g
    for k in range(1, N_val + 1):
        pole = 2 * mp.pi * k / L_mp
        F += coeffs[k, 0] * mp.sqrt(mp.mpf(2)) * g / (g**2 - pole**2)
    return F

# =================================================================
# Hauptanalyse pro (lambda, N)
# =================================================================

def run_config(lam, N_val):
    mp.mp.dps = DPS
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    dim = N_val + 1

    print(f"\n{'='*80}", flush=True)
    print(f"  LAMBDA = {lam}, N = {N_val}, dim = {dim}", flush=True)
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print(f"{'='*80}", flush=True)

    t_total = time.time()

    # QW
    M_full, _ = build_QW_mp(N_val, lam)
    M_even, dim_even = project_to_even(M_full, N_val)
    assert dim_even == dim
    w_sorted, Q_sorted = diagonalize_mp(M_even)

    w_min = w_sorted[0]

    # Cluster
    cluster_tol = mp.mpf("1e-10")
    cluster_idx = [k for k in range(dim) if abs(w_sorted[k] - w_min) < cluster_tol]
    cluster_size = len(cluster_idx)
    gap = w_sorted[cluster_size] - w_min if cluster_size < dim else mp.mpf(0)

    print(f"  w_min={mp.nstr(w_min, 8)}, Cluster={cluster_size}/{dim} ({100*cluster_size/dim:.1f}%), Gap={mp.nstr(gap, 3)}", flush=True)

    # k_lambda
    k_coeffs, _ = project_k_onto_even_basis(lam, N_val)
    k_norm_val = vec_norm(k_coeffs, dim)

    if k_norm_val < mp.mpf("1e-30"):
        print("  WARNUNG: k_lambda ~ 0", flush=True)
        return None

    k_normed = mp.matrix(dim, 1)
    for i in range(dim):
        k_normed[i, 0] = k_coeffs[i, 0] / k_norm_val

    # Cluster-Projektion
    cluster_proj_sq = mp.mpf(0)
    k_cluster = mp.matrix(dim, 1)
    best_ov = mp.mpf(0)
    best_k = 0
    for k in cluster_idx:
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]
        ov = inner_product(k_normed, qk, dim)
        cluster_proj_sq += ov**2
        if abs(ov) > best_ov:
            best_ov = abs(ov)
            best_k = k
        for i in range(dim):
            k_cluster[i, 0] += ov * qk[i, 0]

    cluster_proj = mp.sqrt(cluster_proj_sq)
    leakage = 1 - float(cluster_proj_sq)

    k_perp = mp.matrix(dim, 1)
    for i in range(dim):
        k_perp[i, 0] = k_normed[i, 0] - k_cluster[i, 0]
    perp_norm = vec_norm(k_perp, dim)

    print(f"  ||P_cl k||={mp.nstr(cluster_proj, 12)}, Leckage={leakage:.6e}, ||k_perp||={float(perp_norm):.6e}", flush=True)
    print(f"  Best single EV: k={best_k}, |ov|={mp.nstr(best_ov, 6)}", flush=True)

    # F-Zerlegung
    k_cl_coeffs = mp.matrix(dim, 1)
    k_perp_coeffs = mp.matrix(dim, 1)
    for i in range(dim):
        k_cl_coeffs[i, 0] = k_cluster[i, 0] * k_norm_val
        k_perp_coeffs[i, 0] = k_perp[i, 0] * k_norm_val

    # Cluster-Edge F
    cluster_edge = cluster_idx[-1]
    xi_edge = mp.matrix(dim, 1)
    for i in range(dim):
        xi_edge[i, 0] = Q_sorted[i, cluster_edge]

    print(f"\n  F(gamma) ZERLEGUNG:", flush=True)
    print(f"  {'gamma':>12s}  {'|F_raw|':>12s}  {'|F_cl|':>12s}  {'|F_perp|':>12s}  {'|F_xi_edge|':>12s}", flush=True)

    f_data = []
    for gamma in RIEMANN_ZEROS:
        F_raw = evaluate_F_even(k_coeffs, N_val, L_mp, gamma)
        F_cl = evaluate_F_even(k_cl_coeffs, N_val, L_mp, gamma)
        F_perp = evaluate_F_even(k_perp_coeffs, N_val, L_mp, gamma)
        F_xi = evaluate_F_even(xi_edge, N_val, L_mp, gamma)
        row = {
            "gamma": gamma,
            "F_raw": abs(float(F_raw)),
            "F_cl": abs(float(F_cl)),
            "F_perp": abs(float(F_perp)),
            "F_xi_edge": abs(float(F_xi)),
        }
        f_data.append(row)
        print(f"  {gamma:12.6f}  {row['F_raw']:12.3e}  {row['F_cl']:12.3e}  {row['F_perp']:12.3e}  {row['F_xi_edge']:12.3e}", flush=True)

    # Rayleigh
    QWk = M_even * k_normed
    eps_rayleigh = inner_product(k_normed, QWk, dim)

    dt_total = time.time() - t_total
    print(f"\n  Rayleigh={mp.nstr(eps_rayleigh, 10)}, Diff={mp.nstr(eps_rayleigh - w_min, 4)}", flush=True)
    print(f"  Gesamtzeit: {dt_total:.1f}s ({dt_total/60:.1f} min)", flush=True)

    return {
        "lambda": lam,
        "N": N_val,
        "dim": dim,
        "w_min": float(w_min),
        "cluster_size": cluster_size,
        "cluster_frac": cluster_size / dim,
        "gap": float(gap),
        "leakage": leakage,
        "perp_norm": float(perp_norm),
        "best_overlap": float(best_ov),
        "best_k": best_k,
        "rayleigh": float(eps_rayleigh),
        "k_norm": float(k_norm_val),
        "F_zeros": f_data,
        "time_s": dt_total,
    }


def main():
    mp.mp.dps = DPS
    print("=" * 80, flush=True)
    print("C2/MS2 SERVER-LAUF: LECKAGE-SKALIERUNG MIT N ~ lambda^2", flush=True)
    print(f"  Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print(f"  DPS = {DPS}", flush=True)
    print("=" * 80, flush=True)

    for cfg in CONFIGS:
        print(f"\n  Config: lambda={cfg['lam']}, N={cfg['N']} (dim={cfg['N']+1})", flush=True)

    results = []
    for cfg in CONFIGS:
        r = run_config(cfg["lam"], cfg["N"])
        if r:
            results.append(r)
            # Zwischenspeichern
            with open("/root/ccm_analysis/c2_leakage_results.json", "w") as f:
                json.dump(results, f, indent=2, default=str)

    # Zusammenfassung
    print(f"\n{'='*80}", flush=True)
    print(f"  LECKAGE-SKALIERUNG ZUSAMMENFASSUNG (N ~ lambda^2)", flush=True)
    print(f"  Ende: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"  {'lam':>6s}  {'N':>5s}  {'Cl/Dim':>10s}  {'Gap':>10s}  {'Leckage':>12s}  {'||k_perp||':>12s}  {'|F_raw(g1)|':>12s}  {'|F_cl(g1)|':>12s}  {'|F_perp(g1)|':>12s}  {'Zeit':>8s}", flush=True)
    for r in results:
        cl_dim = f"{r['cluster_size']}/{r['dim']}"
        f_raw = r["F_zeros"][0]["F_raw"] if r["F_zeros"] else 0
        f_cl = r["F_zeros"][0]["F_cl"] if r["F_zeros"] else 0
        f_perp = r["F_zeros"][0]["F_perp"] if r["F_zeros"] else 0
        print(f"  {r['lambda']:6.1f}  {r['N']:5d}  {cl_dim:>10s}  {r['gap']:10.2e}  {r['leakage']:12.4e}  {r['perp_norm']:12.4e}  {f_raw:12.3e}  {f_cl:12.3e}  {f_perp:12.3e}  {r['time_s']/60:6.1f}m", flush=True)

    # Skalierungs-Fit
    if len(results) >= 3:
        lams = [r["lambda"] for r in results]
        leaks = [r["leakage"] for r in results]
        perps = [r["perp_norm"] for r in results]

        valid_l = [(l, e) for l, e in zip(lams, leaks) if e > 0]
        valid_p = [(l, p) for l, p in zip(lams, perps) if p > 0]

        if len(valid_l) >= 2:
            log_lams = [np.log(l) for l, _ in valid_l]
            log_leaks = [np.log(e) for _, e in valid_l]
            n_pts = len(log_lams)
            mean_x = sum(log_lams) / n_pts
            mean_y = sum(log_leaks) / n_pts
            ss_xx = sum((x - mean_x)**2 for x in log_lams)
            ss_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(log_lams, log_leaks))
            alpha_leak = ss_xy / ss_xx if ss_xx > 0 else 0
            print(f"\n  Skalierungs-Exponent Leckage: alpha = {alpha_leak:.2f} (Leckage ~ lam^alpha)", flush=True)

        if len(valid_p) >= 2:
            log_lams = [np.log(l) for l, _ in valid_p]
            log_perps = [np.log(p) for _, p in valid_p]
            n_pts = len(log_lams)
            mean_x = sum(log_lams) / n_pts
            mean_y = sum(log_perps) / n_pts
            ss_xx = sum((x - mean_x)**2 for x in log_lams)
            ss_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(log_lams, log_perps))
            alpha_perp = ss_xy / ss_xx if ss_xx > 0 else 0
            print(f"  Skalierungs-Exponent ||k_perp||: alpha = {alpha_perp:.2f} (||k_perp|| ~ lam^alpha)", flush=True)

        print(f"  (alpha < 0 => verschwindet fuer lambda -> inf)", flush=True)

    print(f"\n  FERTIG.", flush=True)


if __name__ == "__main__":
    main()
