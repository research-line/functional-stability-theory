"""
c2_poisson_decomposition.py — AP2: Poisson-Zerlegung k_lambda = k_n1 + r_n2

Zerlegung des CCM Educated Guess:
  k_lambda(u) = u^{1/2} lambda^{1/2} sum_{n>=1} h(lambda*n*u)

in Hauptterm (n=1) und arithmetischen Tail (n>=2):
  k_n1(u) = u^{1/2} lambda^{1/2} h(lambda*u)          <-- glatter Kern
  r_n2(u) = u^{1/2} lambda^{1/2} sum_{n>=2} h(lambda*n*u) <-- arithm. Korrektur

Poisson-Verifikation:
  sum h(lambda*n*u) = (1/(lambda*u)) sum h_hat(m/(lambda*u))
  (mit h(0)=h_hat(0)=0 per Konstruktion)

h_hat analytisch:
  h_hat(xi) = c4 * h4_hat(xi) - c0 * h0_hat(xi)
  h0_hat(xi) = pi^{-1/4} sqrt(2pi) exp(-2pi^2 xi^2)
  h4_hat(xi) = (384 sqrt(pi))^{-1/2} sqrt(2pi) exp(-2pi^2 xi^2)
               * (12 - 192 pi^2 xi^2 + 256 pi^4 xi^4)

Messungen: L2-Norm, Cluster-Projektion, QW-Residual, F(gamma), lambda-Skalierung

Autor: LG (Opus 4.6, Session 23)
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
)
from c2_approximation_test import (
    k_lambda_value,
    h_educated_guess,
    inner_product,
    norm,
    evaluate_F_even,
)

DPS = int(os.environ.get("DPS", 50))

CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
]

GAMMA_RH = [14.134725141734693, 21.022039638771554, 25.010857580145688]


def h_hat_analytical(xi):
    """Analytische Fouriertransformierte h_hat(xi) = int h(t) e^{-2pi i xi t} dt.

    Herleitung: h = c4*h4 - c0*h0 mit Hermite-Funktionen, alle gerade.
    FT von P(t)*exp(-t^2/2) liefert Q(xi)*exp(-2pi^2 xi^2) analytisch.
    """
    xi = mp.mpf(xi)
    gauss = mp.exp(-2 * mp.pi**2 * xi**2)
    sqrt_2pi = mp.sqrt(2 * mp.pi)

    c4 = mp.sqrt(3) / mp.power(2, mp.mpf("11") / 4)
    c0 = mp.mpf(3) / mp.power(2, mp.mpf("17") / 4)

    h0_hat = mp.power(mp.pi, mp.mpf("-0.25")) * sqrt_2pi * gauss

    pi2xi2 = mp.pi**2 * xi**2
    poly = 12 - 192 * pi2xi2 + 256 * pi2xi2**2
    h4_hat = poly * sqrt_2pi * gauss / mp.sqrt(384 * mp.sqrt(mp.pi))

    return c4 * h4_hat - c0 * h0_hat


def project_to_fourier(f_of_x, L_mp, N_val):
    """Projiziere f(x) auf Even-Fourier-Basis auf [0, L]."""
    inv_sqrt_L = 1 / mp.sqrt(L_mp)
    dim = N_val + 1
    coeffs = mp.matrix(dim, 1)

    coeffs[0, 0] = mp.quad(lambda x: f_of_x(x) * inv_sqrt_L,
                           [0, L_mp], maxdegree=8)
    for m in range(1, dim):
        def integrand(x, mv=m):
            return (f_of_x(x) * mp.sqrt(mp.mpf(2)) * inv_sqrt_L
                    * mp.cos(2 * mp.pi * mv * x / L_mp))
        coeffs[m, 0] = mp.quad(integrand, [0, L_mp], maxdegree=8)
    return coeffs


def run_decomposition(lam, N_val):
    mp.mp.dps = DPS
    dim = N_val + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    print(f"\n{'='*80}")
    print(f"  POISSON-ZERLEGUNG: lam={lam}, N={N_val}, dim={dim}, L={float(L_mp):.4f}")
    print(f"{'='*80}")

    # --- Phase 0: QW Eigensystem ---
    t0 = time.time()
    M_full, _ = build_QW_mp(N_val, lam, chi_func=chi_trivial, q_mod=1,
                            include_W02=True, sign_WR=-1, verbose=False)
    M_even, _ = project_to_parity(M_full, N_val, parity="even")
    w_sorted, Q_sorted = diagonalize_mp(M_even, verbose=False)
    w_min = w_sorted[0]
    print(f"\n  QW in {time.time()-t0:.1f}s, w_min={mp.nstr(w_min, 8)}")

    # --- Phase 1: Analytische Kontrollen ---
    print(f"\n--- Phase 1: Analytische Kontrollen ---")
    h_at_0 = h_educated_guess(mp.mpf(0))
    hh_at_0 = h_hat_analytical(mp.mpf(0))
    print(f"  h(0)     = {mp.nstr(h_at_0, 6)}  (soll 0)")
    print(f"  h_hat(0) = {mp.nstr(hh_at_0, 6)}  (soll 0)")

    # h_hat numerisch vs analytisch bei xi=0.3
    def integrand_check(t):
        return h_educated_guess(t) * mp.cos(2 * mp.pi * mp.mpf("0.3") * t)
    hh_num = 2 * mp.quad(integrand_check, [0, 20], maxdegree=8)
    hh_ana = h_hat_analytical(mp.mpf("0.3"))
    print(f"  h_hat(0.3) num = {mp.nstr(hh_num, 8)}")
    print(f"  h_hat(0.3) ana = {mp.nstr(hh_ana, 8)}")
    print(f"  rel. Fehler     = {float(abs(hh_num - hh_ana) / abs(hh_ana)):.2e}")

    print(f"\n  h_hat-Profil:")
    for xi in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]:
        hh = h_hat_analytical(xi)
        print(f"    h_hat({xi:5.2f}) = {float(hh):14.6e}")

    # --- Phase 2: Poisson-Identitaet ---
    print(f"\n--- Phase 2: Poisson-Identitaet verifizieren ---")
    test_xs = [L_mp * mp.mpf(j) / 6 for j in range(1, 6)]
    max_poi_terms = 0
    for x_test in test_xs:
        u_test = mp.exp(x_test) / lam_mp
        k_dir = k_lambda_value(float(u_test), lam)

        prefac = 1 / (mp.sqrt(lam_mp) * mp.sqrt(u_test))
        s_pois = mp.mpf(0)
        m_used = 0
        for m in range(1, 500):
            xi = mp.mpf(m) / (lam_mp * u_test)
            hh = h_hat_analytical(xi)
            s_pois += hh
            m_used = m
            if m > 3 and abs(hh) < mp.mpf("1e-40"):
                break
        k_poi = prefac * s_pois
        max_poi_terms = max(max_poi_terms, m_used)

        rel = abs(float(k_dir - k_poi)) / max(abs(float(k_dir)), 1e-50)
        print(f"  u={float(u_test):.4f}: direct={float(k_dir):12.4e}  "
              f"poisson={float(k_poi):12.4e}  rel_err={rel:.2e}  m_terms={m_used}")

    print(f"  Max Poisson-Terme benoetigt: {max_poi_terms}")

    # --- Phase 3: n-Split Fourier-Projektion ---
    print(f"\n--- Phase 3: n-Split Fourier-Projektion ---")
    t0 = time.time()

    def k_full_of_x(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    def k_n1_of_x(x):
        ex = mp.exp(x)
        return mp.sqrt(ex) * h_educated_guess(ex)

    c_full = project_to_fourier(k_full_of_x, L_mp, N_val)
    c_n1 = project_to_fourier(k_n1_of_x, L_mp, N_val)

    c_n2 = mp.matrix(dim, 1)
    for i in range(dim):
        c_n2[i, 0] = c_full[i, 0] - c_n1[i, 0]

    norm_full = norm(c_full, dim)
    norm_n1 = norm(c_n1, dim)
    norm_n2 = norm(c_n2, dim)

    print(f"  Projektion in {time.time()-t0:.1f}s")
    print(f"  ||k_full|| = {float(norm_full):.8e}")
    print(f"  ||k_n1||   = {float(norm_n1):.8e}  ({float(norm_n1/norm_full)*100:.4f}%)")
    print(f"  ||r_n2||   = {float(norm_n2):.8e}  ({float(norm_n2/norm_full)*100:.6f}%)")

    # Normiere (alles durch ||k_full||)
    k_normed = mp.matrix(dim, 1)
    n1_part = mp.matrix(dim, 1)
    n2_part = mp.matrix(dim, 1)
    for i in range(dim):
        k_normed[i, 0] = c_full[i, 0] / norm_full
        n1_part[i, 0] = c_n1[i, 0] / norm_full
        n2_part[i, 0] = c_n2[i, 0] / norm_full

    # --- Phase 4: Eigenvektorbasis-Zerlegung ---
    print(f"\n--- Phase 4: Eigenvektorbasis-Zerlegung (erste 15) ---")
    print(f"  {'k':>4s}  {'c_full':>12s}  {'c_n1':>12s}  {'c_n2':>12s}  "
          f"{'w_k-w_min':>12s}")

    ov_full = []
    ov_n1 = []
    ov_n2 = []
    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Q_sorted[i, k]
        of = inner_product(k_normed, qk, dim)
        o1 = inner_product(n1_part, qk, dim)
        o2 = inner_product(n2_part, qk, dim)
        ov_full.append(of)
        ov_n1.append(o1)
        ov_n2.append(o2)
        if k < 15:
            print(f"  {k:4d}  {float(of):12.4e}  {float(o1):12.4e}  "
                  f"{float(o2):12.4e}  {float(w_sorted[k]-w_min):12.4e}")

    # --- Phase 5: Cluster-Analyse ---
    print(f"\n--- Phase 5: Cluster-Projektion ---")

    cluster_results = {}
    for delta_name, delta in [("eng", 1e-10), ("weit", 1.0)]:
        cl_idx = [k for k in range(dim)
                  if abs(float(w_sorted[k] - w_min)) < delta]
        non_cl = [k for k in range(dim) if k not in cl_idx]
        cl_size = len(cl_idx)

        proj_full_cl = sum(float(ov_full[k])**2 for k in cl_idx)
        proj_n1_cl = sum(float(ov_n1[k])**2 for k in cl_idx)
        proj_n2_cl = sum(float(ov_n2[k])**2 for k in cl_idx)
        cross_cl = sum(2 * float(ov_n1[k]) * float(ov_n2[k]) for k in cl_idx)

        proj_full_perp = sum(float(ov_full[k])**2 for k in non_cl)
        proj_n1_perp = sum(float(ov_n1[k])**2 for k in non_cl)
        proj_n2_perp = sum(float(ov_n2[k])**2 for k in non_cl)
        cross_perp = sum(2 * float(ov_n1[k]) * float(ov_n2[k]) for k in non_cl)

        perp_full = np.sqrt(max(0, proj_full_perp))

        print(f"\n  Cluster '{delta_name}' (delta={delta:.0e}, size={cl_size}):")
        print(f"    ||P_cl k||^2     = {proj_full_cl:.10e}")
        print(f"    ||P_perp k||     = {perp_full:.6e}")
        print(f"    ||P_cl k_n1||^2  = {proj_n1_cl:.10e}")
        print(f"    ||P_cl r_n2||^2  = {proj_n2_cl:.10e}")
        print(f"    Cross-Term cl    = {cross_cl:.10e}")
        print(f"    ||P_perp k_n1||^2= {proj_n1_perp:.10e}")
        print(f"    ||P_perp r_n2||^2= {proj_n2_perp:.10e}")
        print(f"    Cross-Term perp  = {cross_perp:.10e}")

        if proj_full_perp > 1e-30:
            n1_frac = proj_n1_perp / proj_full_perp * 100
            n2_frac = proj_n2_perp / proj_full_perp * 100
            cross_frac = cross_perp / proj_full_perp * 100
            print(f"    -> Leckage-Zerlegung: n=1 {n1_frac:.1f}%, "
                  f"n>=2 {n2_frac:.1f}%, Cross {cross_frac:.1f}%")

        cluster_results[delta_name] = {
            "size": cl_size,
            "perp_full": perp_full,
            "perp_n1_sq": proj_n1_perp,
            "perp_n2_sq": proj_n2_perp,
        }

    # --- Phase 6: QW-Residual ---
    print(f"\n--- Phase 6: QW-Residual ---")

    QWk = M_even * k_normed
    rho = inner_product(k_normed, QWk, dim)

    QWn1 = M_even * n1_part
    QWn2 = M_even * n2_part
    res_full = mp.matrix(dim, 1)
    res_n1 = mp.matrix(dim, 1)
    res_n2 = mp.matrix(dim, 1)
    for i in range(dim):
        res_full[i, 0] = QWk[i, 0] - rho * k_normed[i, 0]
        res_n1[i, 0] = QWn1[i, 0] - rho * n1_part[i, 0]
        res_n2[i, 0] = QWn2[i, 0] - rho * n2_part[i, 0]

    r_full = float(norm(res_full, dim))
    r_n1 = float(norm(res_n1, dim))
    r_n2 = float(norm(res_n2, dim))

    print(f"  rho = {mp.nstr(rho, 10)}")
    print(f"  ||(QW-rho)k_full|| = {r_full:.6e}")
    print(f"  ||(QW-rho)k_n1||   = {r_n1:.6e}")
    print(f"  ||(QW-rho)r_n2||   = {r_n2:.6e}")

    print(f"\n  Spektrale Residual-Zerlegung |c_k|^2 (w_k-rho)^2:")
    print(f"  {'k':>4s}  {'full':>14s}  {'n1':>14s}  {'n2':>14s}")
    for k in range(min(15, dim)):
        dw = float(w_sorted[k] - rho)
        cf = float(ov_full[k])**2 * dw**2
        c1 = float(ov_n1[k])**2 * dw**2
        c2 = float(ov_n2[k])**2 * dw**2
        if cf > 1e-25 or k < 10:
            print(f"  {k:4d}  {cf:14.4e}  {c1:14.4e}  {c2:14.4e}")

    # --- Phase 7: F(gamma) ---
    print(f"\n--- Phase 7: F(gamma) an Riemann-Nullstellen ---")
    print(f"  {'gamma':>10s}  {'F_full':>12s}  {'F_n1':>12s}  {'F_n2':>12s}  "
          f"{'|F_n2/F_full|':>14s}")

    f_n2_vals = []
    for g in GAMMA_RH:
        F_full = evaluate_F_even(c_full, N_val, L_mp, g)
        F_n1 = evaluate_F_even(c_n1, N_val, L_mp, g)
        F_n2 = evaluate_F_even(c_n2, N_val, L_mp, g)
        ratio = abs(float(F_n2)) / max(abs(float(F_full)), 1e-50)
        print(f"  {g:10.6f}  {float(F_full):12.4e}  {float(F_n1):12.4e}  "
              f"{float(F_n2):12.4e}  {ratio:14.4e}")
        f_n2_vals.append(abs(float(F_n2)))

    return {
        "lam": lam, "N": N_val,
        "norm_full": float(norm_full),
        "norm_n1": float(norm_n1),
        "norm_n2": float(norm_n2),
        "n2_frac": float(norm_n2 / norm_full),
        "perp_full_eng": cluster_results.get("eng", {}).get("perp_full", 0),
        "perp_n2_sq_eng": cluster_results.get("eng", {}).get("perp_n2_sq", 0),
        "perp_full_weit": cluster_results.get("weit", {}).get("perp_full", 0),
        "perp_n2_sq_weit": cluster_results.get("weit", {}).get("perp_n2_sq", 0),
        "res_full": r_full,
        "res_n1": r_n1,
        "res_n2": r_n2,
        "F_n2_gamma1": f_n2_vals[0] if f_n2_vals else 0,
    }


def main():
    mp.mp.dps = DPS
    print("=" * 80)
    print("AP2: POISSON-ZERLEGUNG k_lambda = k_n1 + r_n2")
    print("=" * 80)
    print(f"  k_n1(u) = u^(1/2) lam^(1/2) h(lam*u)            [n=1 Hauptterm]")
    print(f"  r_n2(u) = u^(1/2) lam^(1/2) sum_{{n>=2}} h(lam*n*u) [arithm. Tail]")
    print(f"  Poisson: sum h(lam*n*u) = (1/(lam*u)) sum h_hat(m/(lam*u))")

    all_results = []
    for cfg in CONFIGS:
        res = run_decomposition(cfg["lam"], cfg["N"])
        all_results.append(res)

    if len(all_results) >= 2:
        print(f"\n{'='*80}")
        print(f"  LAMBDA-SKALIERUNG")
        print(f"{'='*80}")

        lams = np.array([r["lam"] for r in all_results])
        log_lams = np.log(lams)

        metrics = [
            ("||r_n2||/||k||", [r["n2_frac"] for r in all_results]),
            ("||P_perp k|| eng", [r["perp_full_eng"] for r in all_results]),
            ("P_perp_n2^2 eng", [r["perp_n2_sq_eng"] for r in all_results]),
            ("||P_perp k|| weit", [r["perp_full_weit"] for r in all_results]),
            ("P_perp_n2^2 weit", [r["perp_n2_sq_weit"] for r in all_results]),
            ("||(QW-rho)k_full||", [r["res_full"] for r in all_results]),
            ("||(QW-rho)r_n2||", [r["res_n2"] for r in all_results]),
            ("|F_n2(gamma_1)|", [r["F_n2_gamma1"] for r in all_results]),
        ]

        print(f"\n  {'Metrik':>25s}", end="")
        for r in all_results:
            print(f"  {'lam='+str(r['lam']):>12s}", end="")
        print(f"  {'Exponent':>10s}")

        for name, vals in metrics:
            print(f"  {name:>25s}", end="")
            for v in vals:
                print(f"  {v:12.4e}", end="")
            if all(v > 0 for v in vals) and len(vals) >= 2:
                log_vals = np.log(vals)
                exp = (log_vals[1] - log_vals[0]) / (log_lams[1] - log_lams[0])
                print(f"  {exp:10.2f}", end="")
            else:
                print(f"  {'N/A':>10s}", end="")
            print()

        print(f"\n  Vergleich: Server-Leckage ||P_perp k|| ~ lam^{{-1.84}}")
        print(f"  Entscheidungsregel:")
        print(f"    Falls ||r_n2|| denselben Trend traegt -> AP2 ist der analytische Kern")
        print(f"    Falls nicht -> Leckage kommt aus Cluster-Geometrie, nicht Poisson-Tail")


if __name__ == "__main__":
    main()
