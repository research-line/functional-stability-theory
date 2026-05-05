"""
c2bt_server_quasimode_scaling.py

SERVER-LAUF: Quasimode-Skalierung bei N proportional zu L = 2 log(lambda).

Misst das vollstaendige Quasimode-Paket:
  - lambda, N, L, n_cl
  - u_0, u_{n_cl}, g_cl = u_{n_cl} - u_0
  - Rayleigh-Fehler rho_k - u_0
  - Residuen R_0 = ||(A-u0)k||, R_rho = ||(A-rho)k||
  - Direkte Leckage ||(I-P_cl)k||
  - Quotientenprobe R_0/g_cl, R_rho/g_cl
  - mu/alpha

Skalierung: N = max(40, ceil(N_FACTOR * L))
  - Sicherstellt, dass n_cl < dim/2 fuer alle lambda

Laeuft auf ellmos-services (CCX13, 2vCPU, 8GB).
"""

import os, sys, time, json
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import mpmath as mp
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, chi_trivial,
)
from c2_approximation_test import (
    k_lambda_value, norm,
)
from c2_poisson_decomposition import (
    project_to_fourier,
)

DPS = int(os.environ.get("DPS", 35))
N_FACTOR = int(os.environ.get("N_FACTOR", 12))
CLUSTER_TOL = float(os.environ.get("CLUSTER_TOL", 1e-6))

LAMBDAS = [3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0,
           20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0]

OUTPUT_JSON = "c2bt_quasimode_scaling_results.json"


def to_np(M, dim):
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = float(M[i, j])
    return A


def compute_quasimode_full(lam, verbose=True):
    mp.mp.dps = DPS
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

    N = max(40, int(np.ceil(N_FACTOR * L)))
    dim = N + 1

    if verbose:
        print(f"\n  lam={lam:.1f}, L={L:.3f}, N={N}, dim={dim}")

    # Build operators
    t0 = time.time()
    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")
    dt_build = time.time() - t0

    A_np = to_np(Aq, dim)
    H_np = to_np(Ah, dim)
    W_np = A_np - H_np

    # ut (rank-1 direction)
    col0 = W_np[:, 0]
    cn = np.linalg.norm(col0)
    ut = col0 / cn

    # Full eigendecomposition of A
    t0 = time.time()
    ws, vs = np.linalg.eigh(A_np)
    dt_eig = time.time() - t0
    w0 = ws[0]

    # Cluster identification (from A)
    cl_mask = np.abs(ws - w0) < CLUSTER_TOL
    n_cl = int(np.sum(cl_mask))
    cl_idx = np.where(cl_mask)[0]
    off_idx = np.where(~cl_mask)[0]

    # Cluster gap
    if len(off_idx) > 0:
        gap = float(ws[off_idx[0]] - w0)
        u_first_off = float(ws[off_idx[0]])
    else:
        gap = 0.0
        u_first_off = w0

    # Also identify gap from T = P_perp H P_perp (alternative)
    P_perp = np.eye(dim) - np.outer(ut, ut)
    T_np = P_perp @ H_np @ P_perp
    ts_all = np.linalg.eigvalsh(T_np)
    ts_nontrivial = ts_all[np.abs(ts_all) > 1e-10]
    if len(ts_nontrivial) > 0:
        t_min = float(np.min(ts_nontrivial))
        t_gap_from_w0 = t_min - w0
    else:
        t_min = w0
        t_gap_from_w0 = 0.0

    # Poisson kernel k
    t0 = time.time()
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
    dt_k = time.time() - t0

    alpha = float(ut @ kn)

    # Rayleigh quotient
    Ak = A_np @ kn
    rho_k = float(kn @ Ak)
    rayleigh_err = rho_k - w0

    # Quasimode residuals
    res_w0 = Ak - w0 * kn
    res_rho = Ak - rho_k * kn
    R_0 = float(np.linalg.norm(res_w0))
    R_rho = float(np.linalg.norm(res_rho))

    # H-residual for comparison
    Hk = H_np @ kn
    R_H = float(np.linalg.norm(Hk - w0 * kn))

    # Cluster projection and leakage
    V_cl = vs[:, cl_idx]
    k_cl = V_cl @ (V_cl.T @ kn)
    overlap_cl = float(np.dot(kn, k_cl))
    leakage = max(0.0, 1.0 - overlap_cl)
    k_perp_norm = np.sqrt(leakage)

    # Overlap with single ground state
    overlap_q0 = float(kn @ vs[:, 0])**2

    # mu/alpha
    mu_alpha = float(ut @ res_w0) / alpha

    # Quotientenproben
    R0_over_gap = R_0 / gap if gap > 1e-30 else float('inf')
    Rrho_over_gap = R_rho / gap if gap > 1e-30 else float('inf')

    result = {
        'lam': lam, 'L': L, 'N': N, 'dim': dim,
        'n_cl': n_cl, 'n_off': dim - n_cl,
        'w0': w0, 'u_first_off': u_first_off,
        'gap': gap, 't_min': t_min, 't_gap_from_w0': t_gap_from_w0,
        'rho_k': rho_k, 'rayleigh_err': rayleigh_err,
        'R_0': R_0, 'R_rho': R_rho, 'R_H': R_H,
        'overlap_q0': overlap_q0, 'overlap_cl': overlap_cl,
        'leakage': leakage, 'k_perp_norm': k_perp_norm,
        'R0_over_gap': R0_over_gap, 'Rrho_over_gap': Rrho_over_gap,
        'alpha': alpha, 'mu_alpha': mu_alpha, 'C_tilde': cn,
        'dt_build': dt_build, 'dt_eig': dt_eig, 'dt_k': dt_k,
    }

    if verbose:
        print(f"    Build: {dt_build:.0f}s, Eig: {dt_eig:.1f}s, k: {dt_k:.0f}s")
        print(f"    n_cl={n_cl}, n_off={dim-n_cl}, gap={gap:.4e}")
        print(f"    R_0={R_0:.4e}, R_rho={R_rho:.4e}, R_H={R_H:.2f}")
        print(f"    ||k_perp||={k_perp_norm:.4e}, R_0/gap={R0_over_gap:.4e}")
        print(f"    mu/alpha={mu_alpha:.4e}")

    return result


def power_law_fit(x, y):
    mask = (np.array(x) > 0) & (np.array(y) > 0)
    if np.sum(mask) < 3:
        return 0, 0, 0
    lx = np.log(np.array(x)[mask])
    ly = np.log(np.array(y)[mask])
    n = len(lx)
    sx, sy = np.sum(lx), np.sum(ly)
    sxx, sxy = np.sum(lx**2), np.sum(lx * ly)
    denom = n * sxx - sx**2
    if abs(denom) < 1e-30:
        return 0, 0, 0
    alpha = (n * sxy - sx * sy) / denom
    logC = (sy - alpha * sx) / n
    ss_res = np.sum((ly - alpha * lx - logC)**2)
    ss_tot = np.sum((ly - np.mean(ly))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 1e-30 else 0
    return alpha, np.exp(logC), r2


def main():
    print("C2bt SERVER: Quasimode-Skalierung (N proportional L)")
    print(f"DPS={DPS}, N_FACTOR={N_FACTOR}, CLUSTER_TOL={CLUSTER_TOL}")
    print(f"Lambda-Range: {LAMBDAS[0]} .. {LAMBDAS[-1]} ({len(LAMBDAS)} Punkte)")
    print()

    results = []
    t_total = time.time()

    for lam in LAMBDAS:
        try:
            d = compute_quasimode_full(lam)
            results.append(d)
        except Exception as e:
            print(f"  FEHLER bei lam={lam}: {e}")
            import traceback
            traceback.print_exc()

    dt_total = time.time() - t_total
    print(f"\n  Gesamtzeit: {dt_total:.0f}s ({dt_total/60:.1f} min)")

    # Save JSON
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Ergebnisse gespeichert: {OUTPUT_JSON}")

    # Summary Table
    print()
    print("=" * 120)
    print("ZUSAMMENFASSUNG: QUASIMODE-SKALIERUNG (N prop. L)")
    print("=" * 120)
    print(f"{'lam':>5} {'N':>4} {'n_cl':>5} {'n_off':>5} {'gap':>12} "
          f"{'R_0':>12} {'R_rho':>12} {'||k_perp||':>12} "
          f"{'R0/gap':>10} {'Rrho/gap':>10} {'mu/alpha':>12}")
    for d in results:
        print(f"{d['lam']:5.1f} {d['N']:4d} {d['n_cl']:5d} {d['n_off']:5d} {d['gap']:12.4e} "
              f"{d['R_0']:12.4e} {d['R_rho']:12.4e} {d['k_perp_norm']:12.4e} "
              f"{d['R0_over_gap']:10.4e} {d['Rrho_over_gap']:10.4e} {d['mu_alpha']:12.4e}")

    # Power-law fits
    print()
    print("=" * 120)
    print("POWER-LAW FITS (y = C * lam^alpha)")
    print("=" * 120)

    lams = [d['lam'] for d in results]
    quantities = [
        ('R_0 = ||(A-w0)k||', [d['R_0'] for d in results]),
        ('R_rho = ||(A-rho)k||', [d['R_rho'] for d in results]),
        ('gap = g_cl', [d['gap'] for d in results]),
        ('||k_perp||', [d['k_perp_norm'] for d in results]),
        ('R_0/gap', [d['R0_over_gap'] for d in results]),
        ('R_rho/gap', [d['Rrho_over_gap'] for d in results]),
        ('|mu/alpha|', [abs(d['mu_alpha']) for d in results]),
        ('Rayleigh err', [abs(d['rayleigh_err']) for d in results]),
        ('R_H = ||(H-w0)k||', [d['R_H'] for d in results]),
    ]

    print(f"  {'Groesse':30s} {'Exponent':>10} {'C':>12} {'R^2':>8}")
    for name, vals in quantities:
        alpha_fit, C_fit, r2 = power_law_fit(lams, vals)
        print(f"  {name:30s} {alpha_fit:+10.4f} {C_fit:12.4e} {r2:8.4f}")

    # Key diagnostic
    print()
    print("=" * 120)
    print("SCHLUESSELDIAGNOSTIK: R_0/gap -> 0 ?")
    print("=" * 120)
    r0g = [d['R0_over_gap'] for d in results]
    alpha_fit, _, r2 = power_law_fit(lams, r0g)
    trend = "SINKT" if alpha_fit < -0.1 else "WAECHST" if alpha_fit > 0.1 else "FLACH"
    print(f"  R_0/gap ~ lam^({alpha_fit:+.3f}), R^2={r2:.3f} => {trend}")
    if alpha_fit < -0.1 and r2 > 0.5:
        print(f"  => QUASIMODE-ROUTE NUMERISCH BESTAETIGT")
        print(f"  => gap sinkt LANGSAMER als Residual => Quotient -> 0")
    elif alpha_fit > 0.1:
        print(f"  => WARNUNG: Quotient waechst! Gap sinkt SCHNELLER als Residual.")
    else:
        print(f"  => Unklarer Trend. Mehr Datenpunkte noetig.")

    kp = [d['k_perp_norm'] for d in results]
    alpha_fit_kp, _, r2_kp = power_law_fit(lams, kp)
    print(f"\n  ||k_perp|| ~ lam^({alpha_fit_kp:+.3f}), R^2={r2_kp:.3f}")
    if alpha_fit_kp < -0.5 and r2_kp > 0.7:
        print(f"  => LEAKAGE VERSCHWINDET (MS2 numerisch bestaetigt)")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
