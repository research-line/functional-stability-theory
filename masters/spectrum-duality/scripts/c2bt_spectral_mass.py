"""
c2bt_spectral_mass.py

Toleranzfreie Spektraldiagnostik: definiert den Cluster ueber die
von k getragene Spektralmasse, nicht ueber eine starre Eigenwerttol.

Fuer jedes lambda speichert:
  - Volle sortierte Eigenwerte u_j
  - Abstaende Delta_j = u_j - u_0
  - Overlaps m_j = |<q_j, k>|^2
  - Kumulative Masse M(J) = sum_{j<=J} m_j
  - Off-Cluster-Tail 1 - M(J)
  - Effektiver Gap g_eff(eps) = Delta_{J+1} fuer kleinstes J mit M(J) >= 1-eps
  - R_0 / g_eff(eps) fuer verschiedene eps

Skalierung: N = max(40, ceil(N_FACTOR * L)), L = 2*log(lambda)
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

LAMBDAS = [3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0,
           20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0]

EPSILONS = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12]

OUTPUT_JSON = "c2bt_spectral_mass_results.json"


def to_np(M, dim):
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = float(M[i, j])
    return A


def compute_spectral_mass(lam, verbose=True):
    mp.mp.dps = DPS
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

    N = max(40, int(np.ceil(N_FACTOR * L)))
    dim = N + 1

    if verbose:
        print(f"\n  lam={lam:.1f}, L={L:.3f}, N={N}, dim={dim}")

    # Build operator A
    t0 = time.time()
    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")
    dt_build = time.time() - t0

    A_np = to_np(Aq, dim)
    H_np = to_np(Ah, dim)

    # Full eigendecomposition
    t0 = time.time()
    ws, vs = np.linalg.eigh(A_np)
    dt_eig = time.time() - t0

    # Poisson kernel k
    t0 = time.time()
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
    dt_k = time.time() - t0

    # Eigenvalues and overlaps
    u0 = ws[0]
    deltas = ws - u0  # Delta_j = u_j - u_0
    overlaps = np.array([float(kn @ vs[:, j])**2 for j in range(dim)])
    cum_mass = np.cumsum(overlaps)

    # Residual R_0 = ||(A-u0)k||
    Ak = A_np @ kn
    res_w0 = Ak - u0 * kn
    R_0 = float(np.linalg.norm(res_w0))

    # Rayleigh
    rho_k = float(kn @ Ak)
    res_rho = Ak - rho_k * kn
    R_rho = float(np.linalg.norm(res_rho))

    # mu/alpha
    W_np = A_np - H_np
    col0 = W_np[:, 0]
    cn = np.linalg.norm(col0)
    ut = col0 / cn
    alpha = float(ut @ kn)
    mu_alpha = float(ut @ res_w0) / alpha if abs(alpha) > 1e-30 else 0.0

    # Effective gap for various epsilon
    g_eff = {}
    J_eff = {}
    for eps in EPSILONS:
        threshold = 1.0 - eps
        idx = np.searchsorted(cum_mass, threshold)
        if idx < dim - 1:
            J = int(idx)
            g = float(deltas[J + 1])
            g_eff[f"{eps:.0e}"] = g
            J_eff[f"{eps:.0e}"] = J
        else:
            g_eff[f"{eps:.0e}"] = 0.0
            J_eff[f"{eps:.0e}"] = dim - 1

    # R_0/g_eff ratios
    R0_over_geff = {}
    for eps_key, g in g_eff.items():
        if g > 1e-30:
            R0_over_geff[eps_key] = R_0 / g
        else:
            R0_over_geff[eps_key] = float('inf')

    # ||k_perp(eps)|| = sqrt(1 - M(J))
    k_perp_eps = {}
    for eps_key, J in J_eff.items():
        tail = max(0.0, 1.0 - cum_mass[J])
        k_perp_eps[eps_key] = float(np.sqrt(tail))

    # Spectral mass profile: first 20 overlaps + tail summary
    n_report = min(20, dim)
    mass_profile = {
        'overlaps_top20': overlaps[:n_report].tolist(),
        'deltas_top20': deltas[:n_report].tolist(),
        'cum_mass_top20': cum_mass[:n_report].tolist(),
    }

    # Where does 99%, 99.9%, 99.99% of mass sit?
    mass_thresholds = {}
    for pct in [0.99, 0.999, 0.9999, 0.99999, 0.999999]:
        idx = int(np.searchsorted(cum_mass, pct))
        if idx < dim:
            mass_thresholds[f"{pct}"] = {
                'J': idx, 'delta_J': float(deltas[idx]),
                'delta_J1': float(deltas[min(idx+1, dim-1)]),
                'cum_mass_J': float(cum_mass[idx])
            }

    if verbose:
        print(f"    Build: {dt_build:.0f}s, Eig: {dt_eig:.1f}s, k: {dt_k:.0f}s")
        print(f"    R_0={R_0:.4e}, R_rho={R_rho:.4e}, mu/alpha={mu_alpha:.4e}")
        print(f"    m_0={overlaps[0]:.10f}, M(0)={cum_mass[0]:.10f}")
        for eps_key in sorted(g_eff.keys()):
            g = g_eff[eps_key]
            J = J_eff[eps_key]
            ratio = R0_over_geff[eps_key]
            kp = k_perp_eps[eps_key]
            print(f"    eps={eps_key}: J={J}, g_eff={g:.4e}, "
                  f"R0/g_eff={ratio:.4e}, ||k_perp||={kp:.4e}")

    result = {
        'lam': lam, 'L': L, 'N': N, 'dim': dim,
        'u0': float(u0), 'rho_k': float(rho_k),
        'R_0': R_0, 'R_rho': R_rho,
        'mu_alpha': mu_alpha, 'alpha': alpha, 'C_tilde': cn,
        'g_eff': g_eff, 'J_eff': J_eff,
        'R0_over_geff': R0_over_geff,
        'k_perp_eps': k_perp_eps,
        'mass_profile': mass_profile,
        'mass_thresholds': mass_thresholds,
        'dt_build': dt_build, 'dt_eig': dt_eig, 'dt_k': dt_k,
    }

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
    print("C2bt SPECTRAL MASS: Toleranzfreie Clusterdiagnostik")
    print(f"DPS={DPS}, N_FACTOR={N_FACTOR}")
    print(f"Lambda-Range: {LAMBDAS[0]} .. {LAMBDAS[-1]} ({len(LAMBDAS)} Punkte)")
    print(f"Epsilon-Werte: {EPSILONS}")
    print()

    results = []
    t_total = time.time()

    for lam in LAMBDAS:
        try:
            d = compute_spectral_mass(lam)
            results.append(d)
        except Exception as e:
            print(f"  FEHLER bei lam={lam}: {e}")
            import traceback
            traceback.print_exc()

    dt_total = time.time() - t_total
    print(f"\n  Gesamtzeit: {dt_total:.0f}s ({dt_total/60:.1f} min)")

    with open(OUTPUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Ergebnisse gespeichert: {OUTPUT_JSON}")

    # Summary
    print()
    print("=" * 130)
    print("ZUSAMMENFASSUNG: EFFEKTIVER GAP g_eff(eps) UND R_0/g_eff")
    print("=" * 130)

    for eps_key in [f"{e:.0e}" for e in EPSILONS]:
        print(f"\n  --- eps = {eps_key} ---")
        print(f"  {'lam':>5} {'N':>4} {'J':>4} {'g_eff':>12} {'R_0':>12} "
              f"{'R0/g_eff':>12} {'||k_perp||':>12}")
        for d in results:
            J = d['J_eff'].get(eps_key, -1)
            g = d['g_eff'].get(eps_key, 0)
            ratio = d['R0_over_geff'].get(eps_key, float('inf'))
            kp = d['k_perp_eps'].get(eps_key, 0)
            print(f"  {d['lam']:5.1f} {d['N']:4d} {J:4d} {g:12.4e} "
                  f"{d['R_0']:12.4e} {ratio:12.4e} {kp:12.4e}")

        # Power-law fit for R0/g_eff
        lams = [d['lam'] for d in results]
        ratios = []
        for d in results:
            r = d['R0_over_geff'].get(eps_key, 0)
            ratios.append(r if r < 1e10 else 0)
        alpha_fit, C_fit, r2 = power_law_fit(lams, ratios)
        trend = "SINKT" if alpha_fit < -0.1 else "STEIGT" if alpha_fit > 0.1 else "FLACH"
        print(f"  => R0/g_eff ~ lam^({alpha_fit:+.3f}), R^2={r2:.3f} => {trend}")

    # Mass concentration summary
    print()
    print("=" * 130)
    print("MASSENKONZENTRATION: Wo sitzt die Masse von k?")
    print("=" * 130)
    print(f"  {'lam':>5} {'m_0':>12} {'M(0)':>12} {'J(99.99%)':>10} "
          f"{'Delta_J':>12} {'mu/alpha':>12}")
    for d in results:
        m0 = d['mass_profile']['overlaps_top20'][0]
        M0 = d['mass_profile']['cum_mass_top20'][0]
        mt = d['mass_thresholds'].get('0.9999', {})
        J99 = mt.get('J', -1)
        dJ = mt.get('delta_J1', 0)
        print(f"  {d['lam']:5.1f} {m0:12.8f} {M0:12.8f} {J99:10d} "
              f"{dJ:12.4e} {d['mu_alpha']:12.4e}")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
