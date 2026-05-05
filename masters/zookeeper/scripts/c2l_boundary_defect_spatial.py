"""
c2l_boundary_defect_spatial.py — Räumliche Lokalisierung des Poisson-Residuals

Kernfrage: Ist der P₀-Residual (H - w_min)k am Rand x≈0, x≈L konzentriert?
Falls ja, ist dies ein reiner Trunkierungs-Randdefekt (Boundary-Lemma).

Methode:
1. Berechne r_P0 = P₀(H - w_min)k in Fourier-Basis
2. Fourier-Spektrum: niedrige vs hohe Frequenzen
3. Rekonstruiere r_P0(x) auf x-Gitter [0, L]
4. Randanteil: Energie in [0, δ] ∪ [L-δ, L] vs Interior
5. Profil |r_P0(x)|/|k(x)| (relativ zum Poisson-Kernel)
6. λ-Skalierung

Autor: LG (Opus 4.6, Session 27)
Datum: 2026-04-19
"""

import os, sys, time
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, diagonalize_mp, chi_trivial,
)
from c2_approximation_test import (
    k_lambda_value, inner_product, norm,
)
from c2_poisson_decomposition import project_to_fourier

DPS = int(os.environ.get("DPS", 50))
CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
    {"lam": 7.0, "N": 85},
]
NGRID = 400


def run(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

    print(f"\n{'='*72}")
    print(f"  BOUNDARY-DEFECT SPATIAL:  lambda={lam}, N={N}, L={L:.4f}")
    print(f"{'='*72}")

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
    C_tilde = float(cn * cn)

    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
    alpha = float(inner_product(ut, kn, dim))

    ws, _ = diagonalize_mp(Aq, verbose=False)
    w_min = float(ws[0])
    print(f"  Build+Diag: {time.time()-t0:.0f}s")
    print(f"  L = {L:.4f}, alpha = {alpha:.6f}, w_min = {w_min:.6f}, C_tilde = {C_tilde:.2f}")

    # Primes included
    import math
    lam_sq = lam**2
    primes_in = [p for p in range(2, int(lam_sq) + 1)
                 if all(p % d != 0 for d in range(2, int(p**0.5) + 1))]
    max_logp = math.log(max(primes_in)) if primes_in else 0
    print(f"  Primzahlen: {len(primes_in)} (bis p={max(primes_in)}), max log(p) = {max_logp:.3f}")
    print(f"  Randdichte: max(log p)/L = {max_logp/L:.3f}")

    # === Residual ===
    Hk = mp.matrix(dim, 1)
    for i in range(dim):
        Hk[i, 0] = sum(H[i, j] * kn[j, 0] for j in range(dim))

    res_full = mp.matrix(dim, 1)
    for i in range(dim):
        res_full[i, 0] = Hk[i, 0] - mp.mpf(w_min) * kn[i, 0]

    res_u = float(inner_product(ut, res_full, dim))
    res_P0 = mp.matrix(dim, 1)
    for i in range(dim):
        res_P0[i, 0] = res_full[i, 0] - mp.mpf(res_u) * ut[i, 0]

    res_P0_norm = float(norm(res_P0, dim))
    print(f"\n  ||res_full|| = {float(norm(res_full, dim)):.6e}")
    print(f"  |<ũ, res>|   = {abs(res_u):.6e}")
    print(f"  ||res_P0||   = {res_P0_norm:.6e}")
    print(f"  Expected -C̃α = {-C_tilde * alpha:.6e}")
    print(f"  Ratio <ũ,r>/(-C̃α) = {res_u / (-C_tilde * alpha):.6f}")

    # === Fourier-Spektrum ===
    r_coeffs = np.array([float(res_P0[n, 0]) for n in range(dim)])
    k_coeffs = np.array([float(kn[n, 0]) for n in range(dim)])
    u_coeffs = np.array([float(ut[n, 0]) for n in range(dim)])

    r2 = r_coeffs**2
    total_r2 = np.sum(r2)

    print(f"\n  === FOURIER-SPEKTRUM DES P₀-RESIDUALS ===")
    print(f"  {'n-Bereich':>20}  {'Σ|r_n|²':>12}  {'Anteil':>8}")
    print(f"  {'-'*44}")

    nq = dim // 4
    ranges = [(0, nq, "low"), (nq, 2*nq, "mid-low"), (2*nq, 3*nq, "mid-high"), (3*nq, dim, "high")]
    for a, b, label in ranges:
        s = np.sum(r2[a:b])
        print(f"  n={a:3d}..{b-1:3d} ({label:>8})  {s:12.4e}  {100*s/total_r2:7.2f}%")

    top5 = np.argsort(r2)[-5:][::-1]
    print(f"\n  Top-5 Fourier-Moden:")
    for idx in top5:
        print(f"    n={idx:3d}: |r_n|² = {r2[idx]:.4e} ({100*r2[idx]/total_r2:.1f}%)")

    # === Spatial Reconstruction ===
    print(f"\n  === RÄUMLICHE REKONSTRUKTION ===")
    x_grid = np.linspace(0.01 * L, 0.99 * L, NGRID)
    # Even-Basis: phi_0 = 1/sqrt(L), phi_n = sqrt(2/L) cos(2*pi*n*x/L)
    inv_sqrt_L = 1.0 / np.sqrt(L)
    sqrt_2_over_L = np.sqrt(2.0 / L)
    ns = np.arange(dim)
    basis_mat = np.zeros((NGRID, dim))
    basis_mat[:, 0] = inv_sqrt_L
    for n in range(1, dim):
        basis_mat[:, n] = sqrt_2_over_L * np.cos(2.0 * np.pi * n * x_grid / L)

    r_x = basis_mat @ r_coeffs
    k_x = basis_mat @ k_coeffs

    dx = x_grid[1] - x_grid[0]
    r2_x = r_x**2
    total_E = np.sum(r2_x) * dx

    # Boundary fractions
    print(f"\n  Randanteil der P₀-Residual-Energie:")
    print(f"  {'δ/L':>6}  {'δ':>6}  {'E_left':>10}  {'E_right':>10}  {'E_int':>10}  {'Rand%':>8}")
    print(f"  {'-'*58}")

    boundary_fracs = {}
    for frac in [0.10, 0.15, 0.20, 0.25, 0.30, 0.40]:
        delta = frac * L
        left_mask = x_grid < delta
        right_mask = x_grid > L - delta
        interior_mask = ~left_mask & ~right_mask

        if not np.any(interior_mask):
            continue

        E_left = np.sum(r2_x[left_mask]) * dx
        E_right = np.sum(r2_x[right_mask]) * dx
        E_int = np.sum(r2_x[interior_mask]) * dx
        bf = (E_left + E_right) / total_E
        boundary_fracs[frac] = bf
        print(f"  {frac:6.2f}  {delta:6.3f}  {E_left:10.4e}  {E_right:10.4e}  {E_int:10.4e}  {100*bf:7.2f}%")

    # Max and profile
    imax = np.argmax(np.abs(r_x))
    x_max = x_grid[imax]
    print(f"\n  Max |r_P0(x)| bei x = {x_max:.4f} (x/L = {x_max/L:.3f})")

    print(f"\n  Profil |r_P0(x)| vs |k(x)|:")
    print(f"  {'x/L':>6}  {'x':>8}  {'|r(x)|':>12}  {'|k(x)|':>12}  {'|r|/|k|':>12}")
    print(f"  {'-'*56}")
    for frac in np.linspace(0.05, 0.95, 10):
        ix = min(int(frac * NGRID), NGRID - 1)
        x = x_grid[ix]
        rval = abs(r_x[ix])
        kval = abs(k_x[ix])
        ratio = rval / kval if kval > 1e-30 else float("inf")
        print(f"  {frac:6.2f}  {x:8.4f}  {rval:12.4e}  {kval:12.4e}  {ratio:12.4e}")

    # Symmetry
    n_half = NGRID // 2
    left_half = r_x[:n_half]
    right_half = r_x[-n_half:][::-1]
    rms_r = np.sqrt(np.mean(r_x**2))
    if rms_r > 1e-30:
        anti_sym = np.sqrt(np.mean((left_half + right_half)**2)) / rms_r
        sym_ratio = np.sqrt(np.mean((left_half - right_half)**2)) / rms_r
    else:
        anti_sym = sym_ratio = 0.0
    print(f"\n  Symmetrie um L/2:")
    print(f"    r(x) ≈ +r(L-x) (symm):  RMS-Abw = {sym_ratio:.4f}")
    print(f"    r(x) ≈ -r(L-x) (anti):  RMS-Abw = {anti_sym:.4f}")
    if anti_sym < sym_ratio:
        print(f"    => r_P0 ist näher an ANTI-symmetrisch")
    else:
        print(f"    => r_P0 ist näher an SYMMETRISCH")

    # Peak structure
    peaks = []
    for i in range(1, NGRID - 1):
        if abs(r_x[i]) > abs(r_x[i-1]) and abs(r_x[i]) > abs(r_x[i+1]):
            peaks.append((x_grid[i], r_x[i]))
    peaks.sort(key=lambda p: abs(p[1]), reverse=True)
    print(f"\n  Top-5 Peaks von r_P0(x):")
    for x_pk, r_pk in peaks[:5]:
        print(f"    x = {x_pk:.4f} (x/L = {x_pk/L:.3f}): r = {r_pk:+.4e}")

    return {
        "lam": lam,
        "L": L,
        "res_P0_norm": res_P0_norm,
        "boundary_fracs": boundary_fracs,
        "x_max_over_L": x_max / L,
        "max_logp_over_L": max_logp / L,
    }


if __name__ == "__main__":
    results = []
    for cfg in CONFIGS:
        r = run(cfg["lam"], cfg["N"])
        results.append(r)

    if len(results) >= 2:
        print(f"\n{'='*72}")
        print(f"  λ-SKALIERUNG")
        print(f"{'='*72}")

        import math

        print(f"\n  {'Metrik':>25}  ", end="")
        for r in results:
            print(f"{'lam='+str(r['lam']):>12}  ", end="")
        print(f"{'Exponent':>10}")

        metrics = [
            ("||res_P0||", "res_P0_norm"),
            ("max(log p)/L", "max_logp_over_L"),
            ("x_max/L", "x_max_over_L"),
        ]
        for name, key in metrics:
            vals = [r[key] for r in results]
            print(f"  {name:>25}  ", end="")
            for v in vals:
                print(f"{v:12.4e}  ", end="")
            if vals[0] > 0 and vals[-1] > 0:
                lams = [r["lam"] for r in results]
                exp = math.log(vals[-1]/vals[0]) / math.log(lams[-1]/lams[0])
                print(f"{exp:10.2f}")
            else:
                print()

        print(f"\n  Rand-Fraktion bei verschiedenen δ/L:")
        for frac in [0.15, 0.20, 0.25, 0.30]:
            print(f"    δ/L = {frac}:  ", end="")
            for r in results:
                bf = r["boundary_fracs"].get(frac, float("nan"))
                print(f"λ={r['lam']}: {100*bf:.1f}%   ", end="")
            print()

    print(f"\n=== FERTIG ===")
