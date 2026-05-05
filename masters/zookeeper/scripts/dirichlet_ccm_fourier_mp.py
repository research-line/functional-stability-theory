"""
dirichlet_ccm_fourier_mp.py
===========================

CCM 2025 Fourier-Basis (V_n, |n| <= N), aber ALLE Matrix-Eintraege durchgehend
in mpmath (mp.matrix), nicht in Float64. Damit die angegebene Praezision
10^-60 (bei dps=200) erreicht werden kann.

Bugfix gegenueber dirichlet_ccm_mpmath.py (Session 17):
  - build_WR_matrix verwendete `np.zeros + float(val)` -> Matrix in Float64
  - build_W02_matrix verwendete `np.sinh, np.pi`       -> Matrix in Float64
  - build_Wprime_matrix verwendete `np.log, np.cos`    -> Matrix in Float64
Trotz dps=50 war die effektive Praezision nur 10^-15 (Float64).

Neu: mp.matrix durchgehend, Diagonalisierung mit mp.eigsy (hermitesch).

Autor: LG (Opus 4.7 [1M], Session 18 beta)
Datum: 2026-04-18
Ref:   Connes-Consani-Moscovici, "Zeta Spectral Triples", arXiv:2511.22755, §6
"""

import os
import sys
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np

# =================================================================
# Konfiguration
# =================================================================

LAMBDA_DEFAULT = float(np.sqrt(14.0))
N_DEFAULT = 30
DPS_DEFAULT = 100

LAMBDA = float(os.environ.get("LAMBDA", LAMBDA_DEFAULT))
N = int(os.environ.get("N", N_DEFAULT))
DPS = int(os.environ.get("DPS", DPS_DEFAULT))

RIEMANN_ZEROS = [
    14.134725141734693,
    21.022039638771554,
    25.010857580145688,
    30.424876125859513,
    32.935061587739189,
]

# =================================================================
# 1. Grundfunktionen (mpmath-nativ)
# =================================================================


def rho_even_mp(x):
    """EVEN-Parity Kern: rho(x) = e^{x/2} / (e^x - e^{-x}) = e^{x/2}/(2 sinh x).
    Fuer Riemann-zeta (trivial character) und EVEN Dirichlet-Charaktere.
    """
    return mp.exp(x / 2) / (mp.exp(x) - mp.exp(-x))


def rho_odd_mp(x):
    """ODD-Parity Kern: rho(x) = e^{-x/2} / (e^x - e^{-x}) = e^{-x/2}/(2 sinh x).
    Abgeleitet aus Weil-Formel mit psi(3/4+it/2): Kernel e^{(1-2a)x}/(2 sinh x),
    a=3/4 => e^{-x/2}. (Even: a=1/4 => e^{+x/2}.)
    """
    return mp.exp(-x / 2) / (mp.exp(x) - mp.exp(-x))


def rho_mp(x, parity="even"):
    return rho_even_mp(x) if parity == "even" else rho_odd_mp(x)


def alpha_L_mp(n, L, parity="even"):
    if n == 0:
        return mp.mpf(0)

    def integrand(x):
        return mp.sin(2 * mp.pi * n * x / L) * rho_mp(x, parity)

    return mp.quad(integrand, [0, L]) / mp.pi


def beta_L_mp(n, L, parity="even"):
    def integrand(x):
        return x * mp.cos(2 * mp.pi * n * x / L) * rho_mp(x, parity)

    return mp.quad(integrand, [0, L]) / L


def gamma_L_mp(n, L, c_plus_w, parity="even"):
    def integrand(x):
        return (mp.cos(2 * mp.pi * n * x / L) - mp.exp(-x / 2)) * rho_mp(x, parity)

    return mp.quad(integrand, [0, L]) + c_plus_w


def compute_c_plus_w_mp(L, parity="even"):
    """Normierungs-Konstante c(L) + w(L).
    Fuer EVEN: aus CCM Seite 14/15 (Vereinfachung).
    Fuer ODD: eigene Berechnung (hier vereinfacht, nur Korrektur-Differenz).
    """
    eL2 = mp.exp(L / 2)
    gamma_euler = mp.euler
    if parity == "even":
        term1 = mp.mpf("0.5") * mp.log((eL2 - 1) / (eL2 + 1))
        term2 = mp.atan(eL2)
        term3 = -mp.pi / 4
        term4 = gamma_euler / 2
        term5 = mp.mpf("0.5") * mp.log(8 * mp.pi)
        return term1 + term2 + term3 + term4 + term5
    else:
        # ODD: psi(3/4) = psi(1/4) + pi, plus Abschneidekorrektur fuer rho_odd.
        # c_plus_w_odd = c_plus_w_even + 2*atan(eL2) - 3*pi/2 + log(1+e^{-L})
        term1 = mp.mpf("0.5") * mp.log((eL2 - 1) / (eL2 + 1))
        term2 = mp.atan(eL2)
        term3 = -mp.pi / 4
        term4 = gamma_euler / 2
        term5 = mp.mpf("0.5") * mp.log(8 * mp.pi)
        cpw_even = term1 + term2 + term3 + term4 + term5
        return cpw_even + 2 * mp.atan(eL2) - 3 * mp.pi / 2 + mp.log(1 + mp.exp(-L))


# =================================================================
# 2. Matrix-Bausteine (mp.matrix nativ)
# =================================================================


def build_WR_matrix_mp(N, L_mp, parity="even", verbose=True):
    """Archimedische Weil-Matrix (Prop 4.3) in mp.matrix mit parity-Unterscheidung."""
    if verbose:
        print(f"  Berechne alpha/beta/gamma (parity={parity}) fuer n in [-{N}, {N}]...")
    c_plus_w = compute_c_plus_w_mp(L_mp, parity=parity)
    alpha_vals = {}
    beta_vals = {}
    gamma_vals = {}
    for n in range(-N, N + 1):
        alpha_vals[n] = alpha_L_mp(n, L_mp, parity=parity)
        beta_vals[n] = beta_L_mp(n, L_mp, parity=parity)
        gamma_vals[n] = gamma_L_mp(n, L_mp, c_plus_w, parity=parity)

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
    """CCM Eq (4.2): Rang-2-Pol-Matrix (nur Riemann), mp.matrix."""
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


def q_kernel_mp(m, n, y, L_mp):
    """q(U_m, U_n)(y) fuer y in [0, L], mpmath-nativ."""
    if m == n:
        return 2 * (1 - y / L_mp) * mp.cos(2 * mp.pi * n * y / L_mp)
    else:
        return (
            mp.sin(2 * mp.pi * m * y / L_mp) - mp.sin(2 * mp.pi * n * y / L_mp)
        ) / (mp.pi * (n - m))


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


def build_Wprime_matrix_mp(N, L_mp, chi_func, q_mod, verbose=True):
    """Prim-Beitrag per Eq (4.3), mp.matrix."""
    size = 2 * N + 1
    Wp = mp.matrix(size, size)
    k_max = int(mp.exp(L_mp))  # = lambda^2
    if verbose:
        print(f"  Prim-Summe bis k_max = {k_max}")
    for k in range(2, k_max + 1):
        lam_k = mangoldt_mp(k)
        if lam_k == 0:
            continue
        if q_mod > 1:
            p, _ = is_prime_power(k)
            if p and q_mod % p == 0:
                continue
        chi_k = chi_func(k)
        if chi_k == 0:
            continue
        y = mp.log(mp.mpf(k))
        sqrt_k = mp.sqrt(mp.mpf(k))
        pre = mp.mpf(chi_k) * lam_k / sqrt_k
        for i in range(size):
            n = i - N
            for j in range(size):
                m = j - N
                Wp[i, j] += pre * q_kernel_mp(m, n, y, L_mp)
    return Wp


# =================================================================
# 3. Gesamtmatrix
# =================================================================


def chi_trivial(n):
    return 1


def chi_4(n):
    m = n % 4
    return 1 if m == 1 else (-1 if m == 3 else 0)


def build_QW_mp(N, lam, chi_func=chi_trivial, q_mod=1, include_W02=True,
                sign_WR=-1, parity="even", conductor_correction=False,
                verbose=True):
    """
    QW_lambda^{chi,N} = sign_WR * W_R - sum W_p [+ W_02 for chi_0], mp.matrix.

    CCM Konvention: Psi = W_{0,2} - W_R - sum W_p   => sign_WR = -1.

    Fuer Dirichlet chi (q>1):
      - parity = 'even' (chi(-1)=+1) oder 'odd' (chi(-1)=-1)
      - include_W02 = False (L(s,chi) hat keinen Pol)
      - conductor_correction: wenn True, log(q/pi)/L uniform auf WR addieren
        (Diagnostik zeigt: NICHT verwenden — c_plus_w enthaelt Conductor bereits)
    """
    L_mp = 2 * mp.log(mp.mpf(lam))  # L = log(lambda^2) = 2 log(lambda)

    if verbose:
        print(f"\n  Baue QW bei lambda={lam:.6f}, N={N}, dps={mp.mp.dps}")
        print(f"  L = 2 log(lambda) = {float(L_mp):.6f}")
        print(f"  chi parity = {parity}, q_mod = {q_mod}")

    if verbose:
        print(f"  Baue W_R...")
    WR = build_WR_matrix_mp(N, L_mp, parity=parity, verbose=verbose)

    if conductor_correction and q_mod > 1:
        if verbose:
            print(f"  Conductor-Korrektur: log({q_mod}/pi)/L auf W_R")
        conductor = mp.log(mp.mpf(q_mod) / mp.pi)
        delta = conductor / L_mp
        size_tmp = 2 * N + 1
        for i in range(size_tmp):
            for j in range(size_tmp):
                WR[i, j] += delta

    if verbose:
        print(f"  Baue W_prime...")
    Wp = build_Wprime_matrix_mp(N, L_mp, chi_func, q_mod, verbose=verbose)

    size = 2 * N + 1
    M = mp.matrix(size, size)
    for i in range(size):
        for j in range(size):
            M[i, j] = sign_WR * WR[i, j] - Wp[i, j]

    if include_W02:
        if verbose:
            print(f"  Baue W_02 (Pol)...")
        W02 = build_W02_matrix_mp(N, L_mp)
        for i in range(size):
            for j in range(size):
                M[i, j] += W02[i, j]

    # Symmetrisieren (numerisch robust)
    for i in range(size):
        for j in range(i + 1, size):
            avg = (M[i, j] + M[j, i]) / 2
            M[i, j] = avg
            M[j, i] = avg

    return M, L_mp


# =================================================================
# 4. Diagonalisierung + Grundzustand + F(z)-Nullstellen
# =================================================================


def project_to_parity(M, N, parity="even"):
    """
    Projiziert die (2N+1) x (2N+1) Matrix M auf den even- oder odd-Unterraum
    unter der Symmetrie gamma: V_n <-> V_{-n}.

    Even-Basis:  e_0 = V_0,  e_k = (V_k + V_{-k})/sqrt(2)  fuer k=1..N   (Dim N+1)
    Odd-Basis:   o_k = (V_k - V_{-k})/sqrt(2)              fuer k=1..N   (Dim N)

    Rueckgabe: (M_sub, basis_matrix) mit M_sub = U^T M U, U = Einbettungsmatrix.
    """
    size = 2 * N + 1
    if parity == "even":
        dim = N + 1
        U = mp.matrix(size, dim)
        # e_0 = V_0  (Index N in voller Basis)
        U[N, 0] = mp.mpf(1)
        inv_sqrt2 = 1 / mp.sqrt(mp.mpf(2))
        for k in range(1, N + 1):
            # e_k = (V_k + V_{-k})/sqrt(2): Index N+k und N-k
            U[N + k, k] = inv_sqrt2
            U[N - k, k] = inv_sqrt2
    else:
        dim = N
        U = mp.matrix(size, dim)
        inv_sqrt2 = 1 / mp.sqrt(mp.mpf(2))
        for k in range(1, N + 1):
            # o_k = (V_k - V_{-k})/sqrt(2)
            U[N + k, k - 1] = inv_sqrt2
            U[N - k, k - 1] = -inv_sqrt2

    # M_sub = U^T M U  (bei reeller Matrix; conjugate nicht noetig wenn reell)
    UT = U.T
    # mpmath matrix-multiplikation
    M_sub = UT * M * U
    return M_sub, U


def diagonalize_mp(M, verbose=True):
    """Diagonalisiert eine Hermitesche mp.matrix mit mp.eigsy."""
    if verbose:
        print(f"  Diagonalisiere {M.rows}x{M.cols} mit mp.eigsy (dps={mp.mp.dps})...")
    t0 = time.time()
    E, Q = mp.eigsy(M)  # E: eigenvalues (diagonal vec), Q: eigenvectors (spalten)
    dt = time.time() - t0
    if verbose:
        print(f"    Fertig in {dt:.1f}s")

    # E ist mp.matrix (N x 1), Q ist mp.matrix (N x N)
    N = M.rows
    w = [E[i] for i in range(N)]
    # Sortieren (mp.eigsy gibt unsortiert)
    order = sorted(range(N), key=lambda i: float(w[i]))
    w_sorted = [w[i] for i in order]
    # Q_sorted: Spalten in neuer Reihenfolge
    Q_sorted = mp.matrix(N, N)
    for new_idx, old_idx in enumerate(order):
        for row in range(N):
            Q_sorted[row, new_idx] = Q[row, old_idx]
    return w_sorted, Q_sorted


def find_ground_state_mp(w, Q, delta_N, degen_tol_factor=1e-20):
    """
    Sucht den kleinsten Eigenvektor xi mit Overlap zu delta_N.
    Falls Entartung: projiziere auf delta_N-Richtung im entarteten Cluster.

    delta_N: mp-Vektor (L^{-1/2} * ones) in V-Basis
    """
    N = len(w)
    w_min = w[0]

    # Entartete Cluster identifizieren (sehr eng, um echte Entartung von
    # naher-EWs zu unterscheiden)
    degen_idx = [k for k in range(N) if abs(w[k] - w_min) < degen_tol_factor]

    if len(degen_idx) == 1:
        # Einfacher Grundzustand
        xi = mp.matrix(N, 1)
        for i in range(N):
            xi[i, 0] = Q[i, degen_idx[0]]
        overlap = sum(mp.conj(delta_N[i]) * xi[i, 0] for i in range(N))
        return w_min, xi, len(degen_idx), abs(overlap)

    # Mehrere Eigenvektoren entartet: projiziere delta_N in den entarteten Unterraum
    # xi = sum_k <delta_N, V[:,k]> * V[:,k] (innerhalb entarteten Clusters)
    xi = mp.matrix(N, 1)
    total_overlap = mp.mpf(0)
    for k in degen_idx:
        # Overlap <delta_N, V[:,k]>
        ov = sum(mp.conj(delta_N[i]) * Q[i, k] for i in range(N))
        total_overlap += abs(ov) ** 2
        for i in range(N):
            xi[i, 0] += ov * Q[i, k]

    # Normiere xi
    norm2 = sum(abs(xi[i, 0]) ** 2 for i in range(N))
    if norm2 == 0:
        return w_min, None, len(degen_idx), mp.mpf(0)
    norm = mp.sqrt(norm2)
    for i in range(N):
        xi[i, 0] /= norm

    return w_min, xi, len(degen_idx), mp.sqrt(total_overlap)


def find_F_zeros(xi, N, L_mp, z_search=None, n_zeros=5):
    """
    Findet Nullstellen von F(z) = sum_j xi_j / (z - 2*pi*j/L), per CCM Eq (5.25).

    z ist REELL. F(z) ist reell bei reellen xi_j. Nullstellen sind die
    Riemann-Nullstellen gamma_k (im L(s,chi)-Fall entsprechende Imaginaerteile).

    xi: mp.matrix (2N+1) x 1, in V-Basis (Index -N...+N)
    """
    size = 2 * N + 1

    def F(z):
        z_mp = mp.mpf(z) if not isinstance(z, (mp.mpf, mp.mpc)) else z
        s = mp.mpf(0)
        for idx in range(size):
            j = idx - N
            pole = 2 * mp.pi * j / L_mp
            if abs(z_mp - pole) < mp.mpf("1e-100"):
                return mp.inf
            s += xi[idx, 0] / (z_mp - pole)
        return s

    # Scan entlang reeller z-Achse
    if z_search is None:
        # Poles liegen bei 2*pi*j/L. Zwischen den Polen nach Vorzeichen-Wechsel suchen
        L_f = float(L_mp)
        # Ausreichend feines Gitter
        z_search = np.linspace(0.1, 50.0, 2000)

    vals = []
    for z in z_search:
        Fv = F(z)
        vals.append(float(Fv) if not mp.isinf(Fv) else np.nan)

    # Vorzeichenwechsel finden
    zeros = []
    for k in range(1, len(vals) - 1):
        if np.isnan(vals[k]) or np.isnan(vals[k - 1]):
            continue
        if vals[k - 1] * vals[k] < 0:
            # Refinement via bisection (simpel)
            z_lo, z_hi = z_search[k - 1], z_search[k]
            f_lo, f_hi = vals[k - 1], vals[k]
            for _ in range(60):
                z_mid = (z_lo + z_hi) / 2
                f_mid = float(F(z_mid))
                if f_lo * f_mid < 0:
                    z_hi = z_mid
                    f_hi = f_mid
                else:
                    z_lo = z_mid
                    f_lo = f_mid
            zeros.append((z_lo + z_hi) / 2)
        if len(zeros) >= n_zeros:
            break
    return zeros, vals, None


# =================================================================
# 5. Haupt-Routine
# =================================================================


def main(N=N, lam=LAMBDA, dps=DPS, chi_name="trivial"):
    mp.mp.dps = dps
    print("=" * 70)
    print(f"dirichlet_ccm_fourier_mp.py")
    print(f"  mpmath dps = {mp.mp.dps}")
    print(f"  N          = {N}")
    print(f"  lambda     = {lam:.6f}")
    print(f"  chi        = {chi_name}")
    print("=" * 70)

    if chi_name == "trivial":
        chi_func = chi_trivial
        q_mod = 1
        include_W02 = True
    elif chi_name == "chi_4":
        chi_func = chi_4
        q_mod = 4
        include_W02 = False
    else:
        raise ValueError(f"Unknown chi: {chi_name}")

    t_total = time.time()
    M, L_mp = build_QW_mp(N, lam, chi_func=chi_func, q_mod=q_mod,
                          include_W02=include_W02, sign_WR=-1, verbose=True)

    # delta_N in V-Basis: delta_N = L^{-1/2} * sum V_n
    # Eintraege alle = L^{-1/2}
    size = 2 * N + 1
    delta_N = mp.matrix(size, 1)
    inv_sqrt_L = 1 / mp.sqrt(L_mp)
    for i in range(size):
        delta_N[i, 0] = inv_sqrt_L

    # --- Vollraum-Diagonalisierung (Diagnose, nicht fuer Grundzustand genutzt) ---
    print(f"\n  [Diagnose] Vollraum-Diagonalisierung:")
    w_full, Q_full = diagonalize_mp(M, verbose=True)
    print(f"    Niedrigste 10 Eigenwerte (Vollraum):")
    for k in range(min(10, len(w_full))):
        print(f"      k={k:3d}: w = {float(w_full[k]):+.15e}")

    # --- EVEN-Unterraum (CCM Def 5.3 fordert gamma*xi = xi) ---
    print(f"\n  Projektion auf EVEN-Unterraum (Dim N+1 = {N+1})...")
    M_even, U_even = project_to_parity(M, N, parity="even")
    # Symmetrisieren
    for i in range(N + 1):
        for j in range(i + 1, N + 1):
            avg = (M_even[i, j] + M_even[j, i]) / 2
            M_even[i, j] = avg
            M_even[j, i] = avg

    print(f"  Diagonalisiere EVEN-Matrix ({N+1}x{N+1})...")
    w_even, Q_even = diagonalize_mp(M_even, verbose=True)

    print(f"\n  Niedrigste 15 Eigenwerte (EVEN-Unterraum):")
    for k in range(min(15, len(w_even))):
        if k + 1 < len(w_even):
            gap = float(w_even[k + 1] - w_even[k])
            print(f"    k={k:3d}: w = {float(w_even[k]):+.15e}   gap-to-next = {gap:+.3e}")
        else:
            print(f"    k={k:3d}: w = {float(w_even[k]):+.15e}")

    # Cluster-Projektion im EVEN-Unterraum:
    # Bei endlichem N ist der echte Grundzustand eine LK der fast-entarteten
    # Eigenvektoren. Projiziere delta_N (in Even-Basis) in den Cluster.
    w_min = w_even[0]

    # delta_N im Even-Basis: delta_N_even = U_even^T * delta_N
    # delta_N[i] = 1/sqrt(L) in Vollraum. Berechne Projektion.
    delta_N_even = mp.matrix(N + 1, 1)
    for r in range(N + 1):
        s = mp.mpf(0)
        for row in range(size):
            s += U_even[row, r] * delta_N[row, 0]
        delta_N_even[r, 0] = s

    # Entartungs-Toleranz: Cluster mit Gap > tol zum nachsten
    # Default: nimm alle mit (w[k] - w_min) < 1e-3 (deutlich unter Sprung bei k=10)
    cluster_tol_str = os.environ.get("CLUSTER_TOL", "1e-3")
    cluster_tol = float(cluster_tol_str)
    degen_idx = [k for k in range(N + 1) if abs(float(w_even[k] - w_min)) < cluster_tol]
    print(f"\n  Cluster-Projektion im Even-Unterraum (tol={cluster_tol_str}):")
    print(f"    Cluster-Groesse = {len(degen_idx)}  (k=0..{max(degen_idx)})")

    # xi_even = sum_{k in cluster} <delta_N_even, Q_even[:,k]> * Q_even[:,k]
    xi_even = mp.matrix(N + 1, 1)
    coefs = []
    for k in degen_idx:
        ov = mp.mpf(0)
        for r in range(N + 1):
            ov += delta_N_even[r, 0] * Q_even[r, k]  # reell, conj weglassen
        coefs.append(ov)
        for r in range(N + 1):
            xi_even[r, 0] += ov * Q_even[r, k]

    # Normierung
    norm2 = mp.mpf(0)
    for r in range(N + 1):
        norm2 += xi_even[r, 0] ** 2
    if norm2 > 0:
        norm = mp.sqrt(norm2)
        for r in range(N + 1):
            xi_even[r, 0] /= norm

    # Einbettung xi_full = U_even @ xi_even
    xi_full = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for k_inner in range(N + 1):
            s += U_even[row, k_inner] * xi_even[k_inner, 0]
        xi_full[row, 0] = s

    # Overlap pruefen
    overlap = mp.mpf(0)
    for i in range(size):
        overlap += delta_N[i, 0] * xi_full[i, 0]
    overlap_abs = abs(overlap)

    xi = xi_full
    degen_size = len(degen_idx)

    print(f"    Projektions-Koeffizienten: "
          f"{[f'{float(c):+.3e}' for c in coefs[:6]]}{'...' if len(coefs)>6 else ''}")

    print(f"\n  Grundzustand:")
    print(f"    w_min       = {float(w_min):+.15e}")
    print(f"    |<delta_N, xi>| = {float(overlap_abs):.6e}")

    if xi is not None:
        # Summe der xi_j (xi ist mp.matrix size x 1, Eintraege sind mpc oder mpf)
        sum_xi_mp = mp.mpf(0)
        for i in range(size):
            val = xi[i, 0]
            sum_xi_mp += mp.re(val) if isinstance(val, mp.mpc) else val
        print(f"    sum(xi_j)   = {float(sum_xi_mp):.6e}   "
              f"(sqrt(L) = {float(mp.sqrt(L_mp)):.6f})")

        # Normiere xi so dass sum(xi) = sqrt(L)  (CCM-Konvention, Seite 24)
        if abs(sum_xi_mp) > 1e-40:
            scale = mp.sqrt(L_mp) / sum_xi_mp
            xi_norm = mp.matrix(size, 1)
            for i in range(size):
                xi_norm[i, 0] = xi[i, 0] * scale
            print(f"    xi normiert: sum -> sqrt(L)")
        else:
            xi_norm = xi
            print(f"    WARN: sum(xi) zu klein, xi unnormiert")

        # F(z)-Nullstellen suchen
        print(f"\n  Suche F(z)-Nullstellen bei z in [0, 50]...")
        zeros, vre, vim = find_F_zeros(xi_norm, N, L_mp, n_zeros=8)
        print(f"    Gefundene Nullstellen (Grob-Scan): {[f'{z:.4f}' for z in zeros[:8]]}")

        # Feinsuche um gamma_k via brentq auf Im(F)
        def F_imag_real_arg(z):
            z_mpc = mp.mpc(0, z)
            s = mp.mpf(0)
            # F ist komplex; aber fuer reelle z sollte F imaginaer sein.
            # Wir suchen Re(F) == 0  oder  Im(F) == 0
            for idx in range(size):
                j = idx - N
                pole = 2 * mp.pi * j / L_mp
                s = s + xi_norm[idx, 0] / (mp.mpc(0, z) - pole)
            return float(mp.im(s))

        # Direkte Auswertung F(gamma_k) (sanity — sollte klein sein wenn OK)
        def F_direct(z):
            z_mp = mp.mpf(z)
            s = mp.mpf(0)
            for idx in range(size):
                j = idx - N
                pole = 2 * mp.pi * j / L_mp
                s += xi_norm[idx, 0] / (z_mp - pole)
            return s

        print(f"\n  F(gamma_k) direkt (sollte 0 sein):")
        for k, gamma_k in enumerate(RIEMANN_ZEROS[:5]):
            res = float(F_direct(gamma_k))
            print(f"    k={k+1}: gamma_{k+1} = {gamma_k:.10f}, F = {res:+.3e}")

        print(f"\n  Scan nach Vorzeichenwechseln in F(z) auf [0.1, 50]...")
        zeros2, _, _ = find_F_zeros(xi_norm, N, L_mp, n_zeros=10)
        print(f"  Gefundene Nullstellen:")
        for i, z in enumerate(zeros2[:10]):
            if i < len(RIEMANN_ZEROS):
                gamma_ref = RIEMANN_ZEROS[i]
                err = abs(z - gamma_ref)
                print(f"    z_{i+1} = {z:.12f}   gamma_{i+1} = {gamma_ref:.12f}   "
                      f"err = {err:.3e}")
            else:
                print(f"    z_{i+1} = {z:.12f}")

    dt = time.time() - t_total
    print(f"\n  Gesamtzeit: {dt:.1f}s")

    return w_even, Q_even, xi, M


if __name__ == "__main__":
    main(N=N, lam=LAMBDA, dps=DPS, chi_name=os.environ.get("CHI", "trivial"))
