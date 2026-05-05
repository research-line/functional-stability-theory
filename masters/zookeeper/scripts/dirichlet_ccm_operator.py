"""
dirichlet_ccm_operator.py
=========================

Implementation des Dirichlet-CCM-Operators nach CCM 2025 (arXiv:2511.22755).

KONSTRUKTION (CCM Thm 1.1 + Eq 5.5):

  1. Baue Matrix T = QW_lambda^N im Fourier-Basis V_n, n in {-N,...,N}.
  2. Finde Grundzustand xi (kleinster Eigenwert, even Eigenvektor).
  3. Die Eigenwerte von D^(lambda,N)_log = D - |D*xi><eta|  sind die
     NULLSTELLEN der Funktion
        F(s) = sum_{j=-N}^N xi_j / (j - s)
     Das folgt aus Eq (5.5) in CCM 2025.

Fuer Riemann chi=chi_0: CCM 2025 Kap 6 berichtet Fehler 2.5e-55 fuer erste
Nullstelle bei lambda=sqrt(14), Primen p<=13.

Fuer Dirichlet chi: analog, mit
  - arch-Term mit log(q/pi) + parity-korrektem Gamma-Kern
  - chi(p)^m Gewichten in Prim-Summe
  - W_0,2-Term entfaellt (L(s,chi) hat keinen Pol)
  - Grundzustand xi hat parity-bestimmte Symmetrie

TEST-ZIELE:
  (Sanity) chi_0 bei lambda=sqrt(14), N=20: reproduziere 14.134725...
  (Dirichlet) chi_4 bei lambda=sqrt(14), N=20: teste gegen 6.020948...

Autor: LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung)
Datum: 2026-04-18
Ref: Connes-Consani-Moscovici, "Zeta Spectral Triples", arXiv:2511.22755
"""

import numpy as np
from scipy.linalg import eigh
from scipy.special import digamma
from scipy.optimize import brentq
import json
from pathlib import Path

# =================================================================
# 0. Konfiguration
# =================================================================

LAMBDA_DEFAULT = np.sqrt(14.0)
N_DEFAULT = 20
RESULTS_DIR = Path(__file__).parent.parent / "_results"

# Bekannte Nullstellen fuer Vergleich
RIEMANN_ZEROS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062]
CHI4_ZEROS = [6.020948, 10.243766, 12.988096, 16.343297, 18.291996]

# =================================================================
# 1. von Mangoldt-Funktion + Charaktere
# =================================================================

def is_prime_power(k):
    """Gibt (p, m) zurueck wenn k = p^m, sonst (None, 0)."""
    if k < 2: return None, 0
    for p in range(2, int(np.sqrt(k)) + 1):
        if k % p == 0:
            m = 0; r = k
            while r % p == 0:
                r //= p; m += 1
            if r == 1:
                return p, m
            else:
                return None, 0
    return k, 1  # k selbst ist prim

def mangoldt(k):
    """von Mangoldt: log(p) wenn k=p^m, sonst 0."""
    p, m = is_prime_power(k)
    return np.log(p) if p else 0.0

def chi_trivial(n):
    return 1.0

def chi_4(n):
    # Kronecker(-4, n): 1 fuer n=1 mod 4, -1 fuer n=3 mod 4, 0 sonst
    m = n % 4
    if m == 1: return 1.0
    if m == 3: return -1.0
    return 0.0

def chi_5(n):
    # Legendre (n/5): 1,-1,-1,1 fuer n=1,2,3,4 mod 5
    m = n % 5
    table = [0, 1, -1, -1, 1]
    return float(table[m])

def chi_8(n):
    # Kronecker(+8): 1,-1,-1,1 fuer n=1,3,5,7 mod 8
    m = n % 8
    table = [0, 1, 0, -1, 0, -1, 0, 1]
    return float(table[m])

def chi_33(n):
    # Kronecker(+33)
    vals = [0, 1, 1, 0, 1, -1, 0, -1, 1, 0, -1, 0, 0, -1, -1, 0, 1,
            1, 0, -1, -1, 0, 0, -1, 0, 1, -1, 0, -1, 1, 0, 1, 1]
    return float(vals[n % 33])

# =================================================================
# 2. Die q-Funktion (CCM Eq 2.9)
# =================================================================

def q_kernel(m, n, y, L):
    """
    q(U_m, U_n)(y) fuer y in [0, L]:
      m != n: (sin(2*pi*m*y/L) - sin(2*pi*n*y/L)) / (pi*(n-m))
      m == n: 2*(1 - y/L)*cos(2*pi*n*y/L)
    """
    if m == n:
        return 2.0 * (1.0 - y/L) * np.cos(2.0*np.pi*n*y/L)
    else:
        return (np.sin(2.0*np.pi*m*y/L) - np.sin(2.0*np.pi*n*y/L)) / (np.pi*(n-m))

# =================================================================
# 3. Matrix-Beitraege zu QW_lambda^N
# =================================================================

def prim_matrix(lam, N, chi_func, q_mod):
    """
    Prim-Beitrag zur QW-Matrix:
      -sum_{k prime power, (k,q)=1} chi(k) * Lambda(k)/sqrt(k) * q(U_n, U_m)(log k)

    Das Minus kommt aus QW = ... - sum_k Lambda(k) <f|T(k)f>.
    """
    L = 2.0 * np.log(lam)
    k_max = int(lam**2)
    size = 2*N + 1
    M = np.zeros((size, size))
    for k in range(2, k_max + 1):
        lam_k = mangoldt(k)
        if lam_k == 0: continue
        if q_mod > 1:
            p, _ = is_prime_power(k)
            if p is not None and q_mod % p == 0:
                continue
        chi_k = chi_func(k)
        if chi_k == 0: continue
        y = np.log(k)
        if y >= L: continue
        for n_idx in range(size):
            n = n_idx - N
            for m_idx in range(size):
                m = m_idx - N
                M[n_idx, m_idx] -= chi_k * lam_k / np.sqrt(k) * q_kernel(m, n, y, L)
    return M

def arch_matrix_simple(lam, N, parity, q_mod):
    """
    Vereinfachter arch-Term: Rank-1 Conductor-Beitrag log(q/pi)*|V_n(0)|^2
    plus digamma-Integrand grob numerisch.

    Fuer praezise Form siehe CCM 2025 Prop 4.3 (Hurwitz-Lerch + digamma-Formel).
    Diese Vereinfachung erhaelt die dominante Struktur.
    """
    L = 2.0 * np.log(lam)
    size = 2*N + 1

    # digamma-Integrand numerisch ueber t in [-10*L, 10*L]
    arg = 0.25 if parity == +1 else 0.75
    t_grid = np.linspace(-10.0*L, 10.0*L, 4000)
    dt = t_grid[1] - t_grid[0]
    partial_theta = 0.5 * (np.log(q_mod/np.pi) if q_mod > 1 else -np.log(np.pi))
    partial_theta = partial_theta + 0.5 * np.real(digamma(arg + 1j*t_grid/2.0))
    # 2 partial theta / (2 pi) als Gewicht im Mellin-Integral
    weight = 2.0 * partial_theta / (2.0 * np.pi)

    # Fourier-Transformierte von V_n: annaeherungsweise sinc-Form
    # hat_V_n(t) = L^{-1/2} * (-1)^n * 2 sin(pi*n - tL/2) / (2pi*n/L - t)
    # = L^{1/2} * (-1)^n * sinc((n - tL/(2pi)))
    M = np.zeros((size, size))
    for n_idx in range(size):
        n = n_idx - N
        hat_Vn = np.sqrt(L) * (-1)**n * np.sinc(n - t_grid*L/(2*np.pi))
        for m_idx in range(size):
            m = m_idx - N
            if m > n: continue
            hat_Vm = np.sqrt(L) * (-1)**m * np.sinc(m - t_grid*L/(2*np.pi))
            # |hat f|^2-artig, aber hier: hat_V_n * hat_V_m
            integral = np.sum(hat_Vn * hat_Vm * weight) * dt
            M[n_idx, m_idx] = integral
            M[m_idx, n_idx] = integral
    return M

def W02_matrix(lam, N):
    """
    Rank-2 pol-Beitrag (nur fuer Riemann chi_0), CCM Eq 4.2.

    W_{0,2}(V_n, V_m) = 32 L sinh^2(L/4) (L^2 - 16*pi^2*m*n) / ((L^2 + 16*pi^2*m^2)(L^2 + 16*pi^2*n^2))

    Ergibt Beitrag 2 Re(hat_f(i/2) hat_f(-i/2)) in QW.
    """
    L = 2.0 * np.log(lam)
    size = 2*N + 1
    M = np.zeros((size, size))
    pre = 32.0 * L * np.sinh(L/4.0)**2
    for n_idx in range(size):
        n = n_idx - N
        for m_idx in range(size):
            m = m_idx - N
            num = L**2 - 16.0*np.pi**2*m*n
            den = (L**2 + 16.0*np.pi**2*m**2) * (L**2 + 16.0*np.pi**2*n**2)
            M[n_idx, m_idx] = pre * num / den
    return M

def QW_matrix(lam, N, chi_func, q_mod, parity, include_W02=False, include_arch=True):
    """Gesamt-QW-Matrix: arch + prim [+ W_02 falls chi_0]."""
    M = prim_matrix(lam, N, chi_func, q_mod)
    if include_arch:
        M = M + arch_matrix_simple(lam, N, parity, q_mod)
    if include_W02:
        M = M + W02_matrix(lam, N)
    # Symmetrisieren
    M = 0.5 * (M + M.T)
    return M

# =================================================================
# 4. Grundzustand finden
# =================================================================

def find_ground_state(T, target_parity='even'):
    """
    Kleinster Eigenwert mit richtiger Paritaet.

    Paritaet: gamma(V_n) = V_{-n}. Even Eigenvektor: xi_{-n} = xi_n.
    Fuer numerisches: check dass v[-N] * sign == v[+N].
    """
    w, V = eigh(T)
    # Suche unter den niedrigsten Eigenwerten den ersten mit richtiger Paritaet
    N = (len(w) - 1) // 2
    for k in range(len(w)):
        v = V[:, k]
        # Paritaets-Check: vergleiche v[n_idx] mit v[-n_idx=2N-n_idx]
        mid = N
        is_even = np.allclose(v, v[::-1], atol=1e-6)
        is_odd = np.allclose(v, -v[::-1], atol=1e-6)
        if target_parity == 'even' and is_even:
            return w[k], v, k
        if target_parity == 'odd' and is_odd:
            return w[k], v, k
    # Fallback: nimm den kleinsten
    return w[0], V[:, 0], 0

# =================================================================
# 5. F(s) = sum xi_j / (j-s) und Nullstellen
# =================================================================

def F_function(z, xi, N, L):
    """
    F(z) = sum_{j=-N}^N xi[j+N] / (z - 2*pi*j/L)   (CCM 2025 Eq 5.25).

    Nullstellen von F(z) sind die Eigenwerte des pertubierten
    Skalierungs-Operators D^(lambda,N)_log, die die L-Nullstellen approximieren.
    """
    total = 0.0
    for j in range(-N, N+1):
        pole = 2.0 * np.pi * j / L
        if abs(z - pole) < 1e-14:
            return np.inf
        total += xi[j+N] / (z - pole)
    return total

def find_zeros_of_F(xi, N, L, search_range=(-50, 50), n_samples=40000):
    """
    Finde Nullstellen von F(z) = sum_j xi_j / (z - 2*pi*j/L).

    Pole bei z = 2*pi*j/L fuer j in {-N,...,N}. Zwischen den Polen
    sind die Spektral-Nullstellen, die die L-Nullstellen approximieren.
    """
    poles = [2.0 * np.pi * j / L for j in range(-N, N+1)]
    z_test = np.linspace(search_range[0], search_range[1], n_samples)
    # Entferne Punkte zu nah an den Polen
    mask = np.ones_like(z_test, dtype=bool)
    for pole in poles:
        mask &= np.abs(z_test - pole) > 1e-3
    z_test = z_test[mask]

    F_vals = np.array([F_function(z, xi, N, L) for z in z_test])

    zeros = []
    for i in range(len(z_test)-1):
        if F_vals[i] * F_vals[i+1] < 0:
            # Check ob Pol dazwischen
            pole_between = any(z_test[i] < p < z_test[i+1] for p in poles)
            if not pole_between:
                try:
                    root = brentq(lambda z: F_function(z, xi, N, L),
                                  z_test[i], z_test[i+1], xtol=1e-10)
                    zeros.append(root)
                except:
                    pass
    return sorted(zeros)

# =================================================================
# 6. Hauptprogramm
# =================================================================

def test_character(name, chi_func, q_mod, parity, known_zeros, lam, N,
                   include_W02=False, include_arch=True):
    """Komplettes Test: baue T, finde xi, finde Nullstellen, vergleiche."""
    print(f"\n{'='*70}")
    print(f"TEST: {name}  (q={q_mod}, parity={parity:+d}, lambda={lam:.4f}, N={N})")
    print(f"{'='*70}")

    T = QW_matrix(lam, N, chi_func, q_mod, parity,
                  include_W02=include_W02, include_arch=include_arch)
    print(f"QW-Matrix-Shape: {T.shape}, Norm: {np.linalg.norm(T):.4f}")

    # Grundzustand
    target_parity = 'even' if parity == +1 else 'odd'
    eps_N, xi, k_gs = find_ground_state(T, target_parity=target_parity)
    print(f"Grundzustand: kleinster Eigenwert mit {target_parity} Paritaet bei Index {k_gs}")
    print(f"  epsilon_N = {eps_N:.6f}")
    print(f"  xi-Paritaets-Check: xi == reverse(xi)? {np.allclose(xi, xi[::-1], atol=1e-6)}")

    # Normalisierung: <eta|xi> = 1, mit eta = sum V_j (d.h. eta hat Eintraege alle 1 im V_n-Basis)
    eta = np.ones_like(xi)
    xi_sum = np.sum(xi)
    if abs(xi_sum) > 1e-10:
        xi_normed = xi / xi_sum
    else:
        xi_normed = xi
    print(f"  Normalisierung: sum(xi) = {xi_sum:.6f}")

    # Nullstellen von F
    L = 2.0 * np.log(lam)
    z_max = 2.0 * np.pi * N / L  # hoechster Pol
    zeros = find_zeros_of_F(xi_normed, N, L, search_range=(0.5, z_max), n_samples=40000)
    print(f"\nNullstellen von F(z) in (0.5, {z_max:.2f}) (L={L:.4f}, 2pi/L={2*np.pi/L:.4f}):")
    print(f"  Gefunden: {len(zeros)}")
    for i, z in enumerate(zeros[:7]):
        known = known_zeros[i] if i < len(known_zeros) else None
        if known is not None:
            err = z - known
            print(f"  s_{i+1} = {z:.6f}   (bekannt: {known:.6f}, Fehler: {err:+.2e})")
        else:
            print(f"  s_{i+1} = {z:.6f}")

    return {
        'name': name, 'lambda': float(lam), 'N': N,
        'epsilon_N': float(eps_N),
        'zeros_found': [float(z) for z in zeros[:10]],
        'zeros_known': known_zeros[:5],
    }

def main():
    print("="*70)
    print("dirichlet_ccm_operator.py")
    print("CCM 2025 (arXiv:2511.22755) Rang-1-Konstruktion")
    print("="*70)

    lam = LAMBDA_DEFAULT
    N = N_DEFAULT
    results = []

    # TEST 1: Riemann (Sanity-Check)
    res = test_character('chi_0 (Riemann)', chi_trivial, 1, +1,
                         RIEMANN_ZEROS, lam, N,
                         include_W02=True, include_arch=True)
    results.append(res)

    # TEST 2: chi_4 (odd)
    res = test_character('chi_4', chi_4, 4, -1,
                         CHI4_ZEROS, lam, N,
                         include_W02=False, include_arch=True)
    results.append(res)

    # TEST 3: chi_5 (even)
    from pathlib import Path
    atlas_path = Path(__file__).parent.parent.parent / "dirichlet_atlas" / "_results" / "zeros_all_chars.json"
    atlas_zeros = json.load(open(atlas_path))
    res = test_character('chi_5', chi_5, 5, +1,
                         atlas_zeros['chi_5'][:5], lam, N,
                         include_W02=False, include_arch=True)
    results.append(res)

    # TEST 4: chi_33 (even, kritischer Atlas-Charakter)
    res = test_character('chi_33', chi_33, 33, +1,
                         atlas_zeros['chi_33'][:5], lam, N,
                         include_W02=False, include_arch=True)
    results.append(res)

    # Speichern
    out = RESULTS_DIR / "DIRICHLET_CCM_OPERATOR_2026-04-18.json"
    RESULTS_DIR.mkdir(exist_ok=True)
    with open(out, 'w') as f:
        json.dump({'config': {'lambda': float(lam), 'N': N},
                   'results': results}, f, indent=2)
    print(f"\nRohdaten: {out}")

if __name__ == "__main__":
    main()
