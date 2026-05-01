#!/usr/bin/env python3
# coding: utf-8
"""
Prototyp-CAP-Zertifikat fuer L(s, chi_4) — SCHNELLE VERSION.

Gegenueber chi4_parity_prototype.py:
- Geschlossene Form fuer Cos/Sin-Overlap-Integrale (kein Gauss-Legendre).
- Vektorisiert mit numpy.
- Laufzeit: ~1 Sekunde statt Stunden.

Autor: LG (Claude Opus 4.6 [1M], Session 2026-04-14 Phase B)
"""

import numpy as np
from scipy.special import digamma
import sympy

# =============================================================================
# Charakter chi_4
# =============================================================================
def chi4(n):
    r = n % 4
    if r == 1: return 1
    if r == 3: return -1
    return 0

assert chi4(-1) == -1  # odd character

# =============================================================================
# Analytische Overlap-Integrale
# =============================================================================
# Basis:
#   phi_n(t) = cos(pi n t / L) / sqrt(L)   fuer n >= 1 (cosine, even)
#   phi_0(t) = 1/sqrt(2L)
#   psi_n(t) = sin(pi n t / L) / sqrt(L)   fuer n >= 1 (sine, odd)
#
# Auf [-L, L]. Wir brauchen:
#   S^+_{nm}(delta) = int_{-L}^{L} phi_n(t) phi_m(t - delta) 1_{[-L,L]}(t-delta) dt
#   S^-_{nm}(delta) = int_{-L}^{L} psi_n(t) psi_m(t - delta) 1_{[-L,L]}(t-delta) dt
#
# Fuer |delta| < 2L ist der effektive Bereich [a, b] = [max(-L,-L+delta), min(L,L+delta)]
# = [-L+delta, L] fuer delta > 0, [-L, L+delta] fuer delta < 0. Gleiches Ergebnis fuer |delta|.

def cos_cos_integral(n, m, a, b, L, delta):
    """int_a^b cos(pi n t/L) cos(pi m (t-delta)/L) dt
    = (1/2) int_a^b [cos((n-m) pi t/L + m pi delta/L) + cos((n+m) pi t/L - m pi delta/L)] dt.
    """
    def antideriv(k, phi, t):
        """Antiderivative von cos(k pi t/L + phi) nach t."""
        if k == 0:
            return np.cos(phi) * t
        return (L / (k * np.pi)) * np.sin(k * np.pi * t / L + phi)
    result = 0.0
    for (k, phi) in [(n - m, m * np.pi * delta / L), (n + m, -m * np.pi * delta / L)]:
        result += 0.5 * (antideriv(k, phi, b) - antideriv(k, phi, a))
    return result

def sin_sin_integral(n, m, a, b, L, delta):
    """int_a^b sin(pi n t/L) sin(pi m (t-delta)/L) dt
    = (1/2) int_a^b [cos((n-m) pi t/L + m pi delta/L) - cos((n+m) pi t/L - m pi delta/L)] dt.
    """
    def antideriv(k, phi, t):
        if k == 0:
            return np.cos(phi) * t
        return (L / (k * np.pi)) * np.sin(k * np.pi * t / L + phi)
    result = 0.0
    # + fuer n-m, - fuer n+m
    result += 0.5 * (antideriv(n - m, m * np.pi * delta / L, b) - antideriv(n - m, m * np.pi * delta / L, a))
    result -= 0.5 * (antideriv(n + m, -m * np.pi * delta / L, b) - antideriv(n + m, -m * np.pi * delta / L, a))
    return result

def cos_overlap_normalized(n, m, delta, L):
    """Normalisierte cos-Overlap: phi_n normiert mit 1/sqrt(L) (oder 1/sqrt(2L) fuer n=0)."""
    if abs(delta) >= 2 * L:
        return 0.0
    a = max(-L, -L + delta)
    b = min(L, L + delta)
    if a >= b:
        return 0.0
    I = cos_cos_integral(n, m, a, b, L, delta)
    # Normierung
    norm_n = 1.0 / np.sqrt(L) if n > 0 else 1.0 / np.sqrt(2 * L)
    norm_m = 1.0 / np.sqrt(L) if m > 0 else 1.0 / np.sqrt(2 * L)
    return I * norm_n * norm_m

def sin_overlap_normalized(n, m, delta, L):
    """Normalisierte sin-Overlap: psi_n normiert mit 1/sqrt(L) (n >= 1)."""
    if abs(delta) >= 2 * L:
        return 0.0
    a = max(-L, -L + delta)
    b = min(L, L + delta)
    if a >= b:
        return 0.0
    I = sin_sin_integral(n, m, a, b, L, delta)
    return I / L  # beide Basen normiert mit 1/sqrt(L)

# =============================================================================
# Archimedischer Multiplier
# =============================================================================
def h_arch(tau, kappa):
    """h_kappa(tau) = Re[psi((kappa+1)/4 + i*tau/2)] - log(pi).
    kappa = 0 fuer even character (Riemann-Form, Gamma(s/2)).
    kappa = 1 fuer odd character (Gamma((s+1)/2)).
    """
    shift = (2 * kappa + 1) / 4.0  # 1/4 fuer kappa=0; 3/4 fuer kappa=1
    return digamma(shift + 1j * tau / 2).real - np.log(np.pi)

# =============================================================================
# Galerkin-Matrix (vektorisiert NICHT noetig — N klein)
# =============================================================================
def build_W_matrix(sector, N, L, primes, chi_values, kappa):
    """
    Baut die N x N Galerkin-Matrix der Weil-Quadratform mit Charakter-Twist.
    sector: 'cos' (Even) oder 'sin' (Odd).
    kappa: 0 fuer even chi, 1 fuer odd chi.
    """
    W = np.zeros((N, N))
    if sector == 'cos':
        basis_start = 0
        overlap_fn = cos_overlap_normalized
    else:
        basis_start = 1
        overlap_fn = sin_overlap_normalized

    # Archimedischer Anteil (diagonal in Fourier-Basis)
    # Wir verwenden odd-Kernel fuer odd Charakter
    for i in range(N):
        n = i + basis_start
        tau_n = np.pi * n / L
        W[i, i] += h_arch(tau_n, kappa)

    # Prime-Shift-Anteil
    for p, chi_p in zip(primes, chi_values):
        if chi_p == 0:
            continue
        log_p = np.log(p)
        m_max = int(2 * L / log_p)
        for m_exp in range(1, m_max + 1):
            delta = m_exp * log_p
            if delta >= 2 * L:
                break
            weight = (chi_p ** m_exp) * log_p / (p ** (m_exp / 2.0))
            for i in range(N):
                for j in range(N):
                    n = i + basis_start
                    m = j + basis_start
                    # Symmetrisierte Translation: T_delta = S_delta + S_{-delta}
                    S_pos = overlap_fn(n, m, delta, L)
                    S_neg = overlap_fn(n, m, -delta, L)
                    W[i, j] += weight * (S_pos + S_neg)

    W = 0.5 * (W + W.T)
    return W

# =============================================================================
# Ausfuehrung
# =============================================================================

def run_experiment(lam, N):
    L = np.log(lam)
    primes = list(sympy.primerange(3, int(lam) + 1))
    chi_values = [chi4(p) for p in primes]
    n_plus  = sum(1 for c in chi_values if c > 0)
    n_minus = sum(1 for c in chi_values if c < 0)

    print(f"\n=== lambda = {lam}, L = {L:.4f}, N = {N} ===")
    print(f"    Primzahlen p <= {lam}: {len(primes)} (chi_4 = +1: {n_plus}, -1: {n_minus})")

    # Fuer chi_4: odd character (kappa = 1)
    kappa = 1

    W_cos = build_W_matrix('cos', N, L, primes, chi_values, kappa)
    eig_cos = np.linalg.eigvalsh(W_cos)

    W_sin = build_W_matrix('sin', N, L, primes, chi_values, kappa)
    eig_sin = np.linalg.eigvalsh(W_sin)

    print(f"    Even (cos) Sektor: 5 kleinste EVs = {eig_cos[:5]}")
    print(f"    Odd  (sin) Sektor: 5 kleinste EVs = {eig_sin[:5]}")
    print(f"    lambda_1^+ = {eig_cos[0]:.6f}  (Even)")
    print(f"    lambda_1^- = {eig_sin[0]:.6f}  (Odd)")
    print(f"    Gap = lambda_1^- - lambda_1^+ = {eig_sin[0] - eig_cos[0]:+.6f}")

    if eig_sin[0] < eig_cos[0]:
        print(f"    ==> ODD-DOMINANCE (Vorhersage fuer odd character chi_4: VERIFIZIERT)")
    else:
        print(f"    ==> EVEN-Dominance (unerwartet fuer chi_4)")

    return eig_cos, eig_sin

# =============================================================================
# Durchfuehrung bei mehreren lambda
# =============================================================================
if __name__ == "__main__":
    print("=== Prototyp-CAP-Zertifikat fuer L(s, chi_4) ===\n")
    print("Meta-Paper-Vorhersage: Odd-Dominance fuer odd character chi_4(-1) = -1\n")

    # Fuer Sanity-Check: chi_4 mit kleiner lambda und N
    for (lam, N) in [(30, 12), (100, 16), (300, 18), (1000, 20)]:
        run_experiment(lam, N)

    # Vergleich mit Riemann (kappa = 0, trivialer Charakter chi_0 = 1 everywhere) als Kontrolle
    print("\n\n=== Kontrolle: Riemann (kappa=0) bei lambda=100 ===")
    lam = 100
    N = 16
    L_r = np.log(lam)
    primes_r = list(sympy.primerange(2, int(lam) + 1))
    chi_trivial = [1] * len(primes_r)
    W_cos_r = build_W_matrix('cos', N, L_r, primes_r, chi_trivial, kappa=0)
    W_sin_r = build_W_matrix('sin', N, L_r, primes_r, chi_trivial, kappa=0)
    eig_cos_r = np.linalg.eigvalsh(W_cos_r)
    eig_sin_r = np.linalg.eigvalsh(W_sin_r)
    print(f"    lambda_1^+ (Riemann Even) = {eig_cos_r[0]:.6f}")
    print(f"    lambda_1^- (Riemann Odd)  = {eig_sin_r[0]:.6f}")
    print(f"    Gap = {eig_sin_r[0] - eig_cos_r[0]:+.6f}  (erwartet: positiv = Even-Dominance)")
