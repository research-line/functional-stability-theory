#!/usr/bin/env python3
# coding: utf-8
"""
Prototyp-CAP-Zertifikat fuer L(s, chi_4) via v2.1-Twist.

Ziel: Verifikation der Meta-Paper-Vorhersage (Dirichlet-Zeile der C2-Cartography):
"C2^parity closable via v2.1-twist".

Fuer chi_4 (primitiv mod 4, odd character chi_4(-1) = -1) erwarten wir
ODD-DOMINANCE der Weil-Quadratform (Spiegel von Even-Dominance bei Riemann).

Method:
1. Baue Galerkin-Matrix des getwisteten Weil-Operators W_N^{chi_4} bei lambda = 100.
2. Berechne Even-Sektor-Block (cosine-Basis) und Odd-Sektor-Block (sine-Basis).
3. Vergleiche kleinste Eigenwerte.

Dies ist ein FLOATING-POINT-Prototyp (keine Interval-Arithmetic; das kommt erst
in der Paper-Version mit mpmath.iv).

Autor: LG (Claude Opus 4.6 [1M], Session 2026-04-14)
"""

import numpy as np
from scipy.special import digamma
import sympy

# -----------------------------------------------------------------------------
# 1. Charakter chi_4 (odd, primitive, conductor 4)
# -----------------------------------------------------------------------------

def chi4(n):
    """Dirichlet-Charakter mod 4: chi_4(n) = 1 wenn n ≡ 1 (mod 4),
    -1 wenn n ≡ 3 (mod 4), 0 sonst."""
    r = n % 4
    if r == 1:
        return 1
    elif r == 3:
        return -1
    else:
        return 0

# Basic check
assert chi4(1) == 1
assert chi4(3) == -1
assert chi4(5) == 1  # 5 mod 4 = 1
assert chi4(7) == -1  # 7 mod 4 = 3
assert chi4(-1) == chi4(3) == -1  # odd character confirmed

# -----------------------------------------------------------------------------
# 2. Parameter
# -----------------------------------------------------------------------------

LAMBDA = 100.0
L = np.log(LAMBDA)  # L ≈ 4.605
N_MODES = 20  # Galerkin-Dimension pro Sektor

print(f"[Setup] lambda = {LAMBDA}, L = log(lambda) = {L:.4f}, N = {N_MODES}")

# -----------------------------------------------------------------------------
# 3. Primzahlen p <= lambda mit chi_4(p) != 0 (d.h. p ungerade)
# -----------------------------------------------------------------------------

primes = list(sympy.primerange(3, int(LAMBDA) + 1))  # p = 2 ist ausgeschlossen
print(f"[Primes] {len(primes)} ungerade Primzahlen p <= {LAMBDA}: {primes[:10]} ... {primes[-5:]}")

# chi_4(p) pro Primzahl
chi_values = [chi4(p) for p in primes]
print(f"[chi_4 values] {dict(zip(primes[:10], chi_values[:10]))}")

# -----------------------------------------------------------------------------
# 4. Basis-Funktionen und Overlap-Integrale
# -----------------------------------------------------------------------------

def cos_basis(n, t, L):
    """ONB cosine fuer Even-Sektor: phi_n(t) = cos(pi n t / L) / sqrt(L) fuer n > 0;
    phi_0(t) = 1/sqrt(2L)."""
    if n == 0:
        return 1.0 / np.sqrt(2 * L)
    return np.cos(np.pi * n * t / L) / np.sqrt(L)

def sin_basis(n, t, L):
    """ONB sine fuer Odd-Sektor: psi_n(t) = sin(pi n t / L) / sqrt(L) fuer n >= 1."""
    return np.sin(np.pi * n * t / L) / np.sqrt(L)

def shift_overlap_cos(n, m, delta, L):
    """
    S_{nm}^+(delta) = int_{-L+delta}^{L} phi_n(t) phi_m(t - delta) dt.

    Fuer cosine-Basis: geschlossene Form via Produkt-zu-Summe-Identitaet.
    Wir implementieren numerisch via Gauss-Legendre, da wir analytische Formeln
    nicht fuer alle (n, m)-Paare vereinfachen wollen.
    """
    if abs(delta) >= 2 * L:
        return 0.0
    # Integrationsbereich: max(-L, -L+delta) bis min(L, L+delta)
    a = max(-L, -L + delta)
    b = min(L, L + delta)
    if a >= b:
        return 0.0
    # Gauss-Legendre mit 200 Punkten fuer hohe Genauigkeit
    x_gl, w_gl = np.polynomial.legendre.leggauss(200)
    # Skaliere auf [a, b]
    t = 0.5 * (b - a) * x_gl + 0.5 * (a + b)
    w = 0.5 * (b - a) * w_gl
    integrand = cos_basis(n, t, L) * cos_basis(m, t - delta, L)
    return np.sum(w * integrand)

def shift_overlap_sin(n, m, delta, L):
    """S_{nm}^-(delta) fuer sine-Basis."""
    if abs(delta) >= 2 * L:
        return 0.0
    a = max(-L, -L + delta)
    b = min(L, L + delta)
    if a >= b:
        return 0.0
    x_gl, w_gl = np.polynomial.legendre.leggauss(200)
    t = 0.5 * (b - a) * x_gl + 0.5 * (a + b)
    w = 0.5 * (b - a) * w_gl
    integrand = sin_basis(n, t, L) * sin_basis(m, t - delta, L)
    return np.sum(w * integrand)

# -----------------------------------------------------------------------------
# 5. Archimedischer Kern fuer odd character chi_4(-1) = -1
# -----------------------------------------------------------------------------
# Gamma-Faktor: Gamma((s+1)/2) / pi^((s+1)/2)
# Weil-Multiplier h^-(tau) = Re[digamma((1 + i*tau)/2 + 1/2)] - log(pi)
# Fuer odd character: h^-(tau) = Re[digamma(3/4 + i*tau/2)] - log(pi)

def h_arch_chi(tau):
    """Archimedischer Weil-Multiplier fuer odd character (chi_4).
    Verschiebung 3/4 statt 1/4 im Digamma-Argument."""
    return digamma(0.75 + 1j * tau / 2).real - np.log(np.pi)

# Basic sanity check
print(f"[Archimedean multiplier] h^chi(0) = {h_arch_chi(0):.4f}")
print(f"[Archimedean multiplier] h^chi(10) = {h_arch_chi(10):.4f}")

# -----------------------------------------------------------------------------
# 6. Galerkin-Matrix
# -----------------------------------------------------------------------------

def build_galerkin_matrix_chi(sector, N, L, lambda_val, primes, chi_values):
    """
    Baut die N x N Galerkin-Matrix der Weil-Quadratform W^{chi_4} im gegebenen
    Sektor ('cos' = Even oder 'sin' = Odd).

    Struktur: W = W_arch + W_prime
    W_prime_{nm} = sum_{p, m_exp}  chi(p)^{m_exp} * (log p) / p^{m_exp/2}
                 * [S_{nm}(m_exp*log p) + S_{nm}(-m_exp*log p)]
    """
    W = np.zeros((N, N))

    # Wahl des Basis-Startindex und der Overlap-Funktion
    if sector == 'cos':
        basis_start = 0  # cosine-Basis beginnt bei n=0
        overlap = shift_overlap_cos
    elif sector == 'sin':
        basis_start = 1  # sine-Basis beginnt bei n=1
        overlap = shift_overlap_sin
    else:
        raise ValueError(f"Unknown sector: {sector}")

    # Archimedischer Anteil: diagonal in Fourier-Basis mit Multiplikator h(tau_n)
    # tau_n = pi * n / L fuer cosine (n >= 0) oder sine (n >= 1)
    for i in range(N):
        n = i + basis_start
        tau_n = np.pi * n / L
        W[i, i] += h_arch_chi(tau_n)

    # Prime-Shift-Anteil
    for p, chi_p in zip(primes, chi_values):
        if chi_p == 0:
            continue
        log_p = np.log(p)
        m_exp_max = int(2 * L / log_p)  # m_exp * log p <= 2L
        for m_exp in range(1, m_exp_max + 1):
            if m_exp * log_p > 2 * L:
                break
            delta = m_exp * log_p
            weight = (chi_p ** m_exp) * log_p / (p ** (m_exp / 2.0))
            # Symmetrisierte Translation: T_delta = S_delta + S_{-delta}
            for i in range(N):
                for j in range(N):
                    n = i + basis_start
                    m = j + basis_start
                    S_plus = overlap(n, m, delta, L)
                    S_minus = overlap(n, m, -delta, L)
                    W[i, j] += weight * (S_plus + S_minus)

    # Symmetrisierung (numerisch)
    W = 0.5 * (W + W.T)
    return W

# -----------------------------------------------------------------------------
# 7. Ausfuehrung
# -----------------------------------------------------------------------------

print("\n=== Building Galerkin matrix for chi_4, EVEN sector (cosine basis) ===")
W_cos = build_galerkin_matrix_chi('cos', N_MODES, L, LAMBDA, primes, chi_values)
eigvals_cos = np.linalg.eigvalsh(W_cos)
print(f"[Cosine sector] Smallest 5 eigenvalues: {eigvals_cos[:5]}")
print(f"[Cosine sector] lambda_1^+ = {eigvals_cos[0]:.6f}")

print("\n=== Building Galerkin matrix for chi_4, ODD sector (sine basis) ===")
W_sin = build_galerkin_matrix_chi('sin', N_MODES, L, LAMBDA, primes, chi_values)
eigvals_sin = np.linalg.eigvalsh(W_sin)
print(f"[Sine sector] Smallest 5 eigenvalues: {eigvals_sin[:5]}")
print(f"[Sine sector] lambda_1^- = {eigvals_sin[0]:.6f}")

# -----------------------------------------------------------------------------
# 8. Parity-Dominance-Analyse
# -----------------------------------------------------------------------------

gap = eigvals_sin[0] - eigvals_cos[0]  # Odd minus Even
print(f"\n=== Parity Dominance Analysis ===")
print(f"lambda_1^+ (Even/cos) = {eigvals_cos[0]:.6f}")
print(f"lambda_1^- (Odd/sin)  = {eigvals_sin[0]:.6f}")
print(f"Gap = lambda_1^- - lambda_1^+ = {gap:.6f}")

if eigvals_sin[0] < eigvals_cos[0]:
    print("==> ODD-DOMINANCE bestaetigt (lambda_1^- < lambda_1^+)")
    print("    Vorhersage des Meta-Papers fuer odd character chi_4: VERIFIZIERT")
elif eigvals_cos[0] < eigvals_sin[0]:
    print("==> EVEN-Dominance (unerwartet fuer chi_4; Vorhersage widerlegt oder Bug)")
else:
    print("==> Kein klares Dominance-Ergebnis")

# Zum Vergleich: wie ist das Verhaeltnis?
print(f"\nRatio gap / |lambda_1^+| = {gap / abs(eigvals_cos[0]):.4f}")
print(f"Ratio gap / sqrt(lambda) = {gap / np.sqrt(LAMBDA):.4f}")
