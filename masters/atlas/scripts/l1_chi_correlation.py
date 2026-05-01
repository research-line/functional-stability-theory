#!/usr/bin/env python3
# coding: utf-8
"""
Teste Hypothese aus THEORIE_CHI12_KONSTANZ.md:
  Gap_chi(infty) ~ C * L'(1, chi) / L(1, chi)

Berechnet L(1, chi) und L'(1, chi) fuer mehrere primitive reelle Charaktere
via sympy.mpmath und vergleicht mit gemessenen Gap-Werten aus Session 4.

L(s, chi) = sum_{n>=1} chi(n) / n^s
L'(s, chi) = -sum_{n>=1} chi(n) * log(n) / n^s

Wir nutzen asymptotische Reihen-Summation mit hoher Genauigkeit.
"""
import numpy as np
from mpmath import mpf, mpc, mp, log as mplog
import mpmath

mp.dps = 50  # 50 Dezimalstellen Genauigkeit

# Charaktere
def chi5(n):
    if n % 5 == 0: return 0
    return {1: 1, 2: -1, 3: -1, 4: 1}[n % 5]

def chi8a(n):
    if n % 2 == 0: return 0
    return {1: 1, 3: -1, 5: -1, 7: 1}[n % 8]

def chi12(n):
    import math
    if math.gcd(n, 12) != 1: return 0
    return {1: 1, 5: -1, 7: -1, 11: 1}[n % 12]

def chi13(n):
    if n % 13 == 0: return 0
    qr = {pow(i, 2, 13) for i in range(1, 13)}
    return 1 if (n % 13) in qr else -1

def chi4(n):
    # Zum Vergleich (odd character)
    if n % 2 == 0: return 0
    return 1 if (n % 4) == 1 else -1

characters = [
    ('chi_5',  chi5,  5,  'Legendre mod 5'),
    ('chi_8a', chi8a, 8,  'primitive mod 8'),
    ('chi_12', chi12, 12, 'primitive mod 12'),
    ('chi_13', chi13, 13, 'Legendre mod 13'),
    ('chi_4',  chi4,  4,  'Legendre mod 4 (odd)'),
]

# Gemessene Gap-Werte (mean bei N=200 aus Session 4)
measured_gap = {
    'chi_5':  None,    # oszillierend, kein einheitlicher Mean
    'chi_8a': None,    # wird nach Server-Lauf bekannt
    'chi_12': mpf('0.127'),  # konstant
    'chi_13': None,    # wird nach Server-Lauf bekannt
    'chi_4':  None,    # odd, kein klarer Gap
}

# ============================================================
# Via Dirichlet L-Funktion direkt in mpmath
# ============================================================

def L_at_1(chi_fn, q):
    """Gauss-Formel fuer primitive nicht-triviale chi mod q:
       L(1, chi) = -(1/q) sum_{a=1}^{q-1} chi(a) * psi(a/q)
    (wobei psi = digamma, Pol-Kompensation wegen sum chi(a) = 0)"""
    total = mpc(0)
    for a in range(1, q):
        ca = chi_fn(a)
        if ca == 0: continue
        total += ca * mpmath.digamma(mpf(a)/q)
    return -total / mpf(q)

def Lprime_at_1(chi_fn, q):
    """L'(1, chi) via Hurwitz-Zeta-Expansion bei s=1.
    Ableitung der Beziehung L(s,chi) = q^{-s} sum_a chi(a) zeta(s, a/q):
       L'(1, chi) = -log(q) * L(1, chi) - (1/q) sum chi(a) log Gamma(a/q)
    (Lerch-Formel: zeta'(1, x) = -log Gamma(x) + (1/2) log(2pi), der (1/2) log(2pi)-Term
     verschwindet wegen sum chi(a) = 0.)"""
    L1 = L_at_1(chi_fn, q)
    total = mpc(0)
    for a in range(1, q):
        ca = chi_fn(a)
        if ca == 0: continue
        total += ca * mpmath.loggamma(mpf(a)/q)
    return -mpmath.log(q) * L1 - total / mpf(q)

# ============================================================
# Hauptrechnung
# ============================================================

header_ratio = "-Lprime/L"
print(f"{'Charakter':>10}  {'L(1, chi)':>14}  {'L_prime(1, chi)':>16}  "
      f"{header_ratio:>12}  {'Gemessener Gap':>16}  {'Verhaeltnis':>12}")
print("-" * 90)

results = []
for name, chi_fn, q, label in characters:
    try:
        L1 = L_at_1(chi_fn, q)
        Lp1 = Lprime_at_1(chi_fn, q)
        ratio = -Lp1 / L1
        ratio_real = float(ratio.real) if isinstance(ratio, mpc) else float(ratio)
        gap = measured_gap[name]
        if gap is not None:
            proportionality = float(gap) / ratio_real if ratio_real != 0 else None
            gap_str = f"{float(gap):>+14.6f}"
            prop_str = f"{proportionality:>+12.4f}" if proportionality else "---"
        else:
            gap_str = "    (pending)   "
            prop_str = "          ---"
        L1_str = f"{float(L1.real):>+14.8f}" if isinstance(L1, mpc) else f"{float(L1):>+14.8f}"
        Lp1_str = f"{float(Lp1.real):>+16.8f}" if isinstance(Lp1, mpc) else f"{float(Lp1):>+16.8f}"
        print(f"{name:>10}  {L1_str}  {Lp1_str}  {ratio_real:>+12.6f}  {gap_str}  {prop_str}")
        results.append({
            'name': name, 'q': q,
            'L1': float(L1.real) if isinstance(L1, mpc) else float(L1),
            'Lp1': float(Lp1.real) if isinstance(Lp1, mpc) else float(Lp1),
            'ratio': ratio_real,
            'gap': float(gap) if gap is not None else None,
        })
    except Exception as e:
        print(f"{name:>10}  ERROR: {e}")

print()
print("=== Interpretation ===")
print()
print("Falls Hypothese 'Gap = C * (-L'/L)' stimmt, sollte 'Verhaeltnis' konstant sein.")
print("Aktuell nur chi_12 gemessen; nach Server-Lauf (four_even_results.json) koennen")
print("chi_5, chi_8a, chi_13 ergaenzt werden.")
print()
print("=== Zusatz: weitere primitive reelle even Kandidaten ===")
print()

def kronecker_symbol(a, n):
    """Elementare Implementation des Kronecker-Symbols (a/n)."""
    import math
    if n == 0:
        return 1 if abs(a) == 1 else 0
    if n < 0:
        return kronecker_symbol(a, -n) * (1 if a >= 0 else -1)
    # Faktorisiere 2 aus n
    result = 1
    while n % 2 == 0:
        if a % 2 == 0:
            return 0
        if a % 8 in (3, 5):
            result = -result
        n //= 2
    # n ungerade, > 0
    # Jetzt Jacobi-Symbol (a/n)
    a = a % n
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a = a % n
    if n == 1:
        return result
    return 0

def chi_kronecker(D):
    """Konstruiere Kronecker-Charakter chi_D fuer fundamentale Diskriminante D."""
    def chi(n):
        import math
        if math.gcd(n, abs(D)) != 1: return 0
        return kronecker_symbol(D, n)
    return chi

# Fundamentale Diskriminanten, alle EVEN characters (D > 0)
# D = 5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37, 40, 41, 44, 53, 56, 57, 60, ...
D_candidates = [5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37, 40, 41, 44, 57, 60]

print(f"{'D':>4}  {'L(1, chi_D)':>14}  {'L_prime(1)':>14}  {header_ratio:>12}  {'log|D|':>8}")
for D in D_candidates:
    chi_fn = chi_kronecker(D)
    q = abs(D)
    try:
        L1 = L_at_1(chi_fn, q)
        Lp1 = Lprime_at_1(chi_fn, q)
        ratio = -Lp1 / L1
        ratio_real = float(ratio.real) if isinstance(ratio, mpc) else float(ratio)
        L1_val = float(L1.real) if isinstance(L1, mpc) else float(L1)
        Lp1_val = float(Lp1.real) if isinstance(Lp1, mpc) else float(Lp1)
        print(f"{D:>4}  {L1_val:>+14.6f}  {Lp1_val:>+14.6f}  {ratio_real:>+12.6f}  {np.log(D):>8.3f}")
    except Exception as e:
        print(f"{D:>4}  ERROR: {e}")

print()
print("[done]")
