#!/usr/bin/env python3
# coding: utf-8
"""
Vergleiche den gemessenen Gap mit der Primsum-Struktur, die in der
Weil-Form direkt einfließt:

  T_chi(lam) := sum_{p <= lam} chi(p) log(p) / sqrt(p)

Die Hypothese (grob): gap_chi(lam) ~ C * T_chi(lam) bis auf Truncation.

Interessant ist der Unterschied:
- chi_12: T_chi(50000) sollte Konvergenz-Wert haben (nah einer Konstanten).
- chi_5: T_chi oszilliert staerker.

Wir vergleichen T_chi(lam) fuer lam in der gleichen Grid wie Session 4.
"""
import numpy as np
import sympy
import math
import sys
sys.stdout.reconfigure(line_buffering=True)

# Charaktere
def chi5(n):
    if n % 5 == 0: return 0
    return {1: 1, 2: -1, 3: -1, 4: 1}[n % 5]

def chi8a(n):
    if n % 2 == 0: return 0
    return {1: 1, 3: -1, 5: -1, 7: 1}[n % 8]

def chi12(n):
    if math.gcd(n, 12) != 1: return 0
    return {1: 1, 5: -1, 7: -1, 11: 1}[n % 12]

def chi13(n):
    if n % 13 == 0: return 0
    qr = {pow(i, 2, 13) for i in range(1, 13)}
    return 1 if (n % 13) in qr else -1

def chi17(n):
    if n % 17 == 0: return 0
    qr = {pow(i, 2, 17) for i in range(1, 17)}
    return 1 if (n % 17) in qr else -1

def chi29(n):
    if n % 29 == 0: return 0
    qr = {pow(i, 2, 29) for i in range(1, 29)}
    return 1 if (n % 29) in qr else -1

def T_chi(chi_fn, lam):
    """T_chi(lam) = sum_{p <= lam} chi(p) log(p) / sqrt(p)"""
    s = 0.0
    for p in sympy.primerange(2, lam+1):
        cp = chi_fn(p)
        if cp == 0: continue
        s += cp * math.log(p) / math.sqrt(p)
    return s

def T_chi_higher(chi_fn, lam):
    """Mit hoeheren Primpotenzen p^m:
       T_chi^full(lam) = sum_{p,m>=1: p^m <= lam} chi(p)^m log(p) / p^{m/2}
    Das ist der direkte Input in die v2.1 Weil-Matrix."""
    s = 0.0
    for p in sympy.primerange(2, lam+1):
        cp = chi_fn(p)
        if cp == 0: continue
        m = 1
        while p**m <= lam:
            s += (cp**m) * math.log(p) / (p**(m/2.0))
            m += 1
    return s

characters = [
    ('chi_5', chi5),
    ('chi_8a', chi8a),
    ('chi_12', chi12),
    ('chi_13', chi13),
    ('chi_17', chi17),
    ('chi_29', chi29),
]

# Gemessene Gap-Werte bei N=200 aus Session 4
gap_data = {
    'chi_5': {100: 1.5624, 200: -0.0035, 500: 0.0018, 1000: -0.3239, 2000: 0.0541,
              5000: 0.0567, 10000: -0.0836, 20000: 0.0143},
    'chi_8a': {100: 0.2025, 200: 0.2131, 500: 0.0116, 1000: -0.1300, 2000: -0.0546,
               5000: 0.1089, 10000: 0.1369, 20000: -0.0455},
    'chi_12': {100: 0.0248, 200: 0.1500, 500: 0.1178, 1000: 0.2332, 2000: 0.1692,
               5000: 0.1830, 10000: 0.1351, 20000: 0.0325},
    'chi_13': {100: 1.9235, 200: 0.0296, 500: 0.2275, 1000: -0.1507, 2000: 0.0672,
               5000: 0.0532, 10000: 0.1413},
    'chi_17': {100: -0.0098, 200: 0.0148, 500: 0.0251, 1000: 0.0219, 2000: 0.1285,
               5000: 0.0082, 10000: 0.0163, 20000: 0.0041},
    'chi_29': {100: 0.6685, 200: -0.0109, 500: 0.1045, 1000: 0.0319, 2000: 0.0053,
               5000: -0.0114, 10000: -0.0118, 20000: -0.1405},
}

lam_values = [100, 200, 500, 1000, 2000, 5000, 10000, 20000]

print(f"\n{'='*75}")
print("=== Primsum T_chi(lam) vergleiche zu Gap ===")
print(f"{'='*75}")
for name, fn in characters:
    print(f"\n--- {name} ---")
    print(f"{'lam':>6}  {'T_simple':>10}  {'T_full':>10}  {'gap':>10}  {'gap/T_full':>12}")
    for lam in lam_values:
        gap = gap_data[name].get(lam)
        if gap is None: continue
        t_simple = T_chi(fn, lam)
        t_full = T_chi_higher(fn, lam)
        ratio = gap / t_full if abs(t_full) > 1e-6 else None
        ratio_str = f"{ratio:+12.4f}" if ratio is not None else "        ---"
        print(f"{lam:>6}  {t_simple:>+10.5f}  {t_full:>+10.5f}  {gap:>+10.5f}  {ratio_str}")

print(f"\n{'='*75}")
print("=== Asymptotische Primsumme (lam = 20000) ===")
print(f"{'='*75}")
print(f"{'Char':>8}  {'T_simple(20000)':>18}  {'T_full(20000)':>18}  {'mean gap':>10}")
for name, fn in characters:
    t_s = T_chi(fn, 20000)
    t_f = T_chi_higher(fn, 20000)
    gaps = list(gap_data[name].values())
    print(f"{name:>8}  {t_s:>+18.5f}  {t_f:>+18.5f}  {np.mean(gaps):>+10.4f}")

print("\n[done]")
