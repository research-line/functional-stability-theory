#!/usr/bin/env python3
# coding: utf-8
"""
Task #38: chi_33-Anomalie untersuchen.

Kernfrage: Was unterscheidet die 3 Charaktere mit negativen asymptotischen
Gaps (chi_8, chi_13, chi_33) von den 4 mit positiven (chi_5, chi_12,
chi_17, chi_24) bei lam=20000?

Hypothesen:
1. Chebyshev-Bias: Verteilung von chi(p) = +1 vs -1 ueber alle p <= lam.
2. Primsum T_chi(lam) = Sum_{p<=lam} chi(p) log(p) / sqrt(p).
3. Niedrige-Primzahlen-Dominanz: chi-Vorzeichen bei p <= 100.
4. Euler-Produkt-Analog: Produkt(1 - chi(p)/sqrt(p)) vs anti.
5. Faktorisierung des Moduls: chi_33 = chi_3 * chi_11 (Dirichlet-Produkt).
"""
import math
import sympy
import numpy as np

def kronecker_symbol(a, n):
    if n == 0: return 1 if abs(a) == 1 else 0
    if n < 0: return kronecker_symbol(a, -n) * (1 if a >= 0 else -1)
    result = 1
    while n % 2 == 0:
        if a % 2 == 0: return 0
        if a % 8 in (3, 5): result = -result
        n //= 2
    a = a % n
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: result = -result
        a = a % n
    return result if n == 1 else 0

def make_chi(D):
    def chi(n):
        if math.gcd(n, abs(D)) != 1: return 0
        return kronecker_symbol(D, n)
    return chi

# Die 10 Charaktere mit gemessenem Gap bei lam=20000, N=konvergiert
# (N=400 falls nichts besseres, N=600 fuer die 6 oszillierenden)
chars = [
    ('chi_5',  5,  +0.0156),   # N=400 (aus HIGH_N)
    ('chi_8',  8,  -0.04902),  # N=600
    ('chi_12', 12, +0.0337),   # N=400 (schon N-konvergent aus HIGH_N)
    ('chi_13', 13, -0.12504),  # N=600
    ('chi_17', 17, +0.01318),  # N=400
    ('chi_21', 21, -0.00423),  # N=600
    ('chi_24', 24, +0.01076),  # N=400
    ('chi_29', 29, +0.01439),  # N=600
    ('chi_33', 33, -0.14221),  # N=600
    ('chi_60', 60, +0.00521),  # N=600
]

# ---------------------------------------------------------------
# Analyse 1: Primsum T_chi(lam) fuer lam=20000 und Chebyshev-Bias
# ---------------------------------------------------------------
print("="*110)
print("ANALYSE 1: Primsum und Chebyshev-Bias bei lam=20000")
print("="*110)
lam = 20000
primes_up_to_lam = list(sympy.primerange(2, lam+1))
print(f"Anzahl Primzahlen bis {lam}: {len(primes_up_to_lam)}")
print()

print(f"{'Char':>8}  {'D':>3}  {'gap(lam=2e4)':>13}  {'T_chi(lam)':>11}  "
      f"{'#pos':>5}  {'#neg':>5}  {'Bias neg-pos':>13}  {'cov_primes':>10}")

rows = []
for name, D, gap_asymp in chars:
    chi = make_chi(D)
    chi_vals = [chi(p) for p in primes_up_to_lam]
    # Primsum (im "v2.1"-Sinn, m=1)
    T = sum(cp * math.log(p) / math.sqrt(p)
            for p, cp in zip(primes_up_to_lam, chi_vals))
    pos = sum(1 for v in chi_vals if v == 1)
    neg = sum(1 for v in chi_vals if v == -1)
    zero = sum(1 for v in chi_vals if v == 0)
    bias = neg - pos
    cov = (pos + neg) / len(primes_up_to_lam)  # Anteil nicht-ausgeschlossener Primzahlen
    print(f"{name:>8}  {D:>3}  {gap_asymp:>+13.5f}  {T:>+11.4f}  "
          f"{pos:>5}  {neg:>5}  {bias:>+13d}  {cov:>10.3f}")
    rows.append({'name': name, 'D': D, 'gap': gap_asymp, 'T': T,
                 'pos': pos, 'neg': neg, 'zero': zero, 'bias': bias, 'cov': cov})

# ---------------------------------------------------------------
# Analyse 2: Split in negative vs. positive Charaktere
# ---------------------------------------------------------------
print()
print("="*110)
print("ANALYSE 2: Unterscheidet sich das Primsum-Muster?")
print("="*110)
neg_chars = [r for r in rows if r['gap'] < -0.01]
pos_chars = [r for r in rows if r['gap'] > +0.01]
near_zero = [r for r in rows if abs(r['gap']) <= 0.01]
print(f"\nNegative asymptotische Gaps (3): {[r['name'] for r in neg_chars]}")
print(f"Positive asymptotische Gaps (3): {[r['name'] for r in pos_chars]}")
print(f"Nah-Null (4): {[r['name'] for r in near_zero]}")

def stat(group, key):
    vals = [r[key] for r in group]
    if not vals: return "---"
    return f"mean={np.mean(vals):+.3f}, std={np.std(vals):.3f}, min={min(vals):+.3f}, max={max(vals):+.3f}"

print(f"\nT_chi(lam=20000):")
print(f"  negative gap chars: {stat(neg_chars, 'T')}")
print(f"  positive gap chars: {stat(pos_chars, 'T')}")
print(f"  near-zero chars:    {stat(near_zero, 'T')}")

print(f"\nChebyshev-Bias (#neg - #pos) bei lam=20000:")
print(f"  negative gap chars: {stat(neg_chars, 'bias')}")
print(f"  positive gap chars: {stat(pos_chars, 'bias')}")
print(f"  near-zero chars:    {stat(near_zero, 'bias')}")

# Korrelation T_chi vs gap
print(f"\nKorrelation T_chi(lam=2e4) vs gap_asymp: {np.corrcoef([r['T'] for r in rows], [r['gap'] for r in rows])[0,1]:+.3f}")
print(f"Korrelation bias vs gap_asymp: {np.corrcoef([r['bias'] for r in rows], [r['gap'] for r in rows])[0,1]:+.3f}")

# Ratio gap/T
print(f"\nRatio gap/T_chi(lam=2e4) je Charakter:")
for r in rows:
    if abs(r['T']) > 0.01:
        ratio = r['gap']/r['T']
        print(f"  {r['name']:>8}  T={r['T']:+.4f}  gap={r['gap']:+.5f}  ratio={ratio:+.5f}")

# ---------------------------------------------------------------
# Analyse 3: Niedrige Primzahlen (erste 20)
# ---------------------------------------------------------------
print()
print("="*110)
print("ANALYSE 3: Chi-Vorzeichen bei niedrigen Primzahlen (erste 20)")
print("="*110)
first20 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
print(f"{'Char':>8}  " + "  ".join(f"p={p:>2}" for p in first20))
for name, D, gap_asymp in chars:
    chi = make_chi(D)
    row = f"{name:>8}  "
    row += "  ".join(
        ("  0 " if chi(p) == 0 else (" +1 " if chi(p) == 1 else " -1 "))
        for p in first20
    )
    print(row)

# Berechne "signed-sqrt-bias" = Sum_{p in first50, chi(p) != 0} chi(p)/sqrt(p)
# (ohne log, um einfacher zu machen)
print()
print("Signed-sqrt-bias = Sum_{first 50 primes} chi(p)/sqrt(p):")
first50 = list(sympy.primerange(2, 229+1))  # 50 Primzahlen
for name, D, gap_asymp in chars:
    chi = make_chi(D)
    s = sum(chi(p)/math.sqrt(p) for p in first50)
    # Gewichtet: -log(p)/sqrt(p)
    s_w = sum(chi(p)*math.log(p)/math.sqrt(p) for p in first50)
    print(f"  {name:>8}  signed_sqrt={s:+.4f}  signed_weighted={s_w:+.4f}  gap={gap_asymp:+.5f}")

# ---------------------------------------------------------------
# Analyse 4: chi_33 Faktorisierung (chi_33 = chi_{-3} * chi_{-11}?)
# ---------------------------------------------------------------
print()
print("="*110)
print("ANALYSE 4: chi_33 Faktorisierung")
print("="*110)
# 33 = 3 * 11. Beide sind 3 mod 4, also D=-3 und D=-11 sind odd fundamentale Diskriminanten.
# Produkt (-3)*(-11) = 33 ist Even fundamentale Diskriminante.
# Also: chi_33(n) = chi_{-3}(n) * chi_{-11}(n) als primitive Charaktere.
chi_33 = make_chi(33)
chi_m3 = make_chi(-3)
chi_m11 = make_chi(-11)
# Check: chi_33(n) = chi_{-3}(n) * chi_{-11}(n)?
print("Check: chi_33(n) vs chi_{-3}(n)*chi_{-11}(n) fuer n=1..40")
ok = True
for n in range(1, 41):
    a = chi_33(n)
    b = chi_m3(n) * chi_m11(n)
    if a != b:
        print(f"  n={n}: chi_33={a} vs chi_{{-3}}*chi_{{-11}}={b}  MISMATCH")
        ok = False
print(f"  Match: {ok}")
print()

# Berechne Primsummen fuer chi_{-3}, chi_{-11} separat
for name, fn in [('chi_{-3}', chi_m3), ('chi_{-11}', chi_m11), ('chi_33', chi_33)]:
    if fn(-1) == -1:
        parity = "odd"
    elif fn(-1) == 1:
        parity = "even"
    else:
        parity = "?"
    T = sum(fn(p) * math.log(p) / math.sqrt(p) for p in primes_up_to_lam)
    pos = sum(1 for p in primes_up_to_lam if fn(p) == 1)
    neg = sum(1 for p in primes_up_to_lam if fn(p) == -1)
    print(f"  {name:>12} ({parity})  T(lam={lam})={T:+.4f}  pos={pos}  neg={neg}  bias={neg-pos:+d}")

print()
print("[done]")
