#!/usr/bin/env python3
# coding: utf-8
"""
Berechne die ersten nicht-trivialen Null-Stellen von L(s, chi) fuer
verschiedene reelle primitive Dirichlet-Charaktere.

Hypothese: Charaktere mit besonders HOHER erster nicht-trivialer Null-Stelle
(d.h. niedrigliegende Nullstelle weit oben) koennten stabilere Gap-Asymptotik
zeigen (weil der "fuehrende Oszillations-Term" in der Explicit Formula dann
weiter oben liegt).

Via mpmath: finde Nullen von L(1/2 + it, chi) per numerischer Wurzelsuche.
"""
import numpy as np
import math
from mpmath import mp, mpf, mpc, zeta, digamma
mp.dps = 30

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

def L_via_hurwitz(s, chi_fn, q):
    """L(s, chi) = q^{-s} sum_{a=1}^{q-1} chi(a) zeta(s, a/q)
    Vorsicht bei s=1 (Pol-Kompensation), aber fuer s = 1/2 + it gut definiert."""
    total = mpc(0)
    for a in range(1, q):
        ca = chi_fn(a)
        if ca == 0: continue
        total += ca * zeta(s, mpf(a)/q)
    return mpf(q)**(-s) * total

def find_first_zero(chi_fn, q, t_max=50, n_points=500):
    """Suche die erste Null-Stelle von t -> |L(1/2 + i*t, chi)| im Bereich [0, t_max].

    Scanne grob, dann refiniere per mpmath.findroot wenn Vorzeichenwechsel oder Minimum.
    """
    import mpmath
    # Scanne |L|^2 - sehr klein -> Null-Kandidat
    ts = np.linspace(0.5, t_max, n_points)
    abs_vals = []
    for t in ts:
        try:
            val = L_via_hurwitz(mpf(1)/2 + 1j*mpf(t), chi_fn, q)
            abs_vals.append(float(abs(val)))
        except:
            abs_vals.append(1e10)
    abs_vals = np.array(abs_vals)
    # Lokale Minima finden
    minima = []
    for i in range(1, len(ts)-1):
        if abs_vals[i] < abs_vals[i-1] and abs_vals[i] < abs_vals[i+1] and abs_vals[i] < 0.1:
            minima.append((ts[i], abs_vals[i]))
    # Refinieren per findroot um das erste Minimum
    if not minima:
        return None, abs_vals.min()
    t0, v0 = minima[0]
    try:
        root = mpmath.findroot(
            lambda t: L_via_hurwitz(mpf(1)/2 + 1j*t, chi_fn, q),
            mpc(t0, 0),
            tol=mpf('1e-15')
        )
        return float(root.imag), float(abs(L_via_hurwitz(mpf(1)/2 + 1j*root.imag, chi_fn, q)))
    except Exception as e:
        return t0, v0

characters = [
    ('chi_5', chi5, 5),
    ('chi_8a', chi8a, 8),
    ('chi_12', chi12, 12),
    ('chi_13', chi13, 13),
    ('chi_17', chi17, 17),
    ('chi_29', chi29, 29),
]

# Gemessene Mean-Gaps (aus Session 4)
mean_gaps = {
    'chi_5': 0.160,
    'chi_8a': 0.055,
    'chi_12': 0.131,
    'chi_13': 0.271,
    'chi_17': 0.026,
    'chi_29': 0.132,
}

# Stabilitaets-Flag (aus Session 4 Befunden)
stable = {
    'chi_5': False,
    'chi_8a': False,
    'chi_12': True,    # 8/8, slope ~ 0
    'chi_13': False,
    'chi_17': True,    # 7/8, slope ~ 0, aber kleiner Mean
    'chi_29': False,
}

print(f"{'Char':>8}  {'q':>3}  {'erste Nullst. t_1':>18}  "
      f"{'|L|':>10}  {'mean gap':>8}  {'stabil?':>8}")
print("-" * 68)
results = []
for name, fn, q in characters:
    t1, val = find_first_zero(fn, q, t_max=30, n_points=600)
    if t1 is None:
        print(f"{name:>8}  {q:>3}  {'NICHT GEFUNDEN':>18}  {'(large?)':>10}  "
              f"{mean_gaps[name]:>+8.3f}  {'STABIL' if stable[name] else 'osz.':>8}")
    else:
        print(f"{name:>8}  {q:>3}  {t1:>18.6f}  {val:>10.6f}  "
              f"{mean_gaps[name]:>+8.3f}  {'STABIL' if stable[name] else 'osz.':>8}")
        results.append({'name': name, 'q': q, 't1': t1, 'mean_gap': mean_gaps[name],
                        'stable': stable[name]})

print()
print("=== Interpretation ===")
print()
print("Hypothese: hoehere t_1 -> stabileres Gap?")
sorted_by_t1 = sorted(results, key=lambda r: r['t1'])
for r in sorted_by_t1:
    print(f"  t_1 = {r['t1']:.2f}  {r['name']}  "
          f"mean_gap = {r['mean_gap']:+.3f}  "
          f"{'STABIL' if r['stable'] else 'osz.'}")

print()
print("Hypothese Test: Rangkorrelation t_1 <-> Stabilitaet?")
ts = [r['t1'] for r in results]
stabs = [1 if r['stable'] else 0 for r in results]
# Spearman-Korrelation
import scipy.stats
rs, pval = scipy.stats.spearmanr(ts, stabs)
print(f"Spearman rho (t_1, stabil) = {rs:+.4f}, p = {pval:.3f}")
print(f"(N = {len(results)})")
