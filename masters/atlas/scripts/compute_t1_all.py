#!/usr/bin/env python3
# coding: utf-8
"""
Berechne fuer alle getesteten Charaktere die niedrigste nicht-triviale Null-Stelle t_1
der L-Funktion L(s, chi). Fuer reelle Charaktere nutzen wir die Hardy-Z-artige
Funktion:  Z_chi(t) = exp(i theta(t)) * L(1/2 + i*t, chi)  ist reell,
und Null-Stellen von L sind Vorzeichenwechsel von Z_chi.

Wir nutzen die Vereinfachung: fuer reelle Dirichlet-Charaktere gilt
L(conj(s), chi) = conj(L(s, chi)),  also  L(1/2 + it, chi) ist "fast reell"
up to theta-Phase.

Konkreter: die Phase theta_chi(t) fuer den funktionalen-Gleichungs-Part
ist (per analytic conductor Theorie):
  theta_chi(t) = Im(logGamma((s + a_chi)/2)) - (t/2) log(q/pi)
mit a_chi = 0 (even) oder 1 (odd).

Wir plotten einfach |L(1/2+it)| und suchen Minima.
"""
import numpy as np
import math
from mpmath import mp, mpf, mpc, zeta
import mpmath
mp.dps = 30

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

def make_chi_D(D):
    def chi(n):
        if math.gcd(n, abs(D)) != 1: return 0
        return kronecker_symbol(D, n)
    return chi

def L_at(s, chi_fn, q):
    """L(s, chi) via Hurwitz-Zeta (gueltig fuer Re(s) > 0, s != 1)."""
    total = mpc(0)
    for a in range(1, q):
        ca = chi_fn(a)
        if ca == 0: continue
        total += ca * zeta(s, mpf(a)/q)
    return mpf(q)**(-s) * total

def find_first_zero(chi_fn, q, t_min=0.5, t_max=15.0, fine_step=0.02):
    """Finde niedrigstes t > t_min mit L(1/2+it, chi) = 0.
    Methode: scanne |L|, finde lokale Minima, refinere per mpmath.findroot auf L=0."""
    # Scan mit feinem Raster
    ts = np.arange(t_min, t_max, fine_step)
    vals = []
    for t in ts:
        try:
            v = L_at(mpf(1)/2 + 1j*mpf(t), chi_fn, q)
            vals.append(abs(v))
        except Exception:
            vals.append(1e10)
    vals = np.array([float(v) for v in vals])

    # Finde alle lokalen Minima unter Schwelle
    candidates = []
    for i in range(1, len(ts)-1):
        if vals[i] < vals[i-1] and vals[i] < vals[i+1] and vals[i] < 0.3:
            candidates.append((float(ts[i]), float(vals[i])))
    if not candidates:
        return None, min(vals)

    # Refinere erstes Minimum
    t0, v0 = candidates[0]
    # Numerische Suche auf L(1/2 + it) = 0 per Bisection auf dem Realteil * Imagteil-Ansatz:
    # Fuer reelle chi ist  L(1/2 + it, chi) = Re + i*Im mit einer gemeinsamen Phase.
    # Wir nutzen: Null wenn gleichzeitig Re und Im klein.
    # Suche mit feinstem Raster um Kandidat.
    fine_ts = np.arange(t0 - fine_step, t0 + fine_step, fine_step/20)
    fine_vals = [float(abs(L_at(mpf(1)/2 + 1j*mpf(t), chi_fn, q))) for t in fine_ts]
    i_min = int(np.argmin(fine_vals))
    t_best = float(fine_ts[i_min])
    v_best = fine_vals[i_min]
    # Noch feiner
    fine2 = np.arange(t_best - fine_step/20, t_best + fine_step/20, fine_step/500)
    if len(fine2) > 0:
        fine2_vals = [float(abs(L_at(mpf(1)/2 + 1j*mpf(t), chi_fn, q))) for t in fine2]
        i2 = int(np.argmin(fine2_vals))
        t_best = float(fine2[i2])
        v_best = fine2_vals[i2]
    return t_best, v_best

# Alle getesteten Charaktere (reell, primitiv, even)
test_chars = [
    ('chi_5',  5,  'Legendre mod 5'),
    ('chi_8',  8,  'Kronecker 8 (= primitiv mod 8)'),
    ('chi_12', 12, 'Kronecker 12 (= primitiv mod 12)'),
    ('chi_13', 13, 'Legendre mod 13'),
    ('chi_17', 17, 'Legendre mod 17'),
    ('chi_21', 21, 'Kronecker 21 = 3*7'),
    ('chi_24', 24, 'Kronecker 24 = 8*3'),
    ('chi_29', 29, 'Legendre mod 29'),
    ('chi_33', 33, 'Kronecker 33 = 3*11'),
    ('chi_60', 60, 'Kronecker 60 = 4*3*5'),
]

# Gemessene Gap-Werte (mean ueber getestetes lambda-Grid bei N=200)
# aus Session 4 Daten
mean_gaps_N200 = {
    'chi_5':  0.160,
    'chi_8':  0.055,
    'chi_12': 0.131,
    'chi_13': 0.271,
    'chi_17': 0.026,
    'chi_21': 0.094,
    'chi_24': 0.001,
    'chi_29': 0.132,
    'chi_33': -0.011,
    'chi_60': 0.157,
}

# Stabilitaets-Flag (alle pos bei N=200)
pos_fraction_N200 = {
    'chi_5':  5/8, 'chi_8': 5/8, 'chi_12': 8/8, 'chi_13': 6/8,
    'chi_17': 7/8, 'chi_21': 5/8, 'chi_24': 5/8, 'chi_29': 7/8,
    'chi_33': 4/8, 'chi_60': 6/8,
}

# Berechne auch L(1) und Regulator-Info
def L_at_1(chi_fn, q):
    total = mpc(0)
    for a in range(1, q):
        ca = chi_fn(a)
        if ca == 0: continue
        total += ca * mpmath.digamma(mpf(a)/q)
    return -total / mpf(q)

def Lprime_at_1(chi_fn, q):
    L1 = L_at_1(chi_fn, q)
    total = mpc(0)
    for a in range(1, q):
        ca = chi_fn(a)
        if ca == 0: continue
        total += ca * mpmath.loggamma(mpf(a)/q)
    return -mpmath.log(q) * L1 - total / mpf(q)

results = []
H_LOGRATIO = "-Lprime/L"
print(f"{'Char':>8}  {'D':>3}  {'t_1':>8}  {'|L|':>8}  "
      f"{'L(1)':>8}  {H_LOGRATIO:>10}  {'gap':>8}  {'pos_frac':>8}")
print("-" * 80)
for name, D, label in test_chars:
    chi = make_chi_D(D)
    # Sanity: chi(-1) muss +1 sein (even)
    if chi(-1) != 1:
        print(f"{name:>8}  {D:>3}  ODD  (skip)")
        continue

    t1, val = find_first_zero(chi, D)
    L1 = L_at_1(chi, D)
    Lp1 = Lprime_at_1(chi, D)
    ratio = float((-Lp1/L1).real) if isinstance(Lp1, mpc) else float(-Lp1/L1)
    L1_val = float(L1.real) if isinstance(L1, mpc) else float(L1)
    gap = mean_gaps_N200[name]
    pfrac = pos_fraction_N200[name]

    results.append({
        'name': name, 'D': D, 't1': t1, 'L1': L1_val, 'logL_prime_over_L': ratio,
        'gap': gap, 'pos_frac': pfrac, 'log_D': math.log(D),
    })
    t1_str = f"{t1:>8.4f}" if t1 else '     ---'
    val_str = f"{val:>8.4f}" if val else '     ---'
    print(f"{name:>8}  {D:>3}  {t1_str}  {val_str}  "
          f"{L1_val:>+8.4f}  {ratio:>+8.4f}  {gap:>+8.4f}  {pfrac:>8.3f}")

import json
with open('../_results/char_features_all.json', 'w') as f:
    json.dump(results, f, indent=2)

# ========================================================================
# MULTI-VARIATE REGRESSION
# ========================================================================
print()
print("="*70)
print("=== MULTI-VARIATE REGRESSION gap ~ f(t_1, L1, -L'/L, log D) ===")
print("="*70)

# Filtere Zeilen mit t_1 vorhanden
valid = [r for r in results if r['t1'] is not None]
print(f"\nN = {len(valid)} charactere")

# Bau Feature-Matrix
import numpy as np
X = np.array([[r['t1'], r['L1'], r['logL_prime_over_L'], r['log_D']] for r in valid])
y = np.array([r['gap'] for r in valid])
names = [r['name'] for r in valid]
t1s = np.array([r['t1'] for r in valid])

# Einzelne Regressionen (univariate)
print(f"\n{'Feature':>20}  {'Slope':>10}  {'Intercept':>10}  {'R^2':>8}")
features = ['t_1', 'L(1)', H_LOGRATIO, 'log D']
for i, feat in enumerate(features):
    x = X[:, i]
    # Lineare Regression
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1 - ss_res/ss_tot
    print(f"{feat:>20}  {slope:>+10.5f}  {intercept:>+10.5f}  {r2:>8.4f}")

# log-linear (Expo-Fit): log|gap| ~ a*t_1 + b
print(f"\n{'Expo-Fit log|gap| ~ a*t_1 + b':>30}")
mask = y > 1e-3
x_t1 = t1s[mask]
y_log = np.log(y[mask])
slope, intercept = np.polyfit(x_t1, y_log, 1)
y_pred = slope * x_t1 + intercept
r2 = 1 - np.sum((y_log - y_pred)**2) / np.sum((y_log - y_log.mean())**2)
print(f"  Slope = {slope:+.4f}  Intercept = {intercept:+.4f}  R^2 = {r2:.4f}")
print(f"  d.h. gap ~ exp({intercept:+.4f}) * exp({slope:+.4f} * t_1)"
      f" = {math.exp(intercept):.4f} * exp({slope:+.4f} * t_1)")

# Multi-variate Regression: gap = a*t_1 + b*log(D) + c*L(1) + d*(-L'/L) + e
print(f"\nMulti-variate Regression (OLS):")
# X_aug mit Intercept
X_aug = np.column_stack([X, np.ones(len(X))])
coef, residuals, rank, sv = np.linalg.lstsq(X_aug, y, rcond=None)
y_pred_full = X_aug @ coef
r2_full = 1 - np.sum((y - y_pred_full)**2) / np.sum((y - y.mean())**2)
print(f"  Koeffizienten:")
for feat, c in zip(features + ['const'], coef):
    print(f"    {feat:>10}: {c:+.5f}")
print(f"  R^2 = {r2_full:.4f}")

# Stepwise: bestes Ein-Feature-Modell
print(f"\nBeste Ein-Feature-Modelle (Rank nach R^2):")
r2s = []
for i, feat in enumerate(features):
    x = X[:, i]
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1 - ss_res/ss_tot
    r2s.append((r2, feat, slope, intercept))
r2s.sort(reverse=True)
for r2, feat, slope, intercept in r2s:
    print(f"  {feat:>10}: gap = {slope:+.4f} * x + {intercept:+.4f}, R^2 = {r2:.4f}")

# Darstellung der Daten
print(f"\nDaten sortiert nach t_1:")
print(f"{'Name':>8}  {'t_1':>8}  {'gap':>8}  {'log|gap|':>10}")
sorted_r = sorted(valid, key=lambda r: r['t1'])
for r in sorted_r:
    if r['gap'] > 1e-4:
        lg = math.log(r['gap']) if r['gap'] > 0 else None
        lg_s = f"{lg:>+10.3f}" if lg else "   (neg)"
    else:
        lg_s = "  (zero)"
    print(f"{r['name']:>8}  {r['t1']:>8.4f}  {r['gap']:>+8.4f}  {lg_s}")

print("\n[done]")
