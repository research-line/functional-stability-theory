#!/usr/bin/env python3
# coding: utf-8
"""
Task #39: Explicit-Formula-Experiment.

Plan:
1. Berechne fuer alle 10 Charaktere die ersten ~6 nicht-trivialen Null-
   Stellen t_k von L(1/2+it, chi) via mpmath (Minima von |L|).
2. Teste mehrere Test-Funktionen F_lam(t):
   (a) 1/(1/4 + t^2)       - Realteil von -1/rho fuer rho=1/2+it
   (b) cos(L*t)/(1/4+t^2)  - Weil-artig mit Frontier-Kern
   (c) exp(-t^2/(2*sigma^2)) - Gaussian mit Breite ~ 1/L
   (d) (sin(L*t/2) / (t/2))^2 - Fejer-Kern (sinc^2)
   fuer lambda=20000, L = log(lambda) ≈ 9.9.
3. Fuer jede Test-Funktion: berechne Sum_k F(t_k) (mit k fuer Paar k, -k)
   und vergleiche mit gemessenem Gap.
4. Beste Test-Funktion = groesste |Korrelation| und beste Vorzeichen-
   Vorhersage.

Laufzeit: ~10 Min pro Charakter fuer 6 Null-Stellen bei mpmath dps=20.
Total ca. 1.5 h.
"""
import math
import time
import numpy as np
from mpmath import mp, mpf, mpc, zeta
import mpmath
mp.dps = 20

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
    """L(s, chi) via Hurwitz-Zeta. Gueltig fuer Re(s) > 0, s != 1."""
    total = mpc(0)
    for a in range(1, q):
        ca = chi_fn(a)
        if ca == 0: continue
        total += ca * zeta(s, mpf(a)/q)
    return mpf(q)**(-s) * total

def find_zeros(chi_fn, q, t_max=20.0, fine_step=0.02):
    """Finde alle Minima-Kandidaten von |L(1/2+it,chi)| fuer t in (0, t_max]."""
    ts = np.arange(0.5, t_max, fine_step)
    vals = []
    for t in ts:
        try:
            v = L_at(mpf(1)/2 + 1j*mpf(t), chi_fn, q)
            vals.append(float(abs(v)))
        except Exception:
            vals.append(1e10)
    vals = np.array(vals)
    # Lokale Minima unter Schwelle
    candidates = []
    for i in range(1, len(ts)-1):
        if vals[i] < vals[i-1] and vals[i] < vals[i+1] and vals[i] < 0.4:
            candidates.append((float(ts[i]), vals[i]))
    # Pro Kandidat: feinere Verfeinerung
    zeros = []
    for t0, v0 in candidates:
        fine_ts = np.arange(t0 - fine_step, t0 + fine_step, fine_step/20)
        fvs = [float(abs(L_at(mpf(1)/2 + 1j*mpf(t), chi_fn, q))) for t in fine_ts]
        i_min = int(np.argmin(fvs))
        t_best = float(fine_ts[i_min])
        v_best = fvs[i_min]
        # Noch feiner
        fine2 = np.arange(t_best - fine_step/20, t_best + fine_step/20, fine_step/500)
        if len(fine2) > 0:
            fvs2 = [float(abs(L_at(mpf(1)/2 + 1j*mpf(t), chi_fn, q))) for t in fine2]
            i2 = int(np.argmin(fvs2))
            t_best = float(fine2[i2])
            v_best = fvs2[i2]
        if v_best < 0.05:   # echtes Minimum, bei der gewaehlten dps
            zeros.append(t_best)
    return zeros

# Die 10 Charaktere mit gemessenem Gap bei lam=20000, bestes N
chars_with_gaps = [
    ('chi_5',  5,  +0.0156),
    ('chi_8',  8,  -0.04902),
    ('chi_12', 12, +0.0337),
    ('chi_13', 13, -0.12504),
    ('chi_17', 17, +0.01318),
    ('chi_21', 21, -0.00423),
    ('chi_24', 24, +0.01076),
    ('chi_29', 29, +0.01439),
    ('chi_33', 33, -0.14221),
    ('chi_60', 60, +0.00521),
]

# ---------------------------------------------------------------
# Schritt 1: Null-Stellen berechnen
# ---------------------------------------------------------------
print("="*85)
print("Schritt 1: Null-Stellen bis t_max=20 via mpmath")
print("="*85)
zeros_by_char = {}
t0_total = time.time()
for name, D, _ in chars_with_gaps:
    ts_start = time.time()
    chi = make_chi_D(D)
    zeros = find_zeros(chi, D, t_max=20.0)
    dt = time.time() - ts_start
    zeros_by_char[name] = zeros
    print(f"  {name:>8}  D={D:>3}  n_zeros={len(zeros)}  "
          f"t_ks={[f'{z:.3f}' for z in zeros[:8]]}  [{dt:.1f}s]",
          flush=True)
print(f"\nSchritt 1 gesamt: {(time.time()-t0_total)/60:.1f} min\n")

# Speichere Null-Stellen
import json
with open('../_results/zeros_all_chars.json', 'w') as f:
    json.dump({k: v for k, v in zeros_by_char.items()}, f, indent=2)

# ---------------------------------------------------------------
# Schritt 2: Test-Funktionen evaluieren
# ---------------------------------------------------------------
print("="*85)
print("Schritt 2: Test-Funktionen-Summen vs. gemessener Gap")
print("="*85)

lam = 20000
L = math.log(lam)

def F_weil(t):
    """Real(1 / (1/2 + it)) = (1/2) / (1/4 + t^2)."""
    return 0.5 / (0.25 + t*t)

def F_cos_weil(t):
    """cos(L*t) * F_weil(t)."""
    return math.cos(L * t) * F_weil(t)

def F_gauss(t, sigma=1.0/L):
    """exp(-t^2 * sigma^2 / 2), normiert."""
    return math.exp(-(t*sigma)**2 / 2)

def F_fejer(t):
    """sin(Lt/2)^2 / (Lt/2)^2 -- Fejer-Kern."""
    x = L*t/2
    if abs(x) < 1e-9: return 1.0
    return (math.sin(x)/x)**2

def F_inv_t2(t):
    """1/t^2, dominance by low zeros."""
    return 1.0/(t*t) if t > 0.01 else 100.0

def F_sinc(t):
    """sin(Lt)/(Lt)."""
    x = L*t
    if abs(x) < 1e-9: return 1.0
    return math.sin(x)/x

tests = [
    ('F_weil = Re(-1/rho)', F_weil),
    ('F_cos_weil = cos(Lt)*F_weil', F_cos_weil),
    ('F_gauss(sigma=1/L)', F_gauss),
    ('F_fejer = sinc^2(Lt/2)', F_fejer),
    ('F_inv_t2 = 1/t^2', F_inv_t2),
    ('F_sinc = sin(Lt)/(Lt)', F_sinc),
]

print(f"\nlambda = {lam}, L = log(lambda) = {L:.4f}")
print()

# Tabelle: Fuer jede Test-Funktion die Summe S_chi = 2*Sum_k F(t_k) (Paar +t, -t)
results_table = {tname: {} for tname, _ in tests}
for name, D, gap in chars_with_gaps:
    tks = zeros_by_char[name]
    for tname, Ffn in tests:
        S = 2.0 * sum(Ffn(t) for t in tks)
        results_table[tname][name] = S

print(f"{'Char':>8}  {'gap':>10}  " + "  ".join(f"{tn[:22]:>22}" for tn, _ in tests))
for name, D, gap in chars_with_gaps:
    row = f"{name:>8}  {gap:>+10.5f}"
    for tname, _ in tests:
        S = results_table[tname][name]
        row += f"  {S:>+22.5f}"
    print(row)

# ---------------------------------------------------------------
# Korrelationen
# ---------------------------------------------------------------
print(f"\n{'='*85}")
print("Korrelationen und beste lineare Fits: gap ~ a * S + b")
print(f"{'='*85}")
gaps = np.array([g for _, _, g in chars_with_gaps])
for tname, _ in tests:
    Ss = np.array([results_table[tname][name] for name, _, _ in chars_with_gaps])
    # Korrelation
    r = np.corrcoef(Ss, gaps)[0, 1]
    # Fit
    a, b = np.polyfit(Ss, gaps, 1)
    gaps_pred = a * Ss + b
    r2 = 1 - np.sum((gaps - gaps_pred)**2) / np.sum((gaps - gaps.mean())**2)
    # Vorzeichen-Uebereinstimmung
    same_sign = sum(1 for s, g in zip(gaps_pred, gaps) if s*g > 0)
    print(f"  {tname:>28}:  r={r:+.3f}, R^2={r2:.3f}, slope={a:+.6f}, intercept={b:+.4f}, "
          f"sign_ok={same_sign}/10")

# ---------------------------------------------------------------
# Best-Kombinations-Regression (multivariate, alle Test-Funktionen)
# ---------------------------------------------------------------
print(f"\n{'='*85}")
print("Multivariate Regression auf allen Test-Funktionen")
print(f"{'='*85}")
X = np.array([[results_table[tn][name] for tn, _ in tests]
              for name, _, _ in chars_with_gaps])
y = gaps
X_aug = np.column_stack([X, np.ones(len(X))])
coef, _, _, _ = np.linalg.lstsq(X_aug, y, rcond=None)
y_pred = X_aug @ coef
r2 = 1 - np.sum((y - y_pred)**2) / np.sum((y - y.mean())**2)
same_sign = sum(1 for s, g in zip(y_pred, y) if s*g > 0)
print(f"R^2_multivariate = {r2:.4f}, sign_ok = {same_sign}/10")
for (tname, _), c in zip(tests, coef[:-1]):
    print(f"  {tname:>28}: {c:+.6f}")
print(f"  {'const':>28}: {coef[-1]:+.6f}")

# Speichere Ergebnisse
out = {
    'lam': lam, 'L': L,
    'chars': [{'name': n, 'D': D, 'gap': g, 'zeros': zeros_by_char[n]}
              for n, D, g in chars_with_gaps],
    'testfunctions': {tn: results_table[tn] for tn, _ in tests},
}
with open('../_results/explicit_formula_experiment.json', 'w') as f:
    json.dump(out, f, indent=2)
print(f"\n[saved] ../_results/explicit_formula_experiment.json")
