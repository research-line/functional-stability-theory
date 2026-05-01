#!/usr/bin/env python3
# coding: utf-8
"""
Wiederholung der Multivariaten Regression auf N=400-Daten.

Input:
  - all10_high_N_results.json (N=400 fuer chi_8, chi_13, chi_17, chi_21,
    chi_24, chi_29, chi_33, chi_60 -- 8 chars, 4 lambdas)
  - chi12_high_N_results.json (chi_5, chi_12 bei N in {200,300,400})
  - char_features_all.json (t_1, L(1), -L'/L, log D aus Session 4)

Output:
  - N=400-Mean je Charakter
  - univariate und multivariate Regression auf N=400-Daten
"""
import json
import math
import numpy as np
from pathlib import Path

R = Path(__file__).parent.parent / "_results"

# Lade die drei JSON-Dateien
with open(R / "all10_high_N_results.json") as f:
    all10 = json.load(f)
with open(R / "chi12_high_N_results.json") as f:
    chi12_data = json.load(f)
with open(R / "char_features_all.json") as f:
    features = json.load(f)

# Baue Mean-Gap bei N=400 je Charakter (nur lambdas 500, 2000, 10000, 20000)
lambda_set = {500, 2000, 10000, 20000}

def mean_N400(chi_name, records):
    vals = [r['gap'] for r in records
            if r['chi'] == chi_name and r['N'] == 400
            and r['lambda'] in lambda_set]
    return vals

gaps_N400 = {}
for name in ['chi_5', 'chi_12']:
    vals = mean_N400(name, chi12_data)
    if vals:
        gaps_N400[name] = vals
for name in ['chi_8', 'chi_13', 'chi_17', 'chi_21', 'chi_24', 'chi_29', 'chi_33', 'chi_60']:
    vals = mean_N400(name, all10)
    if vals:
        gaps_N400[name] = vals

print("=" * 85)
print("N=400 MEAN-GAP JE CHARAKTER (Lambda in {500, 2000, 10000, 20000})")
print("=" * 85)
print(f"{'Char':>8}  {'D':>3}  {'t_1':>7}  {'N=200 mean':>11}  {'N=400 mean':>11}  "
      f"{'std':>7}  {'pos/4':>5}  {'Skala Shrink':>12}")

# N=200 mean aus MULTIVARIATE_REGRESSION_2026-04-16.md
mean_N200 = {
    'chi_5':  0.160, 'chi_8':  0.055, 'chi_12': 0.131, 'chi_13': 0.271,
    'chi_17': 0.026, 'chi_21': 0.094, 'chi_24': 0.001, 'chi_29': 0.132,
    'chi_33': -0.011, 'chi_60': 0.157,
}

# Erstelle Feature-Dict nach Name
feat_by_name = {f['name']: f for f in features}

rows = []
for name in ['chi_5', 'chi_8', 'chi_12', 'chi_13', 'chi_17',
             'chi_21', 'chi_24', 'chi_29', 'chi_33', 'chi_60']:
    vals = gaps_N400.get(name, [])
    if not vals:
        print(f"{name:>8}  ---  no N=400 data")
        continue
    m = sum(vals)/len(vals)
    s = math.sqrt(sum((v-m)**2 for v in vals) / max(1, len(vals)-1)) if len(vals) > 1 else 0
    pos = sum(1 for v in vals if v > 0)
    f = feat_by_name.get(name, {})
    t1 = f.get('t1', float('nan'))
    L1 = f.get('L1', float('nan'))
    ratio = f.get('logL_prime_over_L', float('nan'))
    logD = f.get('log_D', float('nan'))
    D = f.get('D', '?')
    m200 = mean_N200.get(name, float('nan'))
    shrink = abs(m200/m) if m != 0 else float('inf')
    print(f"{name:>8}  {D:>3}  {t1:>7.3f}  {m200:>+11.4f}  {m:>+11.5f}  "
          f"{s:>7.4f}  {pos}/{len(vals):<2}  {shrink:>12.2f}x")
    rows.append({'name': name, 'D': D, 't1': t1, 'L1': L1,
                 'logL_prime_over_L': ratio, 'log_D': logD,
                 'gap_N400_mean': m, 'gap_N400_std': s, 'pos': pos, 'n': len(vals),
                 'gap_N200_mean': m200})

# ====================================================================
# REGRESSION AUF N=400-MEANS
# ====================================================================
print()
print("=" * 85)
print("MULTIVARIATE REGRESSION auf N=400 MEAN-GAPS")
print("=" * 85)

valid = [r for r in rows if not math.isnan(r['t1'])]
print(f"\nN = {len(valid)} Charaktere")

y = np.array([r['gap_N400_mean'] for r in valid])
X = np.array([[r['t1'], r['L1'], r['logL_prime_over_L'], r['log_D']] for r in valid])

features_names = ['t_1', 'L(1)', '-Lprime/L', 'log D']

# Univariate
print(f"\n{'Feature':>20}  {'Slope':>10}  {'Intercept':>10}  {'R^2':>8}")
for i, feat in enumerate(features_names):
    x = X[:, i]
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
    print(f"{feat:>20}  {slope:>+10.5f}  {intercept:>+10.5f}  {r2:>8.4f}")

# Multivariate
print(f"\nMultivariate Regression (alle 4 Features + intercept):")
X_aug = np.column_stack([X, np.ones(len(X))])
coef, _, _, _ = np.linalg.lstsq(X_aug, y, rcond=None)
y_pred_full = X_aug @ coef
r2_full = 1 - np.sum((y - y_pred_full)**2) / np.sum((y - y.mean())**2)
for feat, c in zip(features_names + ['const'], coef):
    print(f"  {feat:>12}: {c:+.5f}")
print(f"  R^2 = {r2_full:.4f}")

# Vergleich N=200 vs N=400 Regression
print(f"\n{'='*85}")
print("VERGLEICH: wie stark aendert sich die Regression von N=200 zu N=400?")
print("="*85)
y200 = np.array([r['gap_N200_mean'] for r in valid])
print(f"\nBeste N=200 Regression (aus MULTIVARIATE_REGRESSION_2026-04-16.md):")
print(f"  L(1): slope -0.180, R^2 = 0.31")

print(f"\nN=400 Regression:")
for i, feat in enumerate(features_names):
    x = X[:, i]
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    r2 = 1 - np.sum((y - y_pred)**2) / np.sum((y - y.mean())**2)
    print(f"  {feat}: slope {slope:+.4f}, R^2 = {r2:.4f}")

# Korrelation gap_N200 vs gap_N400
from numpy import corrcoef
r = corrcoef(y200, y)[0, 1]
print(f"\nKorrelation(gap_N200, gap_N400) = {r:.3f}")
print(f"  d.h. N=200-Werte sagen nur {r**2*100:.0f}% der N=400-Variabilitaet vorher.")
print(f"  Die Skala 'schrumpft' um Faktor 2-16 je nach Charakter.")

# Konsistenz der Vorzeichen
same_sign = sum(1 for a, b in zip(y200, y) if a*b > 0)
print(f"\nGleiches Vorzeichen (N=200 vs N=400): {same_sign}/{len(y)}")

# Speichern
out = {
    'description': 'N=400 Regression (Session 5, 2026-04-16 Teil 2)',
    'rows': rows,
    'coef_multivar': list(float(c) for c in coef),
    'r2_multivar': float(r2_full),
    'corr_N200_N400': float(r),
}
with open(R / "regression_N400_results.json", 'w') as f:
    json.dump(out, f, indent=2)

print(f"\n[saved] {R / 'regression_N400_results.json'}")
