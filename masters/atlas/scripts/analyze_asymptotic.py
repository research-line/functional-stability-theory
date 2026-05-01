#!/usr/bin/env python3
# coding: utf-8
"""Auswertung asymptotic_results.json."""
import json
import numpy as np

with open("../_results/asymptotic_results.json") as f:
    results = json.load(f)

print(f"Total rows: {len(results)}\n")

for chi in ['chi_5', 'chi_12']:
    print(f"\n{'='*70}")
    print(f"=== {chi} ===")
    print(f"{'='*70}")
    for N in [80, 120, 160, 200]:
        rows = [r for r in results if r['chi'] == chi and r['N'] == N]
        if not rows: continue
        gaps = np.array([r['gap'] for r in rows])
        lams = np.array([r['lambda'] for r in rows])
        pos = (gaps > 0).sum()
        print(f"\n-- N = {N} (n={len(rows)}) --")
        print(f"  Positive gaps: {pos}/{len(rows)}   mean={gaps.mean():+.4f}   "
              f"median={np.median(gaps):+.4f}   std={gaps.std():.4f}")
        # All gaps
        for r in rows:
            sign = '+' if r['gap'] > 0 else '-'
            print(f"    lam={r['lambda']:6d}  gap={r['gap']:+.4f}  "
                  f"|gap|/sqrt(L)={abs(r['gap'])/np.sqrt(r['L']):.4f}  "
                  f"|gap|/L={abs(r['gap'])/r['L']:.4f}")

# Regression: |gap| ~ lambda^alpha, mit allen (nicht nur positiven) Werten
print(f"\n{'='*70}")
print("=== Log-Log Regression |gap| vs lambda (beide Charaktere, N=200) ===")
print(f"{'='*70}")
for chi in ['chi_5', 'chi_12']:
    rows = [r for r in results if r['chi'] == chi and r['N'] == 200]
    if not rows: continue
    xs = np.log([r['lambda'] for r in rows])
    ys = np.log([abs(r['gap']) for r in rows])
    A = np.polyfit(xs, ys, 1)
    print(f"  {chi}: log|gap| = {A[0]:+.4f} * log(lam) + {A[1]:+.4f}")
    print(f"    → |gap| ~ lambda^({A[0]:+.3f}) * exp({A[1]:+.3f}) = "
          f"{np.exp(A[1]):.4f} * lambda^{A[0]:+.3f}")

# Fuer chi_12 separat signed gap (alle positiv)
print("\n=== chi_12 N=200: Signed Gap Konvergenz ===")
rows = [r for r in results if r['chi'] == 'chi_12' and r['N'] == 200]
if rows:
    gaps = [r['gap'] for r in rows]
    lams = [r['lambda'] for r in rows]
    print(f"  Alle Werte positiv: {all(g > 0 for g in gaps)}")
    print(f"  Min gap: {min(gaps):.4f}   Max gap: {max(gaps):.4f}")
    print(f"  Mean gap: {np.mean(gaps):.4f}   Std: {np.std(gaps):.4f}")
    # Test: konstant?
    xs = np.log(lams)
    A, cov = np.polyfit(xs, gaps, 1, cov=True)
    err_slope = np.sqrt(cov[0,0])
    print(f"  Fit gap = a*log(lam) + b:  a = {A[0]:+.4f} ± {err_slope:.4f}")
    print(f"  -> slope konsistent mit 0 (konstanter Gap): "
          f"{'JA' if abs(A[0]) < 2*err_slope else 'NEIN'}")
