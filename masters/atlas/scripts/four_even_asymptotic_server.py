#!/usr/bin/env python3
# coding: utf-8
"""
Server-Lauf: Asymptotik-Test aller 4 even characters bei N=200.

Ziel: Klassifikation "strong even" (chi_12 gezeigt) vs. "weak even" (chi_5)
empirisch festigen. Test-Grid:
  chi in {chi_5, chi_8a, chi_12, chi_13}
  N in {160, 200}
  lambda in {100, 200, 500, 1000, 2000, 5000, 10000, 20000}

Output: four_even_asymptotic_results.json
"""
import numpy as np
import sympy
from scipy.special import digamma
import json
import time
import sys
sys.stdout.reconfigure(line_buffering=True)

def chi5(n):
    if n % 5 == 0: return 0
    return {1: 1, 2: -1, 3: -1, 4: 1}[n % 5]

def chi8a(n):
    if n % 2 == 0: return 0
    return {1: 1, 3: -1, 5: -1, 7: 1}[n % 8]

def chi12(n):
    if np.gcd(n, 12) != 1: return 0
    return {1: 1, 5: -1, 7: -1, 11: 1}[n % 12]

def chi13(n):
    if n % 13 == 0: return 0
    qr = {pow(i, 2, 13) for i in range(1, 13)}
    return 1 if (n % 13) in qr else -1

def build_W(sector, N, L, primes, chi_vals, kappa):
    bs = 0 if sector == 'cos' else 1
    idx = np.arange(N) + bs
    nn = idx[:, None]; mm = idx[None, :]
    W = np.zeros((N, N))
    shift = (2*kappa+1)/4.0
    tau = np.pi * idx / L
    diag = np.array([digamma(shift + 1j*t/2).real - np.log(np.pi) for t in tau])
    np.fill_diagonal(W, diag)
    if sector == 'cos':
        norm = np.where(idx > 0, 1.0/np.sqrt(L), 1.0/np.sqrt(2*L))
    else:
        norm = np.full(N, 1.0/np.sqrt(L))
    NM = np.outer(norm, norm)

    def A(k, phi, t, Lv):
        r = np.zeros_like(k, dtype=float)
        mk = (k != 0)
        r[mk] = (Lv / (k[mk] * np.pi)) * np.sin(k[mk] * np.pi * t / Lv + phi[mk])
        r[~mk] = np.cos(phi[~mk]) * t
        return r

    for p, cp in zip(primes, chi_vals):
        if cp == 0: continue
        lp = np.log(p); mmx = int(2*L/lp)
        for me in range(1, mmx+1):
            d = me*lp
            if d >= 2*L: break
            w = (cp**me) * lp / (p**(me/2.0))
            for ds in [d, -d]:
                a = max(-L, -L+ds); b = min(L, L+ds)
                if a >= b: continue
                k1 = nn-mm; k2 = nn+mm
                ph1 = mm*np.pi*ds/L; ph2 = -mm*np.pi*ds/L
                P1 = np.broadcast_to(ph1, (N,N)).copy()
                P2 = np.broadcast_to(ph2, (N,N)).copy()
                I1 = 0.5*(A(k1, P1, b, L) - A(k1, P1, a, L))
                I2 = 0.5*(A(k2, P2, b, L) - A(k2, P2, a, L))
                ovl = I1+I2 if sector == 'cos' else I1-I2
                W += w * ovl * NM
    return 0.5*(W + W.T)

def gap_at(chi_fn, lam, N, kappa):
    L = np.log(lam)
    primes = [p for p in sympy.primerange(2, int(lam)+1) if chi_fn(p) != 0]
    cv = [chi_fn(p) for p in primes]
    W_c = build_W('cos', N, L, primes, cv, kappa)
    W_s = build_W('sin', N, L, primes, cv, kappa)
    ec = np.linalg.eigvalsh(W_c); es = np.linalg.eigvalsh(W_s)
    return float(ec[0]), float(es[0]), float(es[0]-ec[0])

chars = [('chi_5', chi5), ('chi_8a', chi8a), ('chi_12', chi12), ('chi_13', chi13)]
N_vals = [160, 200]
lam_vals = [100, 200, 500, 1000, 2000, 5000, 10000, 20000]

results = []
t0 = time.time()
for name, fn in chars:
    for N in N_vals:
        for lam in lam_vals:
            L = np.log(lam)
            if N < 0.5*2*L*L:
                print(f"[skip] {name} lam={lam} N={N}", flush=True)
                continue
            ts = time.time()
            ec, es, g = gap_at(fn, lam, N, kappa=0)
            dt = time.time()-ts
            row = {'chi': name, 'lambda': lam, 'N': N, 'L': L,
                   'lam1_cos': ec, 'lam1_sin': es, 'gap': g,
                   'gap_sqrtL': g/np.sqrt(L), 't': dt}
            results.append(row)
            print(f"[{name}] N={N:3d} lam={lam:6d} gap={g:+8.4f} "
                  f"gap/sqrtL={g/np.sqrt(L):+7.4f} t={dt:5.1f}s", flush=True)
            with open('four_even_results.json', 'w') as f:
                json.dump(results, f, indent=2)

print(f"\n[done] {(time.time()-t0)/60:.1f} min, n={len(results)}")

# Auswertung
print(f"\n{'='*70}")
print("=== KONSOLIDIERUNG (N=200) ===")
print(f"{'='*70}")
print(f"{'chi':>8}  {'pos':>5}/{'n':>3}  {'mean':>8}  {'std':>7}  "
      f"{'min':>8}  {'max':>8}  {'slope':>10}  {'slope_err':>10}")
for name, _ in chars:
    sub = [r for r in results if r['chi'] == name and r['N'] == 200]
    if not sub: continue
    gaps = np.array([r['gap'] for r in sub])
    lams = np.array([r['lambda'] for r in sub])
    pos = (gaps > 0).sum()
    # Fit gap = a log(lam) + b
    xs = np.log(lams)
    A, cov = np.polyfit(xs, gaps, 1, cov=True)
    err = np.sqrt(cov[0,0])
    print(f"{name:>8}  {pos:5d}/{len(sub):3d}  {gaps.mean():+8.4f}  "
          f"{gaps.std():7.4f}  {gaps.min():+8.4f}  {gaps.max():+8.4f}  "
          f"{A[0]:+10.4f}  {err:10.4f}")

print("\nKlassifikation:")
print("  STRONG EVEN: alle gap > 0, slope konsistent mit 0.")
print("  WEAK EVEN: oszillatorisch oder abklingend.")
