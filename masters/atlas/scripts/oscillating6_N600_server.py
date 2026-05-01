#!/usr/bin/env python3
# coding: utf-8
"""
Task #37: N=600-Lauf fuer die 6 oszillierenden Charaktere aus der N=400-
Analyse. Ziel: klaeren, ob die Oszillationen bei N=400 Truncation-Resten
sind (dann bei N=600 kleiner) oder echt (dann bleiben sie).

Entscheidend: chi_33 bei N=200 und N=400 beide stabil negativ. Bleibt das
bei N=600?

Kosten: N=600 ist (600/400)^3 = 3.375x teurer als N=400 fuer O(N^3)-
Eigenwerte. Pro Char ca. 34 Min, 6 Chars ca. 3-4 h.
"""
import numpy as np
import sympy
from scipy.special import digamma
import json
import time
import math
import sys
sys.stdout.reconfigure(line_buffering=True)

# -- Kronecker-Symbol --
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

# -- Die 6 oszillierenden Charaktere bei N=400 --
chars = [
    ('chi_8',  make_chi_D(8)),
    ('chi_13', make_chi_D(13)),
    ('chi_21', make_chi_D(21)),
    ('chi_29', make_chi_D(29)),
    ('chi_33', make_chi_D(33)),
    ('chi_60', make_chi_D(60)),
]

N_val = 600
lam_vals = [500, 2000, 10000, 20000]

results = []
t0 = time.time()
for name, fn in chars:
    for lam in lam_vals:
        L = np.log(lam)
        if N_val < L*L:
            print(f"[skip] {name} N={N_val} lam={lam} (L^2={L*L:.0f})", flush=True)
            continue
        ts = time.time()
        try:
            ec, es, g = gap_at(fn, lam, N_val, kappa=0)
            dt = time.time() - ts
            row = {'chi': name, 'N': N_val, 'lambda': lam, 'L': L,
                   'lam1_cos': ec, 'lam1_sin': es, 'gap': g,
                   'gap_sqrtL': g/np.sqrt(L), 't': dt}
            results.append(row)
            print(f"[{name}] N={N_val:3d} lam={lam:5d}  gap={g:+8.5f}  "
                  f"gap/sqrtL={g/np.sqrt(L):+7.5f}  t={dt:7.1f}s",
                  flush=True)
            with open('oscillating6_N600_results.json', 'w') as f:
                json.dump(results, f, indent=2)
        except Exception as e:
            print(f"[error] {name} N={N_val} lam={lam}: {e}", flush=True)

print(f"\n[done] {(time.time()-t0)/60:.1f} min, n={len(results)}", flush=True)

# Mean-Gap je Char plus Vergleich mit N=400
print(f"\n{'='*95}")
print("=== MEAN-GAP BEI N=600 JE CHARAKTER + VERGLEICH N=400 ===")
print(f"{'='*95}")
# N=400-Means aus REVISION_N400_2026-04-16.md
mean_N400 = {'chi_8': 0.00744, 'chi_13': 0.07707, 'chi_21': 0.07302,
             'chi_29': 0.01053, 'chi_33': -0.03361, 'chi_60': 0.03675}
print(f"{'Char':>8}  {'lam=500':>10}  {'lam=2000':>10}  {'lam=10000':>10}  {'lam=20000':>10}  "
      f"{'mean600':>10}  {'mean400':>10}  {'pos':>4}")
for name, _ in chars:
    subs = [r for r in results if r['chi']==name]
    if not subs: continue
    vals = [r['gap'] for r in subs]
    pos = sum(1 for v in vals if v > 0)
    m = sum(vals)/len(vals)
    row = f"{name:>8}"
    for lam in lam_vals:
        sub = [r for r in subs if r['lambda']==lam]
        row += f"  {sub[0]['gap']:+10.5f}" if sub else "     (skip)"
    row += f"  {m:+10.5f}  {mean_N400.get(name,0):+10.5f}  {pos}/{len(vals)}"
    print(row)
