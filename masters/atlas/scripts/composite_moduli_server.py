#!/usr/bin/env python3
# coding: utf-8
"""
Server-Lauf: Test der Hypothese "zusammengesetzte teilerfremde Moduln zeigen Gap-Stabilitaet".

Test-Charaktere (primitive reelle Kronecker-Symbole chi_D mit D = fundamentale
Diskriminante > 0, chi_D(-1) = +1 = EVEN):

  ZUSAMMENGESETZT + TEILERFREMD (Hypothese: "strong"):
    D = 12 = 4·3    (chi_12, Referenz, bereits getestet)
    D = 21 = 3·7    neu
    D = 24 = 8·3    neu
    D = 33 = 3·11   neu
    D = 60 = 4·15 = 4·3·5  neu

  PRIMZAHL oder 2^k (Hypothese: "weak"):
    D = 5, 13       Primzahl (chi_5, chi_13, schon)
    D = 8           Primzahl-Potenz (chi_8a, schon)
    D = 17, 29      neu, Primzahl-Kontrollen

Geschaetzte Laufzeit: ~30 Min bei N=200, lambda bis 20000.

Output: composite_moduli_results.json
"""
import numpy as np
import sympy
from scipy.special import digamma
import json
import time
import math
import sys
sys.stdout.reconfigure(line_buffering=True)

# --- Kronecker-Symbol fuer fundamentale Diskriminanten ---
def kronecker_symbol(a, n):
    if n == 0:
        return 1 if abs(a) == 1 else 0
    if n < 0:
        return kronecker_symbol(a, -n) * (1 if a >= 0 else -1)
    result = 1
    while n % 2 == 0:
        if a % 2 == 0:
            return 0
        if a % 8 in (3, 5):
            result = -result
        n //= 2
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
    return result if n == 1 else 0

def chi_D(D):
    def chi(n):
        if math.gcd(n, abs(D)) != 1: return 0
        return kronecker_symbol(D, n)
    return chi

# --- Weil-Matrix (standard-Basis, kappa=0 fuer EVEN characters) ---
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

# ============================================================
# Test-Liste: (Label, D, Struktur)
# ============================================================
test_list = [
    ('D=17 (prim)',          17, 'primzahl'),
    ('D=21 = 3*7',           21, 'zusammengesetzt_teilerfremd'),
    ('D=24 = 8*3',           24, 'zusammengesetzt_teilerfremd'),
    ('D=29 (prim)',          29, 'primzahl'),
    ('D=33 = 3*11',          33, 'zusammengesetzt_teilerfremd'),
    ('D=60 = 4*3*5',         60, 'zusammengesetzt_teilerfremd'),
]

lam_vals = [100, 200, 500, 1000, 2000, 5000, 10000, 20000]
N = 200

results = []
t0 = time.time()
for label, D, struct in test_list:
    chi_fn = chi_D(D)
    # Quick Sanity: chi_D(-1) == +1?
    if chi_fn(-1) != 1:
        print(f"[skip] {label}: chi_D(-1) = {chi_fn(-1)} != +1 (ODD)", flush=True)
        continue
    print(f"\n=== {label} [{struct}] ===")
    for lam in lam_vals:
        L = np.log(lam)
        if N < 2*L*L*0.5:
            print(f"[skip] lam={lam} (N={N} < 0.5*2*L^2={int(L*L)})", flush=True)
            continue
        ts = time.time()
        try:
            ec, es, g = gap_at(chi_fn, lam, N, kappa=0)
            dt = time.time()-ts
            row = {'label': label, 'D': D, 'struct': struct, 'lambda': lam,
                   'N': N, 'L': L, 'lam1_cos': ec, 'lam1_sin': es, 'gap': g,
                   'gap_sqrtL': g/np.sqrt(L), 't': dt}
            results.append(row)
            print(f"  lam={lam:6d}  gap={g:+8.4f}  gap/sqrtL={g/np.sqrt(L):+7.4f}  t={dt:5.1f}s",
                  flush=True)
            with open('composite_moduli_results.json', 'w') as f:
                json.dump(results, f, indent=2)
        except Exception as e:
            print(f"  lam={lam}: ERROR {e}", flush=True)

print(f"\n[done] {(time.time()-t0)/60:.1f} min, n={len(results)}")

# Auswertung
print(f"\n{'='*70}")
print("=== KONSOLIDIERUNG ===")
print(f"{'='*70}")
print(f"{'Charakter':>16}  {'struct':>26}  {'pos/n':>7}  {'mean':>8}  "
      f"{'std':>7}  {'slope':>10}  {'slope_err':>10}")
for label, D, struct in test_list:
    sub = [r for r in results if r['label'] == label]
    if not sub: continue
    gaps = np.array([r['gap'] for r in sub])
    lams = np.array([r['lambda'] for r in sub])
    pos = (gaps > 0).sum()
    xs = np.log(lams)
    try:
        A, cov = np.polyfit(xs, gaps, 1, cov=True)
        err = np.sqrt(cov[0,0])
    except:
        A = [0, 0]; err = 0
    flag = ' STRONG' if (pos == len(gaps) and abs(A[0]) < 2*err) else ''
    print(f"{label:>16}  {struct:>26}  {pos:3d}/{len(gaps):3d}  "
          f"{gaps.mean():+8.4f}  {gaps.std():7.4f}  "
          f"{A[0]:+10.4f}  {err:10.4f}{flag}")

print("\n[end]")
