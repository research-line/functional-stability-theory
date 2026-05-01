#!/usr/bin/env python3
# coding: utf-8
"""
Server-Lauf: chi_12 bei hohem N (300, 400) - ist die Gap-Konstanz echt?

Hintergrund: Bei N=200 zeigt chi_12 in allen 9 Tests (lambda bis 50000)
positive Gaps mit slope ~ 0. Frage: ist dies ein N-Konvergenz-Artefakt?

Test:
  chi_12 bei N in {200, 300, 400}
  lambda in {500, 2000, 10000, 20000, 50000}

  Plus Kontrolle: chi_5 bei gleichen N und lambda, um N-Sensitivitaet zu
  charakterisieren.

Laufzeit-Schaetzung (bei N=400, lambda=50000): ~30 Min pro Eval (Eigenwert
O(N^3), N^3/200^3 = 8 fach).
Gesamt: 2 chars * 3 N * 5 lambda = 30 evals, geschaetzt 3 h auf CCX13.

Wenn zu lang: setze max lambda = 20000 (Laufzeit ~ 1h).
"""
import numpy as np
import sympy
from scipy.special import digamma
import json
import time
import math
import sys
sys.stdout.reconfigure(line_buffering=True)

def chi12(n):
    if math.gcd(n, 12) != 1: return 0
    return {1: 1, 5: -1, 7: -1, 11: 1}[n % 12]

def chi5(n):
    if n % 5 == 0: return 0
    return {1: 1, 2: -1, 3: -1, 4: 1}[n % 5]

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

# Conservative: lambda bis 20000, N bis 400
chars = [('chi_12', chi12), ('chi_5', chi5)]
N_vals = [200, 300, 400]
lam_vals = [500, 2000, 10000, 20000]

results = []
t0 = time.time()
for name, fn in chars:
    for N in N_vals:
        for lam in lam_vals:
            L = np.log(lam)
            if N < L*L:
                print(f"[skip] {name} N={N} lam={lam} (L^2={L*L:.0f})", flush=True)
                continue
            ts = time.time()
            try:
                ec, es, g = gap_at(fn, lam, N, kappa=0)
                dt = time.time() - ts
                row = {'chi': name, 'N': N, 'lambda': lam, 'L': L,
                       'lam1_cos': ec, 'lam1_sin': es, 'gap': g,
                       'gap_sqrtL': g/np.sqrt(L), 't': dt}
                results.append(row)
                print(f"[{name}] N={N:3d} lam={lam:5d}  gap={g:+8.5f}  "
                      f"gap/sqrtL={g/np.sqrt(L):+7.5f}  t={dt:6.1f}s",
                      flush=True)
                with open('chi12_high_N_results.json', 'w') as f:
                    json.dump(results, f, indent=2)
            except Exception as e:
                print(f"[error] {name} N={N} lam={lam}: {e}", flush=True)

print(f"\n[done] {(time.time()-t0)/60:.1f} min, n={len(results)}")

# Auswertung: N-Konvergenz
print(f"\n{'='*70}")
print("=== N-KONVERGENZ ===")
print(f"{'='*70}")
for chi in ['chi_12', 'chi_5']:
    print(f"\n--- {chi} ---")
    print(f"{'lambda':>7}  " + "  ".join(f"N={n}" for n in N_vals))
    for lam in lam_vals:
        row = [f"{lam:7d}"]
        for N in N_vals:
            sub = [r for r in results if r['chi']==chi and r['N']==N and r['lambda']==lam]
            if sub:
                row.append(f"{sub[0]['gap']:+8.5f}")
            else:
                row.append("  (skip)")
        print("  ".join(row))
