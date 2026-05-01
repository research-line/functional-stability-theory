#!/usr/bin/env python3
# coding: utf-8
"""
Scan-Test: Gap-Stabilitaet bei variierendem lambda und N.
Pruft, ob das Flipping bei lambda=1000 ein Truncation-Artefakt ist.
"""
import numpy as np
import sympy
from scipy.special import digamma

def chi4(n):
    r = n % 4
    if r == 1: return 1
    if r == 3: return -1
    return 0

def cos_cos_int(n, m, a, b, L, delta):
    def F(k, phi, t):
        if k == 0: return np.cos(phi) * t
        return (L / (k * np.pi)) * np.sin(k * np.pi * t / L + phi)
    r = 0.0
    r += 0.5 * (F(n-m, m*np.pi*delta/L, b) - F(n-m, m*np.pi*delta/L, a))
    r += 0.5 * (F(n+m, -m*np.pi*delta/L, b) - F(n+m, -m*np.pi*delta/L, a))
    return r

def sin_sin_int(n, m, a, b, L, delta):
    def F(k, phi, t):
        if k == 0: return np.cos(phi) * t
        return (L / (k * np.pi)) * np.sin(k * np.pi * t / L + phi)
    r = 0.0
    r += 0.5 * (F(n-m, m*np.pi*delta/L, b) - F(n-m, m*np.pi*delta/L, a))
    r -= 0.5 * (F(n+m, -m*np.pi*delta/L, b) - F(n+m, -m*np.pi*delta/L, a))
    return r

def cos_ov(n, m, delta, L):
    if abs(delta) >= 2*L: return 0.0
    a = max(-L, -L+delta); b = min(L, L+delta)
    if a >= b: return 0.0
    I = cos_cos_int(n, m, a, b, L, delta)
    ns = 1/np.sqrt(L) if n>0 else 1/np.sqrt(2*L)
    ms = 1/np.sqrt(L) if m>0 else 1/np.sqrt(2*L)
    return I * ns * ms

def sin_ov(n, m, delta, L):
    if abs(delta) >= 2*L: return 0.0
    a = max(-L, -L+delta); b = min(L, L+delta)
    if a >= b: return 0.0
    return sin_sin_int(n, m, a, b, L, delta) / L

def h_arch(tau, kappa):
    shift = (2*kappa+1)/4.0
    return digamma(shift + 1j*tau/2).real - np.log(np.pi)

def build_W(sector, N, L, primes, chi_vals, kappa):
    W = np.zeros((N, N))
    if sector == 'cos':
        bs = 0; ov = cos_ov
    else:
        bs = 1; ov = sin_ov
    for i in range(N):
        n = i + bs
        W[i, i] += h_arch(np.pi*n/L, kappa)
    for p, cp in zip(primes, chi_vals):
        if cp == 0: continue
        lp = np.log(p)
        mm = int(2*L/lp)
        for me in range(1, mm+1):
            d = me*lp
            if d >= 2*L: break
            w = (cp**me) * lp / (p**(me/2.0))
            for i in range(N):
                for j in range(N):
                    n = i+bs; m = j+bs
                    W[i,j] += w * (ov(n,m,d,L) + ov(n,m,-d,L))
    return 0.5 * (W + W.T)

def gap_at(lam, N, kappa=1):
    L = np.log(lam)
    primes = list(sympy.primerange(3, int(lam)+1))
    chi_vals = [chi4(p) for p in primes]
    W_c = build_W('cos', N, L, primes, chi_vals, kappa)
    W_s = build_W('sin', N, L, primes, chi_vals, kappa)
    ec = np.linalg.eigvalsh(W_c); es = np.linalg.eigvalsh(W_s)
    return ec[0], es[0], es[0] - ec[0]

print("=== Truncation-Stabilitaet: Gap als Funktion von N bei lambda = 1000 ===")
print(f"{'N':>4}  {'lam1+':>10}  {'lam1-':>10}  {'gap':>10}  {'dominance':>12}")
for N in [10, 12, 14, 16, 18, 20, 24, 28, 32, 40]:
    ec, es, g = gap_at(1000, N, kappa=1)
    dom = 'ODD' if g < 0 else 'EVEN'
    print(f"{N:4d}  {ec:+10.4f}  {es:+10.4f}  {g:+10.4f}  {dom:>12}")

print("\n=== Skalierung: Gap ueber lambda bei hinreichend grossem N ===")
print(f"{'lambda':>7}  {'N':>4}  {'lam1+':>10}  {'lam1-':>10}  {'gap':>10}  {'gap/sqrt(L)':>14}")
for lam in [30, 60, 100, 200, 500, 1000, 2000, 5000]:
    # N skaliert mit L^2 fuer Stabilitaet
    L = np.log(lam)
    N = max(16, int(2.5 * L * L))
    N = min(N, 50)  # cap
    ec, es, g = gap_at(lam, N, kappa=1)
    dom = 'ODD' if g < 0 else 'EVEN'
    print(f"{lam:7d}  {N:4d}  {ec:+10.4f}  {es:+10.4f}  {g:+10.4f}  {g/np.sqrt(L):+14.4f}  {dom}")
