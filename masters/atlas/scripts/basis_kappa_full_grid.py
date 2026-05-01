#!/usr/bin/env python3
# coding: utf-8
"""
Voller 2x2x2-Grid: Charakter x Basis x Kappa.

Motivation: twisted_basis_chi4.py zeigte unerwartete Paritaets-Flips
  chi_5 + twisted + kappa=0  zeigte ODD-Dom (wo EVEN erwartet war).

Hypothese zu testen: **Basis und Kappa muessen zusammen geflippt werden**,
  d.h. nicht nur die Basis, sondern auch der Gamma-Faktor muss angepasst werden.

Vier Dimensionen:
  chi in {chi_4 odd, chi_5 even}
  basis in {standard, twisted}
  kappa in {0, 1}

Naturliche Erwartung (Meta-Theorie v0.5):
  'Natural' Matching ist: Paritaet(chi) = Paritaet(basis) = Paritaet(kappa).
  - chi_5 even + standard + kappa=0:  alle EVEN  -> Match (EVEN-Dom)       [bestaetigt]
  - chi_4 odd  + twisted  + kappa=1:  alle ODD   -> Match (EVEN-Dom?)      [zu testen]

Das 'sieht den Gamma-Faktor' aus GAMMA_BASIS_DUALITAET wird hier konkret
geprueft: welche Kombination (Basis, Kappa) ist die natuerliche fuer chi.
"""
import numpy as np
import sympy
from scipy.special import digamma
import sys
sys.stdout.reconfigure(line_buffering=True)

def chi4(n):
    if n % 2 == 0: return 0
    return 1 if (n % 4) == 1 else -1

def chi5(n):
    if n % 5 == 0: return 0
    return {1: 1, 2: -1, 3: -1, 4: 1}[n % 5]

def build_W_matrix(sector, basis, N, L, primes, chi_vals, kappa):
    if basis == 'standard':
        bs = 0 if sector == 'cos' else 1
        idx = np.arange(N) + bs
        omega = np.pi * idx / L
        wav = idx.astype(float)
    else:
        idx = np.arange(N)
        wav = idx.astype(float) + 0.5
        omega = np.pi * wav / L

    W = np.zeros((N, N))
    shift = (2*kappa+1) / 4.0
    diag = np.array([digamma(shift + 1j*t/2).real - np.log(np.pi) for t in omega])
    np.fill_diagonal(W, diag)

    if basis == 'standard' and sector == 'cos':
        norm = np.where(idx > 0, 1.0/np.sqrt(L), 1.0/np.sqrt(2*L))
    else:
        norm = np.full(N, 1.0/np.sqrt(L))
    norm_matrix = np.outer(norm, norm)

    def A(k_arr, phi_arr, t, Lval):
        result = np.zeros_like(k_arr, dtype=float)
        mask = (k_arr != 0)
        kk = k_arr[mask].astype(float)
        pp = phi_arr[mask]
        result[mask] = (Lval / (kk * np.pi)) * np.sin(kk * np.pi * t / Lval + pp)
        result[~mask] = np.cos(phi_arr[~mask]) * t
        return result

    nn = wav[:, None]
    mm = wav[None, :]
    k1f = nn - mm
    k2f = nn + mm
    k1 = np.rint(k1f).astype(int)
    k2 = np.rint(k2f).astype(int)
    assert np.allclose(k1, k1f, atol=1e-10) and np.allclose(k2, k2f, atol=1e-10)

    for p, cp in zip(primes, chi_vals):
        if cp == 0: continue
        lp = np.log(p)
        mm_max = int(2*L / lp)
        for me in range(1, mm_max+1):
            d = me * lp
            if d >= 2*L: break
            weight = (cp**me) * lp / (p**(me/2.0))
            for d_signed in [d, -d]:
                a = max(-L, -L + d_signed)
                b = min(L, L + d_signed)
                if a >= b: continue
                phi1 = mm * np.pi * d_signed / L
                phi2 = -mm * np.pi * d_signed / L
                phi1_full = np.broadcast_to(phi1, (N, N)).copy()
                phi2_full = np.broadcast_to(phi2, (N, N)).copy()
                I1 = 0.5 * (A(k1, phi1_full, b, L) - A(k1, phi1_full, a, L))
                I2 = 0.5 * (A(k2, phi2_full, b, L) - A(k2, phi2_full, a, L))
                if sector == 'cos':
                    ovl = I1 + I2
                else:
                    ovl = I1 - I2
                W += weight * ovl * norm_matrix
    return 0.5 * (W + W.T)

def gap_at(chi_fn, lam, N, kappa, basis):
    L = np.log(lam)
    primes = [p for p in sympy.primerange(2, int(lam)+1) if chi_fn(p) != 0]
    chi_vals = [chi_fn(p) for p in primes]
    W_c = build_W_matrix('cos', basis, N, L, primes, chi_vals, kappa)
    W_s = build_W_matrix('sin', basis, N, L, primes, chi_vals, kappa)
    ec = np.linalg.eigvalsh(W_c); es = np.linalg.eigvalsh(W_s)
    return float(ec[0]), float(es[0]), float(es[0] - ec[0])

lam_values = [30, 50, 100, 200, 500, 1000, 2000, 5000]

chars = [('chi_4', chi4, 'ODD'), ('chi_5', chi5, 'EVEN')]
bases = ['standard', 'twisted']
kappas = [0, 1]

summary = []
for chi_name, chi_fn, chi_parity in chars:
    for basis in bases:
        for kappa in kappas:
            label = f"{chi_name} ({chi_parity}) + {basis:>8} + kappa={kappa}"
            print(f"\n=== {label} ===")
            print(f"{'lambda':>7}  {'gap':>10}  {'gap/sqrtL':>10}  {'dom':>5}")
            data = []
            for lam in lam_values:
                L = np.log(lam)
                N = min(80, max(20, int(2.0 * L * L)))
                ec, es, g = gap_at(chi_fn, lam, N, kappa, basis)
                dom = 'ODD' if g < 0 else 'EVEN'
                print(f"{lam:7d}  {g:+10.4f}  {g/np.sqrt(L):+10.4f}  {dom:>5}", flush=True)
                data.append((lam, g, g/np.sqrt(L), dom))
            n_e = sum(1 for d in data if d[3] == 'EVEN')
            n_o = sum(1 for d in data if d[3] == 'ODD')
            mean_g = np.mean([d[1] for d in data])
            summary.append((chi_name, chi_parity, basis, kappa, n_e, n_o, mean_g, data))
            print(f"[Stats] {n_e} EVEN / {n_o} ODD   mean gap = {mean_g:+.4f}")

print(f"\n{'='*80}")
print("=== KONSOLIDIERTE TABELLE ===")
print(f"{'='*80}")
print(f"{'chi':>6}  {'par(chi)':>9}  {'basis':>9}  {'kappa':>5}  {'EVEN':>4}  {'ODD':>4}  "
      f"{'mean gap':>10}  {'dom':>8}")
for chi_name, parity, basis, kappa, n_e, n_o, mean_g, _ in summary:
    classif = 'EVEN' if n_e > n_o else 'ODD' if n_o > n_e else 'MIXED'
    if abs(mean_g) < 0.1:
        classif += ' (weak)'
    print(f"{chi_name:>6}  {parity:>9}  {basis:>9}  {kappa:>5}  {n_e:4d}  {n_o:4d}  "
          f"{mean_g:+10.4f}  {classif:>8}")

print("\n[done]")
