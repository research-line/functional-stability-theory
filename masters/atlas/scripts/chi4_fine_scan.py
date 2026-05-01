#!/usr/bin/env python3
# coding: utf-8
"""
Feiner Scan: Gap(lambda) fuer chi_4 bei vielen lambda-Werten.
Ziel: charakterisieren des Alternations-Phaenomens.

Zusaetzlich: Korrelation mit Chebyshev-Bias pi(x; 4, 3) - pi(x; 4, 1).
"""
import numpy as np
import sympy
from scipy.special import digamma
import sys
sys.stdout.reconfigure(line_buffering=True)

# --- Copy vectorized build_W function from chi4_vectorized.py ---

def chi4(n):
    r = n % 4
    if r == 1: return 1
    if r == 3: return -1
    return 0

def build_W_matrix_fast(sector, N, L, primes, chi_vals, kappa):
    if sector == 'cos':
        bs = 0
    else:
        bs = 1
    idx = np.arange(N) + bs
    nn = idx[:, None]
    mm = idx[None, :]

    W = np.zeros((N, N))
    shift = (2*kappa+1)/4.0
    tau = np.pi * idx / L
    diag = np.array([digamma(shift + 1j*t/2).real - np.log(np.pi) for t in tau])
    np.fill_diagonal(W, diag)

    if sector == 'cos':
        norm = np.where(idx > 0, 1.0/np.sqrt(L), 1.0/np.sqrt(2*L))
    else:
        norm = np.full(N, 1.0/np.sqrt(L))
    norm_matrix = np.outer(norm, norm)

    def A_manual(k_arr, phi_arr, t, L):
        result = np.zeros_like(k_arr, dtype=float)
        mask = (k_arr != 0)
        kk = k_arr[mask]; pp = phi_arr[mask]
        result[mask] = (L / (kk * np.pi)) * np.sin(kk * np.pi * t / L + pp)
        result[~mask] = np.cos(phi_arr[~mask]) * t
        return result

    for p, cp in zip(primes, chi_vals):
        if cp == 0: continue
        lp = np.log(p)
        mm_max = int(2*L/lp)
        for me in range(1, mm_max+1):
            d = me * lp
            if d >= 2*L: break
            weight = (cp**me) * lp / (p**(me/2.0))
            for d_signed in [d, -d]:
                a = max(-L, -L + d_signed)
                b = min(L, L + d_signed)
                if a >= b: continue
                k1 = nn - mm; k2 = nn + mm
                phi1 = mm * np.pi * d_signed / L
                phi2 = -mm * np.pi * d_signed / L
                phi1_full = np.broadcast_to(phi1, (N, N)).copy()
                phi2_full = np.broadcast_to(phi2, (N, N)).copy()
                I1 = 0.5*(A_manual(k1, phi1_full, b, L) - A_manual(k1, phi1_full, a, L))
                I2 = 0.5*(A_manual(k2, phi2_full, b, L) - A_manual(k2, phi2_full, a, L))
                if sector == 'cos':
                    overlap_matrix = I1 + I2
                else:
                    overlap_matrix = I1 - I2
                W += weight * overlap_matrix * norm_matrix
    return 0.5 * (W + W.T)

def gap_at(lam, N, kappa=1):
    L = np.log(lam)
    primes = list(sympy.primerange(3, int(lam)+1))
    chi_vals = [chi4(p) for p in primes]
    W_c = build_W_matrix_fast('cos', N, L, primes, chi_vals, kappa)
    W_s = build_W_matrix_fast('sin', N, L, primes, chi_vals, kappa)
    ec = np.linalg.eigvalsh(W_c); es = np.linalg.eigvalsh(W_s)
    return ec[0], es[0], es[0] - ec[0]

def chebyshev_bias(x):
    """Zaehlen: #{p <= x, p = 3 mod 4} - #{p <= x, p = 1 mod 4}"""
    primes = list(sympy.primerange(3, int(x)+1))
    n3 = sum(1 for p in primes if p % 4 == 3)
    n1 = sum(1 for p in primes if p % 4 == 1)
    return n3 - n1  # positiv = Chebyshev-Bias in der bekannten Richtung

# ============================================================
# Feiner Scan
# ============================================================
print(f"{'lambda':>7}  {'N':>4}  {'gap':>10}  {'gap/sqrt(L)':>12}  {'pi(3)-pi(1)':>12}  {'dom':>5}")
lam_values = [30, 50, 75, 100, 150, 200, 250, 300, 400, 500, 700, 1000, 1400, 2000, 3000, 5000]
data = []
for lam in lam_values:
    L = np.log(lam)
    N = min(80, max(20, int(2.0 * L * L)))
    ec, es, g = gap_at(lam, N, kappa=1)
    bias = chebyshev_bias(lam)
    dom = 'ODD' if g < 0 else 'EVEN'
    print(f"{lam:7d}  {N:4d}  {g:+10.4f}  {g/np.sqrt(L):+12.4f}  {bias:+12d}  {dom:>5}", flush=True)
    data.append((lam, N, g, bias, dom))

# ============================================================
# Korrelation: Gap ↔ Chebyshev-Bias
# ============================================================
gaps = np.array([d[2] for d in data])
biases = np.array([d[3] for d in data])
lams = np.array([d[0] for d in data])

correlation = np.corrcoef(gaps, biases)[0, 1]
print(f"\n=== Korrelation ===")
print(f"corr(gap, pi_3-pi_1) = {correlation:+.4f}")

# ============================================================
# Zusatz: Dominance-Statistik
# ============================================================
n_odd = sum(1 for d in data if d[4] == 'ODD')
n_even = sum(1 for d in data if d[4] == 'EVEN')
print(f"\n=== Dominance-Statistik ===")
print(f"ODD-Dominance:  {n_odd} / {len(data)} Faellen")
print(f"EVEN-Dominance: {n_even} / {len(data)} Faellen")

# ============================================================
# Log-Log Plot der Gap-Absolutwerte
# ============================================================
print(f"\n=== Gap-Absolutwert vs. lambda ===")
print(f"{'lambda':>7}  {'|gap|':>10}  {'log|gap|/log(lambda)':>25}")
for lam, N, g, bias, dom in data:
    print(f"{lam:7d}  {abs(g):10.4f}  {np.log(abs(g))/np.log(lam):+25.4f}")
