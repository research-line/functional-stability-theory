#!/usr/bin/env python3
# coding: utf-8
"""
Server-Lauf: Asymptotik des EVEN-Dominance-Gap fuer chi_5, chi_12.

Ziel (Handoff Punkt B): ist Gap ~ 1/sqrt(log lambda) (Siegel-Walfisz)
oder ~ 1/log(lambda)?

Dazu: grosser Diskretisierungsparameter N (= 150 oder 200), lambda-Werte
bis 20_000 oder 50_000. Laufzeit pro (chi, lambda) bei N=200 geschaetzt
5-20 Min auf CCX13.

Output: JSON + TXT-Log im aktuellen Verzeichnis.

Aufruf (Server):
  PYTHONIOENCODING=utf-8 python asymptotic_scan_server.py
"""
import numpy as np
import sympy
from scipy.special import digamma
import json
import time
import sys
sys.stdout.reconfigure(line_buffering=True)

# Charaktere
def chi5(n):
    if n % 5 == 0: return 0
    return {1: 1, 2: -1, 3: -1, 4: 1}[n % 5]

def chi12(n):
    if np.gcd(n, 12) != 1: return 0
    return {1: 1, 5: -1, 7: -1, 11: 1}[n % 12]

def build_W_matrix_fast(sector, N, L, primes, chi_vals, kappa):
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
                    ovl = I1 + I2
                else:
                    ovl = I1 - I2
                W += weight * ovl * norm_matrix
    return 0.5 * (W + W.T)

def gap_at(chi_fn, lam, N, kappa):
    L = np.log(lam)
    primes = [p for p in sympy.primerange(2, int(lam)+1) if chi_fn(p) != 0]
    chi_vals = [chi_fn(p) for p in primes]
    W_c = build_W_matrix_fast('cos', N, L, primes, chi_vals, kappa)
    W_s = build_W_matrix_fast('sin', N, L, primes, chi_vals, kappa)
    ec = np.linalg.eigvalsh(W_c); es = np.linalg.eigvalsh(W_s)
    return float(ec[0]), float(es[0]), float(es[0] - ec[0])

# ============================================================
# Lauf-Konfiguration
# ============================================================

lam_values = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]
N_values = [80, 120, 160, 200]

characters = [
    ('chi_5',  chi5),
    ('chi_12', chi12),
]

results = []

t0 = time.time()
for name, chi_fn in characters:
    for N in N_values:
        for lam in lam_values:
            L = np.log(lam)
            # Kriterium: N >= 2L^2 (aus Session-3-Erfahrung)
            required_N = int(2 * L * L)
            if N < required_N * 0.5:
                print(f"[skip] {name} lam={lam} N={N} (req~{required_N})", flush=True)
                continue
            t_start = time.time()
            try:
                ec, es, g = gap_at(chi_fn, lam, N, kappa=0)
                dt = time.time() - t_start
                row = {
                    'chi': name, 'lambda': lam, 'N': N, 'L': L,
                    'lam1_cos': ec, 'lam1_sin': es, 'gap': g,
                    'gap_norm_sqrtL': g / np.sqrt(L),
                    'gap_norm_L': g / L,
                    'runtime_s': dt,
                }
                results.append(row)
                print(f"[{name}] N={N:3d} lam={lam:6d} gap={g:+10.5f} "
                      f"gap/sqrtL={g/np.sqrt(L):+8.5f} gap/L={g/L:+8.5f} t={dt:5.1f}s",
                      flush=True)
                # Zwischensichern
                with open('asymptotic_results.json', 'w') as f:
                    json.dump(results, f, indent=2)
            except Exception as e:
                print(f"[error] {name} N={N} lam={lam}: {e}", flush=True)

print(f"\n[done] total elapsed = {(time.time()-t0)/60:.1f} min; "
      f"n_rows = {len(results)}", flush=True)

# ============================================================
# Asymptotik-Auswertung: gap ~ lambda^alpha * (log lambda)^beta
# ============================================================

print("\n=== Asymptotik-Regression (log-log) ===")
for name in ['chi_5', 'chi_12']:
    print(f"\n--- {name} ---")
    for N in N_values:
        sub = [r for r in results if r['chi'] == name and r['N'] == N and r['gap'] > 0]
        if len(sub) < 3: continue
        xs = np.log([r['lambda'] for r in sub])
        ys = np.log([r['gap'] for r in sub])
        # Fit: log gap = a * log lambda + b
        A = np.polyfit(xs, ys, 1)
        print(f"  N={N}: log(gap) = {A[0]:+.4f} * log(lam) + {A[1]:+.4f}   "
              f"(n={len(sub)})")

print("\n[End]")
