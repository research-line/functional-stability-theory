#!/usr/bin/env python3
# coding: utf-8
"""
Twisted-Basis-Test fuer odd characters (chi_4) und Kontrolle (chi_5 even).

Motivation: Die falsifizierte Vorhersage in SHIFT_PARITY_TWIST.md wird ersetzt
durch eine neue: odd characters zeigen EVEN-Dominance in der TWISTED basis
mit halbzahligen Wellenzahlen omega_n = pi(n+1/2)/L.

Bases:
  standard:  psi_n^+ = cos(pi n t/L), psi_n^- = sin(pi n t/L)
  twisted:   psi_n^+ = cos(pi(n+1/2) t/L), psi_n^- = sin(pi(n+1/2) t/L)

Gamma-Faktoren (archimedische Diagonale):
  kappa = 0 (even-Gamma):  h(tau) = Re psi_dig(1/4 + i tau/2) - log pi
  kappa = 1 (odd-Gamma):   h(tau) = Re psi_dig(3/4 + i tau/2) - log pi

Theorie-Vorhersage:
  chi_4 (odd) + twisted basis + kappa=1  =>  Even-Dominance (stabil, wachsend).
  chi_5 (even) + twisted basis + kappa=0 =>  kein Muster (Kontrolle).
  chi_5 (even) + standard basis + kappa=0 => Even-Dom (bereits bestaetigt).
  chi_4 (odd) + standard basis + kappa=1 => kein Muster (bereits beobachtet).
"""
import numpy as np
import sympy
from scipy.special import digamma
import sys
sys.stdout.reconfigure(line_buffering=True)

# --- Charaktere ---
def chi4(n):
    if n % 2 == 0: return 0
    return 1 if (n % 4) == 1 else -1

def chi5(n):
    if n % 5 == 0: return 0
    return {1: 1, 2: -1, 3: -1, 4: 1}[n % 5]

# --- Basis-unabhaengiger Baustein ---
def build_W_matrix(sector, basis, N, L, primes, chi_vals, kappa):
    """
    sector: 'cos' | 'sin'
    basis:  'standard' | 'twisted'
      standard: omega_n = pi * n / L  (cos: n >= 0, sin: n >= 1)
      twisted:  omega_n = pi * (n+1/2) / L  (both sectors: n >= 0)

    Weil-Matrix auf Galerkin-Raum Spann(psi_0, ..., psi_{N-1}).
    """
    if basis == 'standard':
        # cos: n = 0, 1, ..., N-1;  sin: n = 1, 2, ..., N
        bs = 0 if sector == 'cos' else 1
        idx = np.arange(N) + bs        # Wellenzahl-Indizes
        omega = np.pi * idx / L
    elif basis == 'twisted':
        idx = np.arange(N)             # n = 0, 1, ..., N-1
        half = np.arange(N) + 0.5      # Wellenzahl-Indizes (halbzahlig)
        omega = np.pi * half / L
    else:
        raise ValueError(basis)

    # --- Diagonale (archimedisch) ---
    W = np.zeros((N, N))
    shift = (2*kappa+1) / 4.0
    diag = np.array([digamma(shift + 1j*t/2).real - np.log(np.pi) for t in omega])
    np.fill_diagonal(W, diag)

    # --- Norm: 1/sqrt(L) fuer sin/twisted, 1/sqrt(2L) fuer standard cos n=0 ---
    if basis == 'standard' and sector == 'cos':
        norm = np.where(idx > 0, 1.0/np.sqrt(L), 1.0/np.sqrt(2*L))
    else:
        norm = np.full(N, 1.0/np.sqrt(L))
    norm_matrix = np.outer(norm, norm)

    def A(k_arr, phi_arr, t, Lval):
        """Antiderivative von cos(pi k t / L + phi): L/(pi k) sin(.) fuer k!=0, sonst cos(phi)*t."""
        result = np.zeros_like(k_arr, dtype=float)
        mask = (k_arr != 0)
        kk = k_arr[mask].astype(float)
        pp = phi_arr[mask]
        result[mask] = (Lval / (kk * np.pi)) * np.sin(kk * np.pi * t / Lval + pp)
        result[~mask] = np.cos(phi_arr[~mask]) * t
        return result

    # --- Prime-Overlaps ---
    if basis == 'standard':
        idx_n = idx  # n-Wellenzahlen
    else:
        idx_n = idx + 0  # fuer Twisted: idx_n = n, wir verwenden aber halbzahlig

    # nn, mm als 2D-Indizes:
    if basis == 'standard':
        nn = idx[:, None]  # Wellenzahlen (n)
        mm = idx[None, :]  # Wellenzahlen (m)
    else:
        nn = (idx[:, None]).astype(float) + 0.5   # halbzahlig
        mm = (idx[None, :]).astype(float) + 0.5

    # In beiden Faellen: bei Produkt-zu-Summe sind die k-Werte ganzzahlig:
    #   standard: k1 = n - m, k2 = n + m
    #   twisted:  k1 = n - m (ganzzahlig, da Diff halbzahlig - halbzahlig = ganz),
    #             k2 = n + m (= (n_int + m_int + 1), ganzzahlig)
    # Wir koennen daher direkt mit float rechnen und auf int runden:
    k1f = nn - mm
    k2f = nn + mm
    k1 = np.rint(k1f).astype(int)
    k2 = np.rint(k2f).astype(int)
    # Sanity check: alle Differenzen < 1e-10
    assert np.allclose(k1, k1f, atol=1e-10) and np.allclose(k2, k2f, atol=1e-10), \
        f"k-Werte nicht ganzzahlig: {np.max(np.abs(k1-k1f))}, {np.max(np.abs(k2-k2f))}"

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
                # Phase-Faktoren: Overlap integriert psi_n(t) * psi_m(t - d).
                # cos-cos: 1/2 [cos((omega_n-omega_m) t + omega_m d) + cos((omega_n+omega_m) t - omega_m d)]
                # sin-sin: 1/2 [cos((omega_n-omega_m) t + omega_m d) - cos((omega_n+omega_m) t - omega_m d)]
                phi1 = mm * np.pi * d_signed / L       # = omega_m * d (in pi/L units)
                phi2 = -mm * np.pi * d_signed / L      # = -omega_m * d
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

# ============================================================
# Test-Matrix (4 Kombinationen)
# ============================================================

lam_values = [30, 50, 100, 200, 500, 1000, 2000, 5000]

tests = [
    ('chi_4 odd  + twisted  + kappa=1', chi4, 1, 'twisted'),
    ('chi_5 even + twisted  + kappa=0', chi5, 0, 'twisted'),
    ('chi_4 odd  + standard + kappa=1', chi4, 1, 'standard'),
    ('chi_5 even + standard + kappa=0', chi5, 0, 'standard'),
]

summary = {}
for label, chi_fn, kappa, basis in tests:
    print(f"\n{'='*78}")
    print(f"=== {label} ===")
    print(f"{'='*78}")
    print(f"{'lambda':>7}  {'N':>4}  {'lam1+':>10}  {'lam1-':>10}  "
          f"{'gap':>10}  {'gap/sqrt(L)':>14}  {'dom':>5}")
    data = []
    for lam in lam_values:
        L = np.log(lam)
        N = min(80, max(20, int(2.0 * L * L)))
        ec, es, g = gap_at(chi_fn, lam, N, kappa, basis)
        dom = 'ODD' if g < 0 else 'EVEN'
        print(f"{lam:7d}  {N:4d}  {ec:+10.4f}  {es:+10.4f}  {g:+10.4f}  "
              f"{g/np.sqrt(L):+14.4f}  {dom:>5}", flush=True)
        data.append((lam, g, g/np.sqrt(L), dom))
    n_even = sum(1 for d in data if d[3] == 'EVEN')
    n_odd  = sum(1 for d in data if d[3] == 'ODD')
    print(f"[Stats] {n_even} EVEN / {n_odd} ODD")
    summary[label] = (n_even, n_odd, data)

print(f"\n{'='*78}")
print("=== SELEKTIONSREGEL-TEST (Basis x Charakter) ===")
print(f"{'='*78}")
print(f"Theorie: EVEN-Dom gdw 'matching' Basis.")
print(f"  chi_4 (odd)  natural basis = TWISTED")
print(f"  chi_5 (even) natural basis = STANDARD\n")
print(f"{'Kombination':>40}  {'EVEN':>5}  {'ODD':>5}  {'mean gap':>10}")
for label, (ne, no, data) in summary.items():
    mean_gap = np.mean([d[1] for d in data])
    print(f"{label:>40}  {ne:5d}  {no:5d}  {mean_gap:+10.4f}")

print("\n[done]")
