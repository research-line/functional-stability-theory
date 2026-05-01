#!/usr/bin/env python3
# coding: utf-8
"""
Breit-Test: Dominance-Transfer fuer mehrere EVEN Dirichlet-Charaktere.

Getestete Charaktere:
  - chi_5   (Legendre mod 5,  EVEN)       [bereits in chi5_even_test.py, hier als Referenz]
  - chi_8a  (primitive mod 8,  EVEN)      Kronecker (2/n)-artig: 1,-1,-1,1 auf {1,3,5,7}
  - chi_12  (primitive mod 12, EVEN)      1,-1,-1,1 auf {1,5,7,11}
  - chi_13  (Legendre mod 13,  EVEN)      (p-1)/2 = 6, chi(-1) = +1

Ziel: Konsolidierung der Meta-Paper-Vorhersage:
  "Fuer even characters zeigt v2.1 stabile EVEN-Dominance (Gap > 0)."

Falls alle vier Charaktere EVEN-Dominance zeigen, ist die revidierte
Cartography (§5.5 Meta-Paper) empirisch gefestigt.

Autor: LG, Session 4, 2026-04-15.
"""
import numpy as np
import sympy
from scipy.special import digamma
import sys
sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# Charakter-Definitionen
# ============================================================

def chi5(n):
    """Legendre-Symbol mod 5 (EVEN, chi(-1) = +1)."""
    if n % 5 == 0: return 0
    r = n % 5
    return {1: 1, 2: -1, 3: -1, 4: 1}[r]

def chi8a(n):
    """Primitive mod 8, EVEN: chi(1)=1, chi(3)=-1, chi(5)=-1, chi(7)=1.
    Entspricht Kronecker (2/n) fuer ungerade n (Legendre-Fortsetzung)."""
    if n % 2 == 0: return 0
    r = n % 8
    return {1: 1, 3: -1, 5: -1, 7: 1}[r]

def chi12(n):
    """Primitive mod 12, EVEN: chi(1)=1, chi(5)=-1, chi(7)=-1, chi(11)=1.
    Entspricht Kronecker (12/.) = (3/.)(4/.)."""
    if np.gcd(n, 12) != 1: return 0
    r = n % 12
    return {1: 1, 5: -1, 7: -1, 11: 1}[r]

def chi13(n):
    """Legendre-Symbol mod 13 (EVEN, chi(-1) = (-1)^6 = +1).
    Quadratische Reste mod 13: {1,3,4,9,10,12}."""
    if n % 13 == 0: return 0
    # Quadratische Reste berechnen
    qr = {pow(i, 2, 13) for i in range(1, 13)}
    return 1 if (n % 13) in qr else -1

# Sanity checks
assert chi5(-1) == 1, "chi_5 muss EVEN sein"
assert chi8a(-1) == chi8a(7) == 1, "chi_8a muss EVEN sein"
assert chi12(-1) == chi12(11) == 1, "chi_12 muss EVEN sein"
assert chi13(-1) == chi13(12) == 1, "chi_13 muss EVEN sein"
# chi_13: pruefe (3/13): 3 ist QR? 3 = ?^2 mod 13: 4^2=16=3, also ja.
assert chi13(3) == 1
assert chi13(2) == -1  # 2 ist kein QR mod 13
print("[Sanity] Alle vier Charaktere sind verifiziert EVEN.")

# ============================================================
# Galerkin-Matrix-Aufbau (identisch zu chi5_even_test.py)
# ============================================================

def build_W_matrix_fast(sector, N, L, primes, chi_vals, kappa):
    if sector == 'cos':
        bs = 0
    else:
        bs = 1
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
                    overlap_matrix = I1 + I2
                else:
                    overlap_matrix = I1 - I2
                W += weight * overlap_matrix * norm_matrix
    return 0.5 * (W + W.T)

def gap_at(chi_fn, lam, N, kappa):
    L = np.log(lam)
    primes = [p for p in sympy.primerange(2, int(lam)+1) if chi_fn(p) != 0]
    chi_vals = [chi_fn(p) for p in primes]
    W_c = build_W_matrix_fast('cos', N, L, primes, chi_vals, kappa)
    W_s = build_W_matrix_fast('sin', N, L, primes, chi_vals, kappa)
    ec = np.linalg.eigvalsh(W_c); es = np.linalg.eigvalsh(W_s)
    return ec[0], es[0], es[0] - ec[0]

# ============================================================
# Main Loop: alle vier Charaktere bei selben lambda-Werten
# ============================================================

characters = [
    ('chi_5',  chi5,  'Legendre mod 5'),
    ('chi_8a', chi8a, 'primitive mod 8'),
    ('chi_12', chi12, 'primitive mod 12'),
    ('chi_13', chi13, 'Legendre mod 13'),
]

lam_values = [30, 50, 100, 200, 500, 1000, 2000, 5000]

summary = {}

for name, chi_fn, label in characters:
    print(f"\n{'='*70}")
    print(f"=== {name} ({label}, EVEN) ===")
    print(f"{'='*70}")
    print(f"{'lambda':>7}  {'N':>4}  {'lam1+':>10}  {'lam1-':>10}  {'gap':>10}  {'gap/sqrt(L)':>14}  {'dom':>5}")
    data = []
    for lam in lam_values:
        L = np.log(lam)
        N = min(80, max(20, int(2.0 * L * L)))
        ec, es, g = gap_at(chi_fn, lam, N, kappa=0)
        dom = 'ODD' if g < 0 else 'EVEN'
        print(f"{lam:7d}  {N:4d}  {ec:+10.4f}  {es:+10.4f}  {g:+10.4f}  {g/np.sqrt(L):+14.4f}  {dom:>5}",
              flush=True)
        data.append((lam, N, g, g/np.sqrt(L), dom))
    n_even = sum(1 for d in data if d[4] == 'EVEN')
    n_odd = sum(1 for d in data if d[4] == 'ODD')
    print(f"\n[Stats {name}] {n_even} EVEN / {n_odd} ODD von {len(data)}")
    summary[name] = (n_even, n_odd, data)

# ============================================================
# Konsolidierter Vergleich
# ============================================================

print(f"\n{'='*70}")
print("=== KONSOLIDIERTER VERGLEICH ===")
print(f"{'='*70}")
print(f"{'Charakter':>10}  {'EVEN':>5}  {'ODD':>5}  {'mean gap':>10}  {'median gap':>10}  {'min gap':>10}")
for name, (ne, no, data) in summary.items():
    gaps = np.array([d[2] for d in data])
    print(f"{name:>10}  {ne:5d}  {no:5d}  {np.mean(gaps):+10.4f}  {np.median(gaps):+10.4f}  {np.min(gaps):+10.4f}")

print("\n--- Gap-Skalierung mit lambda (normierter Gap / sqrt(log lambda)) ---")
print(f"{'lambda':>7}  " + "  ".join(f"{n:>10}" for n, _, _ in characters))
for i, lam in enumerate(lam_values):
    row = [f"{lam:7d}"]
    for name, _, _ in characters:
        gn = summary[name][2][i][3]  # gap/sqrt(L)
        row.append(f"{gn:+10.4f}")
    print("  ".join(row))

print("\n[Done]")
