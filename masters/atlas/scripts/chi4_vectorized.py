#!/usr/bin/env python3
# coding: utf-8
"""
Vektorisierte Version: Gap-Stabilitaet ueber N und lambda fuer chi_4.
Nutzt numpy-Broadcasting fuer Overlap-Matrix-Berechnung.
Laufzeit: ~10s fuer grosse Parameter.
"""
import numpy as np
import sympy
from scipy.special import digamma
import sys
sys.stdout.reconfigure(line_buffering=True)

def chi4(n):
    r = n % 4
    if r == 1: return 1
    if r == 3: return -1
    return 0

def build_W_matrix_fast(sector, N, L, primes, chi_vals, kappa):
    """Vektorisierte Galerkin-Matrix fuer W^{chi}.
    Overlap-Integrale werden matrix-weise via numpy-Broadcasting berechnet."""
    if sector == 'cos':
        bs = 0
    else:
        bs = 1
    # Basis-Indizes
    idx = np.arange(N) + bs  # Shape (N,)
    nn = idx[:, None]  # Shape (N, 1): n-Indizes als Spalten
    mm = idx[None, :]  # Shape (1, N): m-Indizes als Zeilen

    # Archimedischer Diagonal-Anteil
    W = np.zeros((N, N))
    shift = (2*kappa+1)/4.0
    tau = np.pi * idx / L
    diag = np.array([digamma(shift + 1j*t/2).real - np.log(np.pi) for t in tau])
    np.fill_diagonal(W, diag)

    # Normierung
    norm = np.where(idx > 0, 1.0/np.sqrt(L), 1.0/np.sqrt(2*L))  # Shape (N,)
    # Fuer sin-basis: immer 1/sqrt(L) (weil idx startet bei 1)
    if sector == 'sin':
        norm = np.full(N, 1.0/np.sqrt(L))
    norm_matrix = np.outer(norm, norm)  # Shape (N, N)

    # Prime-Shift-Schleife: ueber (p, m_exp, Vorzeichen) → 3D
    for p, cp in zip(primes, chi_vals):
        if cp == 0:
            continue
        lp = np.log(p)
        mm_max = int(2*L/lp)
        for me in range(1, mm_max+1):
            d = me * lp
            if d >= 2*L: break
            weight = (cp**me) * lp / (p**(me/2.0))
            # Fuer beide Vorzeichen delta und -delta:
            for d_signed in [d, -d]:
                # Integrationsbereich
                a = max(-L, -L + d_signed)
                b = min(L, L + d_signed)
                if a >= b: continue

                # Schluessel-Berechnung: Overlap-Matrix
                if sector == 'cos':
                    # S_{n,m}(delta) = (1/2) int_a^b [cos((n-m) pi t/L + m*pi*d/L) + cos((n+m) pi t/L - m*pi*d/L)] dt
                    # = (1/2) [A_1(n-m, m*d; a->b) + A_2(n+m, -m*d; a->b)]
                    # wobei A(k, phi; a->b) = int_a^b cos(k pi t/L + phi) dt
                    # Fuer k != 0: (L/(k pi)) (sin(k pi b/L + phi) - sin(k pi a/L + phi))
                    # Fuer k == 0: cos(phi) (b - a)
                    k1 = nn - mm  # Shape (N, N), n - m
                    k2 = nn + mm  # Shape (N, N), n + m
                    phi1 = mm * np.pi * d_signed / L  # depends on m only
                    phi2 = -mm * np.pi * d_signed / L

                    def A(k, phi, t):
                        # k, phi, t can be arrays or scalars
                        k = np.asarray(k); phi = np.asarray(phi)
                        result = np.zeros_like(k, dtype=float)
                        mask = k != 0
                        result[mask] = (L / (k[mask] * np.pi)) * np.sin(k[mask] * np.pi * t / L + phi[mask] if phi.shape == k.shape else k[mask] * np.pi * t / L + phi)
                        result[~mask] = np.cos(phi[~mask] if phi.shape == k.shape else phi) * t
                        return result

                    # Simplify: broadcast manually
                    def A_manual(k_arr, phi_arr, t):
                        # k_arr, phi_arr: (N, N) arrays; t: scalar
                        result = np.zeros_like(k_arr, dtype=float)
                        mask = (k_arr != 0)
                        # fuer k != 0: (L/(k pi)) sin(k pi t/L + phi)
                        kk = k_arr[mask]
                        pp = phi_arr[mask]
                        result[mask] = (L / (kk * np.pi)) * np.sin(kk * np.pi * t / L + pp)
                        # fuer k == 0: cos(phi) * t
                        result[~mask] = np.cos(phi_arr[~mask]) * t

                        return result

                    # phi1 ist m-abhaengig (broadcast). Erweitere:
                    phi1_full = np.broadcast_to(phi1, (N, N)).copy()
                    phi2_full = np.broadcast_to(phi2, (N, N)).copy()

                    I1 = 0.5 * (A_manual(k1, phi1_full, b) - A_manual(k1, phi1_full, a))
                    I2 = 0.5 * (A_manual(k2, phi2_full, b) - A_manual(k2, phi2_full, a))
                    overlap_matrix = I1 + I2

                else:
                    # sin-basis: sin(A)sin(B) = (1/2)[cos(A-B) - cos(A+B)]
                    k1 = nn - mm
                    k2 = nn + mm
                    phi1 = mm * np.pi * d_signed / L
                    phi2 = -mm * np.pi * d_signed / L

                    def A_manual(k_arr, phi_arr, t):
                        result = np.zeros_like(k_arr, dtype=float)
                        mask = (k_arr != 0)
                        kk = k_arr[mask]
                        pp = phi_arr[mask]
                        result[mask] = (L / (kk * np.pi)) * np.sin(kk * np.pi * t / L + pp)
                        result[~mask] = np.cos(phi_arr[~mask]) * t
                        return result

                    phi1_full = np.broadcast_to(phi1, (N, N)).copy()
                    phi2_full = np.broadcast_to(phi2, (N, N)).copy()

                    I1 = 0.5 * (A_manual(k1, phi1_full, b) - A_manual(k1, phi1_full, a))
                    I2 = 0.5 * (A_manual(k2, phi2_full, b) - A_manual(k2, phi2_full, a))
                    overlap_matrix = I1 - I2

                # Normalisierung
                overlap_matrix = overlap_matrix * norm_matrix
                W += weight * overlap_matrix

    return 0.5 * (W + W.T)

def gap_at(lam, N, kappa=1):
    L = np.log(lam)
    primes = list(sympy.primerange(3, int(lam)+1))
    chi_vals = [chi4(p) for p in primes]
    W_c = build_W_matrix_fast('cos', N, L, primes, chi_vals, kappa)
    W_s = build_W_matrix_fast('sin', N, L, primes, chi_vals, kappa)
    ec = np.linalg.eigvalsh(W_c); es = np.linalg.eigvalsh(W_s)
    return ec[0], es[0], es[0] - ec[0]

print("=== Truncation-Stabilitaet: Gap als Funktion von N bei lambda = 1000 ===")
print(f"{'N':>4}  {'lam1+':>10}  {'lam1-':>10}  {'gap':>10}  {'dominance':>12}")
for N in [10, 14, 18, 22, 26, 30, 35, 40, 50, 60, 75]:
    ec, es, g = gap_at(1000, N, kappa=1)
    dom = 'ODD' if g < 0 else 'EVEN'
    print(f"{N:4d}  {ec:+10.4f}  {es:+10.4f}  {g:+10.4f}  {dom:>12}", flush=True)

print("\n=== Skalierung: Gap ueber lambda bei N ~ 2L^2 ===", flush=True)
print(f"{'lambda':>7}  {'N':>4}  {'lam1+':>10}  {'lam1-':>10}  {'gap':>10}  {'gap/sqrt(L)':>14}  {'dominance':>10}")
for lam in [30, 60, 100, 200, 500, 1000, 2000, 5000]:
    L = np.log(lam)
    N = max(20, int(2.0 * L * L))
    N = min(N, 80)
    ec, es, g = gap_at(lam, N, kappa=1)
    dom = 'ODD' if g < 0 else 'EVEN'
    print(f"{lam:7d}  {N:4d}  {ec:+10.4f}  {es:+10.4f}  {g:+10.4f}  {g/np.sqrt(L):+14.4f}  {dom:>10}", flush=True)
