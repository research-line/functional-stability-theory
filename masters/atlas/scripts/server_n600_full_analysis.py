#!/usr/bin/env python3
# coding: utf-8
"""
Session 7 Woche 2: N=600 Full Arch-Term-Analyse aller 10 Charaktere.

Standalone (keine Cross-Imports), zum Lauf auf ellmos-services.

Output: ARCH_TERM_N600_SERVER.json + Fortschrittslog.

Geschaetzte Laufzeit: ~15-20 min pro Sektor bei N=600 (lambda=20000),
10 Charaktere × 2 Sektoren = 5-7 Stunden Gesamtlaufzeit.
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import sympy
from scipy.special import digamma

sys.stdout.reconfigure(line_buffering=True)

# --------------------------------------------------------------------------
# Konfiguration
# --------------------------------------------------------------------------
LAM = 20000
N_GALERKIN = 600
OUT_DIR = Path(__file__).resolve().parent / "results"
OUT_JSON = OUT_DIR / "ARCH_TERM_N600_SERVER.json"
PROGRESS = OUT_DIR / "progress.log"

CHARS = [
    ("chi_5", 5),
    ("chi_8", 8),
    ("chi_12", 12),
    ("chi_13", 13),
    ("chi_17", 17),
    ("chi_21", 21),
    ("chi_24", 24),
    ("chi_29", 29),
    ("chi_33", 33),
    ("chi_60", 60),
]


# --------------------------------------------------------------------------
# Kronecker-Symbol / Dirichlet-Charakter
# --------------------------------------------------------------------------
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


def make_chi_D(D):
    def chi(n):
        if math.gcd(n, abs(D)) != 1:
            return 0
        return kronecker_symbol(D, n)
    return chi


# --------------------------------------------------------------------------
# Weil-Kern-Matrix build_W (kappa=0)
# --------------------------------------------------------------------------
def build_W(sector, N, L, primes, chi_vals, kappa=0):
    bs = 0 if sector == 'cos' else 1
    idx = np.arange(N) + bs
    nn = idx[:, None]
    mm = idx[None, :]
    W = np.zeros((N, N))
    shift = (2 * kappa + 1) / 4.0
    tau = np.pi * idx / L
    diag = np.array([digamma(shift + 1j * t / 2).real - np.log(np.pi) for t in tau])
    np.fill_diagonal(W, diag)
    if sector == 'cos':
        norm = np.where(idx > 0, 1.0 / np.sqrt(L), 1.0 / np.sqrt(2 * L))
    else:
        norm = np.full(N, 1.0 / np.sqrt(L))
    NM = np.outer(norm, norm)

    def A(k, phi, t, Lv):
        r = np.zeros_like(k, dtype=float)
        mk = (k != 0)
        r[mk] = (Lv / (k[mk] * np.pi)) * np.sin(k[mk] * np.pi * t / Lv + phi[mk])
        r[~mk] = np.cos(phi[~mk]) * t
        return r

    for p, cp in zip(primes, chi_vals):
        if cp == 0:
            continue
        lp = np.log(p)
        mmx = int(2 * L / lp)
        for me in range(1, mmx + 1):
            d = me * lp
            if d >= 2 * L:
                break
            w = (cp ** me) * lp / (p ** (me / 2.0))
            for ds in [d, -d]:
                a = max(-L, -L + ds)
                b = min(L, L + ds)
                if a >= b:
                    continue
                k1 = nn - mm
                k2 = nn + mm
                ph1 = mm * np.pi * ds / L
                ph2 = -mm * np.pi * ds / L
                P1 = np.broadcast_to(ph1, (N, N)).copy()
                P2 = np.broadcast_to(ph2, (N, N)).copy()
                I1 = 0.5 * (A(k1, P1, b, L) - A(k1, P1, a, L))
                I2 = 0.5 * (A(k2, P2, b, L) - A(k2, P2, a, L))
                ovl = I1 + I2 if sector == 'cos' else I1 - I2
                W += w * ovl * NM
    return 0.5 * (W + W.T)


def build_W_arch_only(sector, N, L, kappa=0):
    bs = 0 if sector == 'cos' else 1
    idx = np.arange(N) + bs
    tau = np.pi * idx / L
    shift = (2 * kappa + 1) / 4.0
    W = np.zeros((N, N))
    diag = np.array([digamma(shift + 1j * t / 2).real - np.log(np.pi) for t in tau])
    np.fill_diagonal(W, diag)
    return W


def compute_ground_state(chi_fn, lam, N=200, sector='cos'):
    L = np.log(lam)
    primes = [p for p in sympy.primerange(2, int(lam) + 1) if chi_fn(p) != 0]
    cv = [chi_fn(p) for p in primes]
    W = build_W(sector, N, L, primes, cv)
    eigvals, eigvecs = np.linalg.eigh(W)
    return float(eigvals[0]), eigvecs[:, 0]


def phi_from_coefs(coefs, L, sector='cos', n_grid=1600):
    xs = np.linspace(-L, L, n_grid)
    bs = 0 if sector == 'cos' else 1
    N = len(coefs)
    idx = np.arange(N) + bs
    phi = np.zeros_like(xs)
    for n, c in zip(idx, coefs):
        if sector == 'cos':
            if n == 0:
                phi += c / math.sqrt(2 * L)
            else:
                phi += c * np.cos(n * np.pi * xs / L) / math.sqrt(L)
        else:
            phi += c * np.sin(n * np.pi * xs / L) / math.sqrt(L)
    return xs, phi


def autocorr_at(phi, xs, ts):
    dx = xs[1] - xs[0]
    ac = np.correlate(phi, phi, mode='full') * dx
    N = len(phi)
    lags = np.arange(-(N - 1), N) * dx
    return np.interp(ts, lags, ac, left=0.0, right=0.0)


def compute_gap_exact(chi_fn, lam, phi_plus, xs_plus, phi_minus, xs_minus,
                     max_m=10, coeff=-2.0):
    L = math.log(lam)
    twoL = 2.0 * L
    S = 0.0
    for p in sympy.primerange(2, int(lam) + 1):
        cp = chi_fn(p)
        if cp == 0:
            continue
        log_p = math.log(p)
        if log_p > twoL:
            break
        for m in range(1, max_m + 1):
            t = m * log_p
            if t > twoL:
                break
            phase = cp ** m
            weight = phase * log_p / (p ** (m / 2.0))
            ac_plus = float(autocorr_at(phi_plus, xs_plus, np.array([t]))[0])
            ac_minus = float(autocorr_at(phi_minus, xs_minus, np.array([t]))[0])
            delta_t = ac_plus - ac_minus
            S += weight * delta_t
    return coeff * S


def arch_contrib(coefs, sector, L, kappa=0):
    N = len(coefs)
    W_arch = build_W_arch_only(sector, N, L, kappa)
    return float(coefs @ W_arch @ coefs)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(exist_ok=True)
    start = time.time()
    L = math.log(LAM)

    def log(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line)
        with open(PROGRESS, "a") as f:
            f.write(line + "\n")

    log("=" * 75)
    log(f"Server Full Arch-Term-Analysis: lambda={LAM}, N={N_GALERKIN}")
    log("=" * 75)

    results = []
    for i, (name, D) in enumerate(CHARS, 1):
        log(f"\n[{i}/{len(CHARS)}] chi={name}, D={D}")
        chi_fn = make_chi_D(D)

        t0 = time.time()
        ev_plus, coefs_plus = compute_ground_state(chi_fn, LAM, N=N_GALERKIN, sector='cos')
        t_plus = time.time() - t0
        log(f"  cos done: ev+={ev_plus:+.5f} [{t_plus:.0f}s]")

        t0 = time.time()
        ev_minus, coefs_minus = compute_ground_state(chi_fn, LAM, N=N_GALERKIN, sector='sin')
        t_minus = time.time() - t0
        log(f"  sin done: ev-={ev_minus:+.5f} [{t_minus:.0f}s]")

        arch_plus = arch_contrib(coefs_plus, 'cos', L)
        arch_minus = arch_contrib(coefs_minus, 'sin', L)
        arch_diff = arch_minus - arch_plus

        xs_p, phi_p = phi_from_coefs(coefs_plus, L, sector='cos', n_grid=1600)
        xs_m, phi_m = phi_from_coefs(coefs_minus, L, sector='sin', n_grid=1600)
        S_exact = compute_gap_exact(chi_fn, LAM, phi_p, xs_p, phi_m, xs_m, coeff=-2.0)

        gap_galerkin = ev_minus - ev_plus
        S_with_arch = S_exact + arch_diff

        log(f"  gap_gal={gap_galerkin:+.5f}, S={S_exact:+.4e}, arch={arch_diff:+.4e}, S+arch={S_with_arch:+.4e}")

        results.append({
            "chi": name,
            "D": D,
            "lambda": LAM,
            "N_Galerkin": N_GALERKIN,
            "ev_plus": float(ev_plus),
            "ev_minus": float(ev_minus),
            "gap_galerkin": float(gap_galerkin),
            "arch_plus": float(arch_plus),
            "arch_minus": float(arch_minus),
            "arch_diff": float(arch_diff),
            "S_exact": float(S_exact),
            "S_with_arch": float(S_with_arch),
            "walltime_plus_s": float(t_plus),
            "walltime_minus_s": float(t_minus),
        })

        # Inkrementell speichern, falls Abbruch
        with open(OUT_JSON, "w") as f:
            json.dump({
                "parameters": {"lambda": LAM, "N_Galerkin": N_GALERKIN, "coeff": -2.0},
                "results": results,
                "completed": len(results),
                "total": len(CHARS),
                "start_time": start,
                "elapsed_s": time.time() - start,
            }, f, indent=2)

    total = time.time() - start
    log(f"\nGesamtlaufzeit: {total:.0f}s = {total/60:.1f} min = {total/3600:.2f} h")
    log(f"Ergebnisse: {OUT_JSON}")


if __name__ == "__main__":
    main()
