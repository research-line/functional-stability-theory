#!/usr/bin/env python3
# coding: utf-8
"""
Woche-4-Validierung v2: Teste mehrere Fenster-Formen plus **empirischen Grundzustand**.

Nach v1-Scheitern des naiven Gauss-Ansatzes (nur 6/10 sign_ok, R^2=0.15) teste:
  1. Gauss-Fenster verschiedener Breiten und Zentren (Erweiterung von v1)
  2. Rechteck-Fenster [0, T] fuer verschiedene T (Cutoff-Analyse)
  3. **Empirischer Grundzustand**: berechne Eigenvektor via Galerkin, autokorreliere, nutze als Phi_inf
  4. Raw chi-Summe als Baseline: S = sum_{p<=x} chi(p) log p / sqrt p
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import sympy
from scipy.special import digamma

sys.stdout.reconfigure(line_buffering=True)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "_results"

# Reuse aus v1
from analytic_gap_formula_test import (
    kronecker_symbol, make_chi_D, CHARS, load_empirical_gaps, S_chi, gauss_window
)


# --------------------------------------------------------------------------
# Alternative Fenster-Formen
# --------------------------------------------------------------------------
def S_chi_window(chi_fn, lam, window_fn, max_m=10):
    """Generalized S_chi fuer beliebige Fensterfunktion window_fn(t)."""
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
        p_half = math.sqrt(p)
        for m in range(1, max_m + 1):
            t = m * log_p
            if t > twoL:
                break
            phase = cp ** m
            weight = phase * log_p / (p ** (m / 2.0))
            S += weight * window_fn(t)
    return S


def rect_window(t, T):
    """Rechteck-Fenster [0, T]."""
    return 1.0 if 0 <= t <= T else 0.0


def gauss_offset(t, sigma, t0):
    """Gauss um t0."""
    return math.exp(-((t - t0) ** 2) / (4 * sigma ** 2)) / (math.sqrt(2 * math.pi) * sigma)


def cos_window(t, T):
    """Cosine / prolate-like: cos^2(pi t / (2T)) auf [-T, T]."""
    if abs(t) > T:
        return 0.0
    return math.cos(math.pi * t / (2 * T)) ** 2


# --------------------------------------------------------------------------
# Empirischer Grundzustand via Galerkin-Eigenvektor
# --------------------------------------------------------------------------
def build_W(sector, N, L, primes, chi_vals, kappa=0):
    """Identisch zu all10_high_N_server.py."""
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


def compute_ground_state(chi_fn, lam, N=200, sector='cos'):
    """Berechne Grundzustand (Eigenvektor zum niedrigsten EW) fuer Sektor."""
    L = np.log(lam)
    primes = [p for p in sympy.primerange(2, int(lam) + 1) if chi_fn(p) != 0]
    cv = [chi_fn(p) for p in primes]
    W = build_W(sector, N, L, primes, cv)
    eigvals, eigvecs = np.linalg.eigh(W)
    # Niedrigster EW in eigvals[0], zugehoeriger EV in eigvecs[:, 0]
    return float(eigvals[0]), eigvecs[:, 0]


def phi_from_coefs(coefs, L, sector='cos', n_grid=800):
    """Rekonstruiere phi(x) auf [-L, L] aus Fourier-Koeffizienten."""
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


def autocorr_from_phi(xs, phi):
    """Autokorrelation (phi * phi)(t) via FFT."""
    # Numerische Faltung ueber die Gitter-Darstellung
    n = len(phi)
    dx = xs[1] - xs[0]
    # (phi*phi)(t) = int phi(y) phi(y-t) dy, fuer echtes phi und Support in [-L, L]:
    #  Ergebnis hat Support [-2L, 2L]
    ac = np.convolve(phi, phi[::-1], mode='full') * dx
    ts = np.linspace(-2 * (xs[-1] - xs[0]) / 2, 2 * (xs[-1] - xs[0]) / 2, len(ac))
    # Einfacher: ts als Mapping
    ts = np.linspace(xs[0] - xs[-1], xs[-1] - xs[0], len(ac))
    return ts, ac


def S_chi_empirical(chi_fn, lam, phi_ac_fn, max_m=10):
    """S_chi mit empirischem Autokorrelations-Fenster (callable phi_ac_fn(t))."""
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
            S += weight * phi_ac_fn(t)
    return S


# --------------------------------------------------------------------------
# Baseline: Raw chi-Summe ohne Fenster
# --------------------------------------------------------------------------
def S_chi_raw(chi_fn, lam):
    """Rohe Summe sum_{p<=lam} chi(p) log p / sqrt p (keine Fenster-Gewichtung)."""
    S = 0.0
    for p in sympy.primerange(2, int(lam) + 1):
        cp = chi_fn(p)
        if cp == 0:
            continue
        S += cp * math.log(p) / math.sqrt(p)
    return S


# --------------------------------------------------------------------------
# Analyse-Funktion
# --------------------------------------------------------------------------
def compare(label, predictions, empirics, chi_names):
    sign_match = np.sign(predictions) == np.sign(empirics)
    sign_ok = int(np.sum(sign_match))
    if np.std(predictions) > 1e-12 and np.std(empirics) > 1e-12:
        r = float(np.corrcoef(predictions, empirics)[0, 1])
        r2 = r * r
    else:
        r = r2 = float('nan')
    print(f"\n[{label}] sign_ok = {sign_ok}/{len(empirics)},  R = {r:+.4f},  R^2 = {r2:.4f}")
    for name, S, g in zip(chi_names, predictions, empirics):
        ok = "✓" if np.sign(S) == np.sign(g) else "✗"
        print(f"  {name:>8}: pred = {S:+12.4e},  emp = {g:+8.5f}   {ok}")
    return {"label": label, "sign_ok": sign_ok, "total": len(empirics),
            "r": r, "r2": r2,
            "predictions": [float(v) for v in predictions],
            "empirics": [float(v) for v in empirics],
            "chi_names": chi_names}


def main():
    print("=" * 75)
    print("VALIDATION v2 — Multiple Windows + Empirical Ground State")
    print("=" * 75)

    LAM = 20000
    L = math.log(LAM)
    empirical = load_empirical_gaps()

    # Filter Charaktere, fuer die wir empirische Werte bei LAM haben
    chi_names = []
    emp_gaps = []
    for name, D in CHARS:
        if LAM in empirical[name]:
            g, N = empirical[name][LAM]
            chi_names.append(name)
            emp_gaps.append(g)

    emp_arr = np.array(emp_gaps)
    print(f"\n[N Charaktere]: {len(chi_names)}")
    for name, g in zip(chi_names, emp_gaps):
        print(f"  {name:>8}: emp gap = {g:+8.5f}")

    results_all = []

    # --- Test 1: Raw chi-Summe (keine Fenster-Gewichtung) ---
    print("\n" + "=" * 75)
    print("TEST 1: Baseline — Raw chi-Summe")
    print("=" * 75)
    preds = np.array([S_chi_raw(make_chi_D(next(d for n, d in CHARS if n == name)), LAM) for name in chi_names])
    results_all.append(compare("Raw chi-Summe", preds, emp_arr, chi_names))

    # --- Test 2: Gauss-Fenster um t=0 (v1-Bestfall, sigma=0.5) ---
    print("\n" + "=" * 75)
    print("TEST 2: Gauss um t=0 (v1 best: sigma=0.5)")
    print("=" * 75)
    preds = np.array([S_chi(make_chi_D(next(d for n, d in CHARS if n == name)), next(d for n, d in CHARS if n == name), LAM, 0.5) for name in chi_names])
    results_all.append(compare("Gauss t=0, sigma=0.5", preds, emp_arr, chi_names))

    # --- Test 3: Gauss-Fenster um t=L (Frontier-Skala) ---
    print("\n" + "=" * 75)
    print(f"TEST 3: Gauss um t=L={L:.3f} (Frontier-Skala)")
    print("=" * 75)
    for sigma_f in [0.5, 1.0, 2.0, 3.0]:
        preds = np.array([
            S_chi_window(
                make_chi_D(next(d for n, d in CHARS if n == name)),
                LAM,
                lambda t, s=sigma_f: gauss_offset(t, s, L),
            ) for name in chi_names
        ])
        results_all.append(compare(f"Gauss t=L, sigma={sigma_f}", preds, emp_arr, chi_names))

    # --- Test 4: Rechteck-Fenster [0, T] ---
    print("\n" + "=" * 75)
    print("TEST 4: Rechteck [0, T] fuer verschiedene T")
    print("=" * 75)
    for T in [L/4, L/2, L, 3*L/2, 2*L]:
        preds = np.array([
            S_chi_window(
                make_chi_D(next(d for n, d in CHARS if n == name)),
                LAM,
                lambda t, TT=T: rect_window(t, TT),
            ) for name in chi_names
        ])
        results_all.append(compare(f"Rect [0, {T:.2f}]", preds, emp_arr, chi_names))

    # --- Test 5: Cosine-Window (prolate-ish) ---
    print("\n" + "=" * 75)
    print("TEST 5: Cosine^2 Fenster auf [-T, T]")
    print("=" * 75)
    for T in [L/2, L, 3*L/2, 2*L]:
        preds = np.array([
            S_chi_window(
                make_chi_D(next(d for n, d in CHARS if n == name)),
                LAM,
                lambda t, TT=T: cos_window(t, TT),
            ) for name in chi_names
        ])
        results_all.append(compare(f"Cos^2 auf [-{T:.2f}, {T:.2f}]", preds, emp_arr, chi_names))

    # --- Test 6: Empirischer Grundzustand ---
    print("\n" + "=" * 75)
    print("TEST 6: Empirischer Grundzustand (Galerkin N=200)")
    print("=" * 75)
    N_gal = 200  # Reduziert auf 200 fuer Laufzeit
    emp_preds = []
    for name in chi_names:
        D = next(d for n, d in CHARS if n == name)
        chi_fn = make_chi_D(D)
        try:
            # Grundzustand im Cos-Sektor
            mu_cos, vec_cos = compute_ground_state(chi_fn, LAM, N=N_gal, sector='cos')
            xs, phi = phi_from_coefs(vec_cos, L, sector='cos')
            ts, ac = autocorr_from_phi(xs, phi)
            # Interpolator
            from scipy.interpolate import interp1d
            ac_fn = interp1d(ts, ac, bounds_error=False, fill_value=0.0)
            S = S_chi_empirical(chi_fn, LAM, ac_fn)
            emp_preds.append(S)
            print(f"  {name:>8}: mu_cos = {mu_cos:+8.4f}, S = {S:+8.4e}")
        except Exception as e:
            print(f"  {name:>8}: ERROR {e}")
            emp_preds.append(0.0)
    results_all.append(compare("Empirical cos-GS", np.array(emp_preds), emp_arr, chi_names))

    # --- Gesamt-Bestimmung ---
    print("\n" + "=" * 75)
    print("ZUSAMMENFASSUNG: Alle Tests")
    print("=" * 75)
    print(f"  {'Label':<40}  {'sign_ok':>8}  {'R^2':>8}")
    for r in results_all:
        print(f"  {r['label']:<40}  {r['sign_ok']:>3}/{r['total']:<3}     {r['r2']:>6.3f}")

    best = max(results_all, key=lambda r: (r["sign_ok"], r["r2"] if not math.isnan(r["r2"]) else -1))
    print(f"\n[BEST] {best['label']}: sign_ok = {best['sign_ok']}/{best['total']}, R^2 = {best['r2']:.4f}")

    # Verdikt
    if best["sign_ok"] >= 9 and best["r2"] >= 0.8:
        verdict = "FALL A: ERFOLGREICH"
    elif best["sign_ok"] >= 7:
        verdict = "FALL B: GEMISCHT (7-8 Vorzeichen)"
    else:
        verdict = "FALL C: SCHEITERT (Plan B aktivieren)"
    print(f"VERDIKT: {verdict}")

    # Speichern
    out = {"lam": LAM, "L": L, "chi_names": chi_names, "emp_gaps": list(emp_arr),
           "all_tests": results_all, "best": best, "verdict": verdict}
    out_json = RES / "ANALYTIC_GAP_FORMULA_TEST_V2.json"
    with open(out_json, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] {out_json}")


if __name__ == "__main__":
    main()
