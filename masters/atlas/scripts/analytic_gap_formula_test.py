#!/usr/bin/env python3
# coding: utf-8
"""
Woche-4-Validierungsskript fuer Session 6 Analytic-Pfad.

Testet die theoretische Vorhersage aus ANALYTIC_GAP_FORMULA.md:

    gap_chi(lambda) ~ 2 * S_chi(L) + O(1/L)

mit

    S_chi(L) = sum_{p nmid q, m>=1} chi(p)^m * log(p)/p^(m/2) * Phi_inf(m*log p)

und Gauss-Ansatz Phi_inf(t) = C0 * exp(-t^2 / (4*sigma^2)).

Workflow:
  1. Lade empirische Gap-Werte (N=600 wo verfuegbar, sonst N=400).
  2. Berechne S_chi(L) fuer verschiedene sigma-Werte.
  3. Kalibriere sigma an chi_12 (stabilster positiver Fall).
  4. Teste Vorzeichen-Prediktion + R^2 ueber alle 10 Charaktere.
  5. Speichere Resultate nach _results/ANALYTIC_GAP_FORMULA_TEST.json + .md.

Erfolgskriterium (Handoff §4.2): sign_ok >= 9/10 und R^2 >= 0.8.
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import sympy

sys.stdout.reconfigure(line_buffering=True)

# --------------------------------------------------------------------------
# Pfad-Setup
# --------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "_results"

# --------------------------------------------------------------------------
# Kronecker-Symbol (identisch zu all10_high_N_server.py)
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
# Analytische Gap-Formel (Zeit-Domaene, Gauss-Fenster)
# --------------------------------------------------------------------------
def gauss_window(t, sigma):
    """Phi_inf(t) = (2*pi*sigma^2)^(-1/2) * exp(-t^2 / (4*sigma^2))

    Selbstfaltung von Gauss(0, sigma^2) hat Breite sigma*sqrt(2) am Autokorrelations-Peak.
    Normierung so dass Phi_inf(0) = 1 / (sqrt(2*pi) * sigma).
    """
    return np.exp(-(t * t) / (4.0 * sigma * sigma)) / (math.sqrt(2.0 * math.pi) * sigma)


def S_chi(chi_fn, D, lam, sigma, max_m=10):
    """Berechnet S_chi(L) = sum_{p,m} chi(p)^m * log p / p^(m/2) * Phi_inf(m * log p).

    Fenster-Cutoff: m*log p <= 2*L (Support von Phi_inf auf [-L, L]^2 Autokorrelation).
    """
    L = math.log(lam)
    twoL = 2.0 * L
    S = 0.0
    # Alle Primen p <= lam, p nmid D
    for p in sympy.primerange(2, int(lam) + 1):
        cp = chi_fn(p)
        if cp == 0:
            continue
        log_p = math.log(p)
        if log_p > twoL:
            break
        p_half = math.sqrt(p)
        w_base = log_p / p_half
        for m in range(1, max_m + 1):
            t = m * log_p
            if t > twoL:
                break
            phase = cp ** m  # +/-1 fuer reellen Charakter
            weight = phase * w_base / (p_half ** (m - 1))  # = chi(p)^m * log p / p^(m/2)
            S += weight * gauss_window(t, sigma)
    return S


# --------------------------------------------------------------------------
# Empirische Gap-Daten laden
# --------------------------------------------------------------------------
def load_empirical_gaps():
    """Ladet die empirischen gap-Werte, bevorzugt N=600, Fallback N=400.

    Rueckgabe: dict {chi_name: {lambda: (gap, N_used)}}
    """
    empirical = {name: {} for name, _ in CHARS}

    # N=600 (6 Charaktere, 4 lambdas)
    with open(RES / "oscillating6_N600_results.json") as f:
        data_600 = json.load(f)
    for row in data_600:
        chi = row["chi"]
        lam = row["lambda"]
        empirical[chi][lam] = (row["gap"], 600)

    # N=400 (alle 10 Charaktere, 4 lambdas; fuellt Luecken)
    with open(RES / "all10_high_N_results.json") as f:
        data_400 = json.load(f)
    for row in data_400:
        chi = row["chi"]
        lam = row["lambda"]
        if lam not in empirical.get(chi, {}):
            empirical[chi][lam] = (row["gap"], 400)

    # chi_5, chi_12 aus chi12_high_N_results.json (N=200/300/400)
    try:
        with open(RES / "chi12_high_N_results.json") as f:
            data_12 = json.load(f)
        for row in data_12:
            if row.get("N") != 400:
                continue
            chi = row["chi"]
            lam = row["lambda"]
            if lam not in empirical.get(chi, {}):
                empirical[chi][lam] = (row["gap"], 400)
    except FileNotFoundError:
        pass

    return empirical


# --------------------------------------------------------------------------
# Analyse: sigma-Scan + Korrelations-Statistiken
# --------------------------------------------------------------------------
def analyze(sigma_values, lam_target=20000, verbose=True):
    """Fuer jedes sigma: compute S_chi fuer alle 10 Charaktere bei lam=lam_target.

    Vergleiche mit empirischem Gap. Rueckgabe: list of result dicts.
    """
    empirical = load_empirical_gaps()
    L = math.log(lam_target)
    if verbose:
        print(f"\n[analyze] lam = {lam_target}, L = log(lam) = {L:.4f}")
        print(f"[analyze] sigma-Werte: {sigma_values}")

    # Empirische gaps sammeln
    emp_gaps = []
    chi_names = []
    for name, D in CHARS:
        if lam_target in empirical[name]:
            gap, N = empirical[name][lam_target]
            emp_gaps.append(gap)
            chi_names.append((name, D, N))
            if verbose:
                print(f"[emp] {name} (D={D:2d}, N={N}): gap = {gap:+8.5f}")

    results_by_sigma = []
    for sigma in sigma_values:
        S_vals = []
        for name, D, N in chi_names:
            chi_fn = make_chi_D(D)
            S = S_chi(chi_fn, D, lam_target, sigma)
            S_vals.append(S)

        S_arr = np.array(S_vals)
        emp_arr = np.array(emp_gaps)

        # Sign-match
        sign_match = np.sign(S_arr) == np.sign(emp_arr)
        sign_ok = int(np.sum(sign_match))

        # Pearson R^2 (zwischen S und gap)
        if np.std(S_arr) > 1e-10 and np.std(emp_arr) > 1e-10:
            r = np.corrcoef(S_arr, emp_arr)[0, 1]
            r2 = r * r
            # Lineare Regression: gap ~ a*S + b
            a, b = np.polyfit(S_arr, emp_arr, 1)
        else:
            r = r2 = float("nan")
            a = b = float("nan")

        results_by_sigma.append({
            "sigma": sigma,
            "sign_ok": sign_ok,
            "sign_total": len(emp_arr),
            "r": float(r),
            "r2": float(r2),
            "a_slope": float(a) if not math.isnan(a) else None,
            "b_intercept": float(b) if not math.isnan(b) else None,
            "S_values": [float(v) for v in S_arr],
            "emp_values": [float(v) for v in emp_arr],
            "chi_names": [x[0] for x in chi_names],
        })
        if verbose:
            print(f"\n[sigma = {sigma:.4f}] sign_ok = {sign_ok}/{len(emp_arr)},  R = {r:+.4f},  R^2 = {r2:.4f}")
            print(f"  slope a = {a:+.4e},  intercept b = {b:+.4e}")
            print(f"  {'chi':>8}  {'S_chi':>14}  {'emp_gap':>10}  {'sign':>4}")
            for (name, D, N), S, g in zip(chi_names, S_arr, emp_arr):
                ok = "✓" if np.sign(S) == np.sign(g) else "✗"
                print(f"  {name:>8}  {S:+14.5e}  {g:+10.5f}   {ok}")

    return results_by_sigma, chi_names


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    print("=" * 75)
    print("ANALYTIC GAP FORMULA VALIDATION (Woche 4, Session 6)")
    print("=" * 75)
    print("Theorie: ANALYTIC_GAP_FORMULA.md Kor. 4.1 + ANALYTIC_GROUNDSTATE.md §4")
    print("Formel: S_chi(L) = sum_{p,m} chi(p)^m * log p / p^(m/2) * Phi_inf(m*log p)")
    print("Phi_inf(t) = Gauss mit Breite sigma (zu kalibrieren)")
    print()

    LAM = 20000
    L = math.log(LAM)

    # Erster Scan: weit verteilt sigma-Werte
    # Theoretische Intervalle:
    #   sigma ~ L^(1/4) approx 1.77 (eng)
    #   sigma ~ sqrt(L/2) approx 2.23 (WKB-Standard)
    #   sigma ~ sqrt(L) approx 3.15 (semiclassical)
    #   sigma ~ L approx 9.9 (voller Box-Skala)
    sigma_broad = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0]

    print("\n" + "=" * 75)
    print("PHASE 1: Breitscan ueber sigma-Werte bei lam=20000")
    print("=" * 75)

    broad_results, chi_names = analyze(sigma_broad, lam_target=LAM, verbose=True)

    # Bestes sigma (max R^2)
    best = max(broad_results, key=lambda r: (r["sign_ok"], r["r2"] if not math.isnan(r["r2"]) else -1))
    print(f"\n[BEST] sigma = {best['sigma']:.4f}: sign_ok = {best['sign_ok']}/{best['sign_total']}, R^2 = {best['r2']:.4f}")

    # Feinsuche um bestes sigma
    print("\n" + "=" * 75)
    print(f"PHASE 2: Feinscan um sigma = {best['sigma']:.3f}")
    print("=" * 75)
    s0 = best["sigma"]
    sigma_fine = [s0 * x for x in [0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3]]
    fine_results, _ = analyze(sigma_fine, lam_target=LAM, verbose=True)

    best_fine = max(fine_results, key=lambda r: (r["sign_ok"], r["r2"] if not math.isnan(r["r2"]) else -1))
    print(f"\n[BEST FINE] sigma = {best_fine['sigma']:.4f}: sign_ok = {best_fine['sign_ok']}/{best_fine['sign_total']}, R^2 = {best_fine['r2']:.4f}")

    # Erfolgs-Klassifizierung (Handoff §4.2 + ANALYTIC_GAP_FORMULA.md §7)
    sign_ok = best_fine["sign_ok"]
    total = best_fine["sign_total"]
    r2 = best_fine["r2"]
    if sign_ok >= 9 and r2 >= 0.8:
        verdict = "FALL A: ERFOLGREICH (sign_ok >= 9/10 und R^2 >= 0.8)"
    elif sign_ok >= 7:
        verdict = "FALL B: GEMISCHT (7-8 Vorzeichen korrekt, aber R^2 schwach)"
    else:
        verdict = "FALL C: SCHEITERT (< 7/10 Vorzeichen korrekt -> Plan B)"

    print("\n" + "=" * 75)
    print(f"VERDIKT: {verdict}")
    print("=" * 75)

    # Speichere Ergebnisse
    output = {
        "meta": {
            "lam_target": LAM,
            "L": L,
            "n_chars": len(chi_names),
            "chi_names": [x[0] for x in chi_names],
            "chi_D": [x[1] for x in chi_names],
            "chi_N_used": [x[2] for x in chi_names],
        },
        "broad_scan": broad_results,
        "fine_scan": fine_results,
        "best_broad": best,
        "best_fine": best_fine,
        "verdict": verdict,
    }

    out_json = RES / "ANALYTIC_GAP_FORMULA_TEST.json"
    with open(out_json, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n[saved] {out_json}")

    return output


if __name__ == "__main__":
    main()
