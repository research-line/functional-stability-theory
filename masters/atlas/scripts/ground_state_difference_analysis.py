#!/usr/bin/env python3
# coding: utf-8
"""
Session 7 Woche 1: Ground-State-Differenz-Analyse.

Direkter Test der EXAKTEN Weil-Kern-Gap-Formel aus ANALYTIC_KERNEL.md Thm 4.1
OHNE die in Session 6 falsifizierte Reduktion phi^+ ~ phi^-:

    gap_chi(lambda) ~ 2 * sum_{p,m} chi(p)^m * log p / p^(m/2) * Delta_chi(m * log p)

mit

    Delta_chi(t) := (phi^+ * phi^+)(t) - (phi^- * phi^-)(t)

wobei phi^+ der Grundzustand des cos-Sektors und phi^- der Grundzustand des sin-Sektors
sind (beide numerisch via Galerkin-Diagonalisierung, N=200).

Zwei Erfolgskriterien:
  (a) Full-formula Test: sign_ok >= 9/10, R^2 >= 0.8 -> Option B validiert
  (b) Invarianten-Scan: welche skalaren Kenngroessen von Delta_chi / phi^+ / phi^-
      korrelieren mit dem empirischen Gap?

Ausgaben:
  _results/GROUND_STATE_DIFFERENCE_ANALYSIS.json   (raw data)
  _results/GROUND_STATE_DIFFERENCE_ANALYSIS.md     (summary report)
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import sympy

sys.stdout.reconfigure(line_buffering=True)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "_results"

# Reuse Galerkin-Infrastruktur aus v2
from analytic_gap_formula_test_v2 import (
    build_W, compute_ground_state, phi_from_coefs,
)
from analytic_gap_formula_test import (
    make_chi_D, CHARS, load_empirical_gaps,
)


# --------------------------------------------------------------------------
# Autokorrelations-Primitive
# --------------------------------------------------------------------------
def autocorr(phi, dx):
    """Diskrete Autokorrelation (phi * phi)(t) bei lag t = (k - (N-1)) * dx.

    Liefert ac[k] fuer k=0..2N-2 und entsprechende lags.
    """
    ac = np.correlate(phi, phi, mode='full') * dx
    N = len(phi)
    lags = np.arange(-(N - 1), N) * dx
    return lags, ac


def autocorr_at(phi, xs, ts):
    """Autokorrelation phi ausgewertet an gegebenen Lag-Werten ts (interpoliert)."""
    dx = xs[1] - xs[0]
    lags, ac = autocorr(phi, dx)
    return np.interp(ts, lags, ac, left=0.0, right=0.0)


# --------------------------------------------------------------------------
# Exakte Gap-Formel (ohne phi^+ ~ phi^- Reduktion)
# --------------------------------------------------------------------------
def compute_gap_exact(chi_fn, lam, phi_plus, xs_plus, phi_minus, xs_minus,
                     max_m=10, coeff=-2.0):
    """S_chi = coeff * sum_{p,m} chi(p)^m * log p / p^(m/2) * Delta_chi(m * log p)
    mit Delta_chi(t) = (phi^+ * phi^+)(t) - (phi^- * phi^-)(t).

    Vorzeichen-Konvention (validiert 2026-04-17 Session 7):
    Mit build_W-Konvention (Prim-Term +w in der Matrix) ist
      <phi^pm, K^pm phi^pm> = arch + 2 * sum w_p * (phi^pm * phi^pm)(m log p)
    Daher gap = ev_minus - ev_plus = 2 * sum w_p * (AC_minus - AC_plus)
              = -2 * sum w_p * Delta_chi
    Also coeff = -2. Ein initialer coeff=+2 gab anti-korrelierte Predictions
    (R=-0.70, 2/10); Flip auf -2 liefert R=+0.70, 8/10.
    """
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
            # Delta_chi(t):
            ac_plus = float(autocorr_at(phi_plus, xs_plus, np.array([t]))[0])
            ac_minus = float(autocorr_at(phi_minus, xs_minus, np.array([t]))[0])
            delta_t = ac_plus - ac_minus
            S += weight * delta_t
    return coeff * S


# --------------------------------------------------------------------------
# Invarianten von phi^+, phi^-, Delta_chi
# --------------------------------------------------------------------------
def invariants(chi_name, chi_fn, lam, phi_plus_coefs, phi_minus_coefs, N_Galerkin=200):
    """Extrahiere skalar-Invarianten von den Grundzustaenden fuer Korrelationsanalyse."""
    L = math.log(lam)
    xs_p, phi_p = phi_from_coefs(phi_plus_coefs, L, sector='cos', n_grid=1600)
    xs_m, phi_m = phi_from_coefs(phi_minus_coefs, L, sector='sin', n_grid=1600)

    # Autokorrelationen
    dx_p = xs_p[1] - xs_p[0]
    dx_m = xs_m[1] - xs_m[0]
    lags_p, ac_p = autocorr(phi_p, dx_p)
    lags_m, ac_m = autocorr(phi_m, dx_m)

    # Werte an Prim-Log-Punkten
    max_m = 6
    prime_ts = []
    delta_at_primes = []
    weighted_sum = 0.0
    weighted_sum_p = 0.0
    weighted_sum_m = 0.0
    for p in sympy.primerange(2, int(lam) + 1):
        cp = chi_fn(p)
        if cp == 0:
            continue
        log_p = math.log(p)
        if log_p > 2 * L:
            break
        for me in range(1, max_m + 1):
            t = me * log_p
            if t > 2 * L:
                break
            prime_ts.append(t)
            a_p = float(np.interp(t, lags_p, ac_p, left=0.0, right=0.0))
            a_m = float(np.interp(t, lags_m, ac_m, left=0.0, right=0.0))
            d = a_p - a_m
            delta_at_primes.append(d)
            w = (cp ** me) * log_p / (p ** (me / 2.0))
            weighted_sum += w * d
            weighted_sum_p += w * a_p
            weighted_sum_m += w * a_m

    # Skalare Invarianten
    inv = {
        "chi": chi_name,
        "lambda": lam,
        "L": L,
        "N_Galerkin": N_Galerkin,
        # Fourier-Schwerpunkt (naeherungsweise: ersten 20 Koeffizienten)
        "phi_plus_centroid_index": float(np.sum(np.arange(len(phi_plus_coefs)) * phi_plus_coefs**2)
                                          / max(np.sum(phi_plus_coefs**2), 1e-30)),
        "phi_minus_centroid_index": float(np.sum(np.arange(len(phi_minus_coefs)) * phi_minus_coefs**2)
                                           / max(np.sum(phi_minus_coefs**2), 1e-30)),
        "phi_plus_L2": float(np.sum(phi_plus_coefs**2)),
        "phi_minus_L2": float(np.sum(phi_minus_coefs**2)),
        "ac_plus_at_0": float(np.interp(0.0, lags_p, ac_p)),
        "ac_minus_at_0": float(np.interp(0.0, lags_m, ac_m)),
        "delta_at_0": float(np.interp(0.0, lags_p, ac_p) - np.interp(0.0, lags_m, ac_m)),
        "delta_L1_prime_norm": float(np.sum(np.abs(delta_at_primes))),
        "delta_L2_prime_norm": float(np.sqrt(np.sum(np.array(delta_at_primes)**2))),
        "delta_signed_sum_primes": float(np.sum(delta_at_primes)),
        "weighted_sum": float(weighted_sum),
        "weighted_sum_plus": float(weighted_sum_p),
        "weighted_sum_minus": float(weighted_sum_m),
        "n_prime_terms": len(prime_ts),
        # Erste 10 prime-t Werte fuer Inspection
        "delta_first10": [float(d) for d in delta_at_primes[:10]],
        "prime_ts_first10": [float(t) for t in prime_ts[:10]],
    }
    return inv


# --------------------------------------------------------------------------
# Haupt-Analyse
# --------------------------------------------------------------------------
def analyze_all(lam_target=20000, N_Galerkin=200, verbose=True):
    """Fuer alle 10 Charaktere: Ground-States + exakte Formel + Invarianten."""
    empirical = load_empirical_gaps()
    results = []

    print(f"\n[analyze] lam = {lam_target}, L = {math.log(lam_target):.4f}, N = {N_Galerkin}")
    print(f"[analyze] 10 Charaktere, je 2 Sektoren -> {20} Galerkin-Diagonalisierungen")

    for name, D in CHARS:
        if lam_target not in empirical.get(name, {}):
            print(f"[skip] {name}: keine empirischen Daten bei lam={lam_target}")
            continue
        gap_emp, N_emp = empirical[name][lam_target]
        chi_fn = make_chi_D(D)
        L = math.log(lam_target)

        t0 = time.time()
        # Grundzustaende: cos (phi^+) und sin (phi^-)
        ev_plus, coefs_plus = compute_ground_state(chi_fn, lam_target, N=N_Galerkin, sector='cos')
        ev_minus, coefs_minus = compute_ground_state(chi_fn, lam_target, N=N_Galerkin, sector='sin')

        # Rekonstruktion
        xs_p, phi_p = phi_from_coefs(coefs_plus, L, sector='cos', n_grid=1600)
        xs_m, phi_m = phi_from_coefs(coefs_minus, L, sector='sin', n_grid=1600)

        # Exakte Formel mit korrektem Vorzeichen (coeff=-2)
        S_exact = compute_gap_exact(chi_fn, lam_target, phi_p, xs_p, phi_m, xs_m, coeff=-2.0)

        # Invarianten
        inv = invariants(name, chi_fn, lam_target, coefs_plus, coefs_minus, N_Galerkin)

        dt = time.time() - t0

        sign_ok = np.sign(S_exact) == np.sign(gap_emp)
        print(f"[{name:>8}] D={D:>2}, gap_emp={gap_emp:+8.5f} (N={N_emp}),  "
              f"ev+={ev_plus:+.4f}, ev-={ev_minus:+.4f},  "
              f"S_exact={S_exact:+9.5e}   {'OK' if sign_ok else 'FAIL'}   [{dt:.1f}s]")

        results.append({
            "chi": name,
            "D": D,
            "lambda": lam_target,
            "gap_emp": float(gap_emp),
            "N_emp": int(N_emp),
            "eigenvalue_plus": float(ev_plus),
            "eigenvalue_minus": float(ev_minus),
            "S_exact": float(S_exact),
            "gap_ev_diff": float(ev_minus - ev_plus),
            "sign_ok": bool(sign_ok),
            **{k: v for k, v in inv.items() if k not in ("chi", "lambda")},
        })

    return results


def correlation_stats(results, field, emp_field="gap_emp"):
    """Korrelation zwischen results[field] und results[emp_field]."""
    x = np.array([r[field] for r in results])
    y = np.array([r[emp_field] for r in results])
    if np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return {"r": float("nan"), "r2": float("nan"), "sign_ok": 0}
    r = float(np.corrcoef(x, y)[0, 1])
    sign_ok = int(np.sum(np.sign(x) == np.sign(y)))
    return {"r": r, "r2": r * r, "sign_ok": sign_ok, "n": len(x)}


def main():
    print("=" * 75)
    print("GROUND-STATE-DIFFERENCE ANALYSIS (Session 7 Woche 1)")
    print("=" * 75)
    print("Direkter Test der exakten Weil-Kern-Gap-Formel (ANALYTIC_KERNEL Thm 4.1)")
    print("OHNE Reduktion phi^+ ~ phi^- (die in Session 6 falsifiziert wurde).")
    print()

    LAM = 20000
    N_G = 200
    print(f"Parameter: lambda = {LAM}, N_Galerkin = {N_G}")
    print()

    results = analyze_all(lam_target=LAM, N_Galerkin=N_G, verbose=True)

    print("\n" + "=" * 75)
    print("KORRELATIONS-ANALYSE (Invarianten vs. empirischer Gap)")
    print("=" * 75)

    fields_to_test = [
        "S_exact",             # Exakte Gap-Formel (Hauptkandidat)
        "gap_ev_diff",         # ev(sin) - ev(cos), klassischer Gap-Proxy
        "weighted_sum",        # Prim-gewichtete Delta-Summe, = S_exact / 2
        "weighted_sum_plus",   # Nur cos-Sektor-Beitrag
        "weighted_sum_minus",  # Nur sin-Sektor-Beitrag
        "delta_signed_sum_primes",
        "delta_L1_prime_norm",
        "delta_L2_prime_norm",
        "delta_at_0",
        "phi_plus_centroid_index",
        "phi_minus_centroid_index",
        "eigenvalue_plus",
        "eigenvalue_minus",
    ]
    corr = {}
    for field in fields_to_test:
        stats = correlation_stats(results, field)
        corr[field] = stats
        print(f"  {field:>28}:  R = {stats['r']:+.4f},  R^2 = {stats['r2']:.4f},  "
              f"sign_ok = {stats['sign_ok']:>2}/{stats.get('n', len(results))}")

    # Save
    RES.mkdir(exist_ok=True)
    out_json = RES / "GROUND_STATE_DIFFERENCE_ANALYSIS.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump({
            "parameters": {"lambda": LAM, "N_Galerkin": N_G},
            "results": results,
            "correlations": corr,
        }, f, indent=2)
    print(f"\n[out] {out_json}")

    # Markdown Summary
    main_field = "S_exact"
    s = corr[main_field]
    verdict = "SUCCESS" if (s["sign_ok"] >= 9 and s["r2"] >= 0.8) else "FAIL"

    md_lines = [
        "# Ground-State-Differenz-Analyse — Session 7 Woche 1",
        "",
        f"**Datum:** {time.strftime('%Y-%m-%d')}",
        f"**Parameter:** lambda = {LAM}, N_Galerkin = {N_G}",
        "**Formel:** gap_chi ~ 2 * sum_{p,m} chi(p)^m * log p / p^(m/2) * [(phi^+*phi^+)(m log p) - (phi^-*phi^-)(m log p)]",
        "",
        "## Hauptresultat",
        "",
        f"**S_exact vs. gap_emp:** sign_ok = {s['sign_ok']}/10, R^2 = {s['r2']:.4f} ({verdict})",
        "",
        "## Tabelle",
        "",
        "| chi | D | gap_emp | S_exact | sign | ev+ | ev- | ev_diff |",
        "|-----|---|---------|---------|------|-----|-----|---------|",
    ]
    for r in results:
        ok = "OK" if r["sign_ok"] else "FAIL"
        md_lines.append(
            f"| {r['chi']} | {r['D']} | {r['gap_emp']:+.5f} | {r['S_exact']:+.4e} | {ok} | "
            f"{r['eigenvalue_plus']:+.4f} | {r['eigenvalue_minus']:+.4f} | {r['gap_ev_diff']:+.4f} |"
        )
    md_lines += ["", "## Korrelations-Scan aller Invarianten", "",
                 "| Invariante | R | R^2 | sign_ok |",
                 "|------------|---|-----|---------|"]
    for k, c in corr.items():
        md_lines.append(f"| {k} | {c['r']:+.4f} | {c['r2']:.4f} | {c['sign_ok']} |")

    out_md = RES / "GROUND_STATE_DIFFERENCE_ANALYSIS.md"
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines) + "\n")
    print(f"[out] {out_md}")


if __name__ == "__main__":
    main()
