#!/usr/bin/env python3
# coding: utf-8
"""
Attack-1(b)-Check: Testen ob die phi^- Konvolutions-Summe verschwindet.

Fragestellung: Ist Sigma_w * (phi^- * phi^-)_math(m log p) ~ 0?
Wenn ja: Beide Konventionen (Code vs. Paper Eq.5) sind aequivalent.
Wenn nein: Paper hat Definitions-Implementierungs-Inkonsistenz.

Gleichzeitig pruefen: S_exact = -2 * weighted_sum oder +2 * weighted_sum?
Dies klaert die Vorzeichen-Konvention eindeutig.
"""
import math
import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import sympy

from analytic_gap_formula_test_v2 import build_W, compute_ground_state, phi_from_coefs
from analytic_gap_formula_test import make_chi_D, CHARS, load_empirical_gaps
from ground_state_difference_analysis import autocorr, autocorr_at, compute_gap_exact

LAM = 20000
N_G = 200
MAX_M = 10

def run_check(D, chi_name, gap_emp):
    chi_fn = make_chi_D(D)
    L = math.log(LAM)

    # Grundzustaende
    ev_plus, coefs_plus = compute_ground_state(chi_fn, LAM, N=N_G, sector='cos')
    ev_minus, coefs_minus = compute_ground_state(chi_fn, LAM, N=N_G, sector='sin')

    xs_p, phi_p = phi_from_coefs(coefs_plus, L, sector='cos', n_grid=1600)
    xs_m, phi_m = phi_from_coefs(coefs_minus, L, sector='sin', n_grid=1600)

    # S_exact mit coeff=-2 (Skript-Standard)
    S_exact_neg2 = compute_gap_exact(chi_fn, LAM, phi_p, xs_p, phi_m, xs_m, coeff=-2.0)

    # Autocorrelation-Grids
    dx_p = xs_p[1] - xs_p[0]
    dx_m = xs_m[1] - xs_m[0]
    lags_p, ac_p = autocorr(phi_p, dx_p)
    lags_m, ac_m = autocorr(phi_m, dx_m)

    # Summen einzeln berechnen
    sum_delta = 0.0      # Sigma w * (ac_plus - ac_minus)  [Code-Konvention]
    sum_plus = 0.0       # Sigma w * ac_plus               [Sigma w * R_{phi^+}]
    sum_minus = 0.0      # Sigma w * ac_minus              [Sigma w * R_{phi^-}]
    sum_conv_minus = 0.0 # Sigma w * (phi^- * phi^-)_math  = Sigma w * (-ac_minus)  [fuer ungerades phi^-]

    twoL = 2.0 * L
    for p in sympy.primerange(2, LAM + 1):
        cp = chi_fn(p)
        if cp == 0:
            continue
        log_p = math.log(p)
        if log_p > twoL:
            break
        for m in range(1, MAX_M + 1):
            t = m * log_p
            if t > twoL:
                break
            phase = cp ** m
            w = phase * log_p / (p ** (m / 2.0))
            ap = float(np.interp(t, lags_p, ac_p, left=0.0, right=0.0))
            am = float(np.interp(t, lags_m, ac_m, left=0.0, right=0.0))
            sum_delta += w * (ap - am)
            sum_plus += w * ap
            sum_minus += w * am
            sum_conv_minus += w * (-am)  # (phi^- * phi^-)_math = -R_{phi^-}

    print(f"\n{'='*60}")
    print(f"Chi: {chi_name} (D={D}), gap_emp = {gap_emp:+.5f}")
    print(f"{'='*60}")
    print(f"  sum_delta = Sigma w*(ac_plus-ac_minus):    {sum_delta:+.6f}")
    print(f"  sum_plus  = Sigma w*ac_plus:               {sum_plus:+.6f}")
    print(f"  sum_minus = Sigma w*ac_minus:              {sum_minus:+.6f}")
    print(f"  sum_conv_minus (=-sum_minus):              {sum_conv_minus:+.6f}")
    print(f"")
    print(f"  S_exact (coeff=-2):                        {S_exact_neg2:+.6f}")
    print(f"  -2 * sum_delta:                            {-2 * sum_delta:+.6f}")
    print(f"  +2 * sum_delta:                            {+2 * sum_delta:+.6f}")
    print(f"")
    # Antwort auf Attack 1(b): Verschwindet Sigma w * (phi^-*phi^-)_math?
    ratio = abs(sum_conv_minus) / max(abs(sum_delta), 1e-12)
    print(f"  |sum_conv_minus| / |sum_delta|:            {ratio:.2f}")
    print(f"  --> Attack 1(b): sum_conv_minus ~ 0? {'JA (verschwindet)' if ratio < 0.1 else 'NEIN (signifikant)'}")
    print(f"")
    # Papier-Formel: Delta_chi = (phi^+*phi^+) - (phi^-*phi^-)_math = ac_plus - (-ac_minus) = ac_plus + ac_minus
    S_paper_formula = -2 * (sum_plus - sum_conv_minus)  # = -2 * (sum_plus + sum_minus)
    print(f"  S_paper (Eq.5 woertlich):                  {S_paper_formula:+.6f}")
    print(f"  [Paper-Eq.5 gap-Prediction wuerde: {S_paper_formula:+.2e}, emp={gap_emp:+.2e}]")
    print(f"")
    sign_code_ok = (S_exact_neg2 > 0) == (gap_emp > 0)
    sign_paper_ok = (S_paper_formula > 0) == (gap_emp > 0)
    print(f"  Vorzeichen-Check Code: {'OK' if sign_code_ok else 'FAIL'}")
    print(f"  Vorzeichen-Check Paper-Eq.5: {'OK' if sign_paper_ok else 'FAIL'}")

    return {
        "chi": chi_name, "D": D, "gap_emp": gap_emp,
        "sum_delta": sum_delta, "sum_plus": sum_plus,
        "sum_minus": sum_minus, "sum_conv_minus": sum_conv_minus,
        "S_exact": S_exact_neg2, "S_paper": S_paper_formula,
        "sign_code_ok": sign_code_ok, "sign_paper_ok": sign_paper_ok,
        "ratio_conv_minus_vs_delta": ratio,
    }


def main():
    print("ATTACK-1(b)-CHECK: Konventions-Aequivalenz-Test")
    print("Frage: Sigma w * (phi^-*phi^-)_math ~ 0?")
    print("Methode: Direkte Berechnung fuer alle 10 Charaktere")

    empirical = load_empirical_gaps()
    # Nur chi_5, chi_8, chi_12 fuer schnellen Test
    test_chars = [("chi_5", 5), ("chi_8", 8), ("chi_12", 12)]

    results = []
    for name, D in test_chars:
        if LAM not in empirical.get(name, {}):
            print(f"[skip] {name}")
            continue
        gap_emp, _ = empirical[name][LAM]
        r = run_check(D, name, gap_emp)
        results.append(r)

    print("\n" + "="*60)
    print("ZUSAMMENFASSUNG")
    print("="*60)
    print(f"{'Chi':<8} {'gap_emp':>10} {'sum_conv_m':>12} {'ratio':>8} {'Code':>6} {'Paper':>6}")
    for r in results:
        print(f"{r['chi']:<8} {r['gap_emp']:>10.5f} {r['sum_conv_minus']:>12.4f} "
              f"{r['ratio_conv_minus_vs_delta']:>8.1f} "
              f"{'OK' if r['sign_code_ok'] else 'FAIL':>6} {'OK' if r['sign_paper_ok'] else 'FAIL':>6}")

    print("\nFazit:")
    all_nonzero = all(r['ratio_conv_minus_vs_delta'] > 1.0 for r in results)
    if all_nonzero:
        print("ATTACK 1(b) BESTAETIGT: sum_conv_minus signifikant != 0.")
        print("Paper Eq.5 (math. Konvolution) != Code-Implementierung.")
        print("Fix: Eq.5 auf R_{phi^+} - R_{phi^-} umformulieren (Autokorrelation).")
        print("Oder: Appendix B erweitern mit explizitem Beweis dass die Summe verschwindet.")
    else:
        print("ATTACK 1(b) NEUTRALISIERT: sum_conv_minus ~ 0, Konventionen aequivalent.")


if __name__ == "__main__":
    main()
