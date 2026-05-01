#!/usr/bin/env python3
# coding: utf-8
"""
Session 7 Woche 1 Schritt (a): chi_21 N-Konvergenz-Analyse.

Validiert die N-Truncation-Hypothese aus Session 7 Teil 3: chi_21 zeigt bei
N=200 gap_gal=+0.28, bei N=600 gap_emp=-0.004 — das ist bekannte Oszillation.
Wir laufen nun chi_21 bei N ∈ {400, 600} und messen gap_galerkin + S_with_arch
an jedem N, um die Konvergenz explizit zu dokumentieren.

Nur chi_21 (Ziel) und chi_12 (Kontrolle, stabil bereits bei N=200) — mit
chi_12 als sanity check dass die gesamte Infrastruktur bei hoeherem N
dieselben stabilen Zahlen liefert.
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "_results"

from analytic_gap_formula_test_v2 import (
    build_W, compute_ground_state, phi_from_coefs,
)
from analytic_gap_formula_test import (
    make_chi_D, load_empirical_gaps,
)
from ground_state_difference_analysis import (
    autocorr, autocorr_at, compute_gap_exact,
)
from arch_term_analysis import (
    build_W_arch_only, arch_contrib,
)


# Nur zwei Charaktere
TARGETS = [
    ("chi_21", 21),   # Oszillations-Kandidat
    ("chi_12", 12),   # stabile Kontrolle
]


def analyze_convergence(chi_name, D, lam=20000, N_values=(200, 400, 600)):
    chi_fn = make_chi_D(D)
    L = math.log(lam)
    results = []

    for N in N_values:
        t0 = time.time()
        print(f"\n[{chi_name}] N={N} — starte Galerkin...")
        ev_plus, coefs_plus = compute_ground_state(chi_fn, lam, N=N, sector='cos')
        ev_minus, coefs_minus = compute_ground_state(chi_fn, lam, N=N, sector='sin')

        arch_plus = arch_contrib(coefs_plus, 'cos', L)
        arch_minus = arch_contrib(coefs_minus, 'sin', L)
        arch_diff = arch_minus - arch_plus

        xs_p, phi_p = phi_from_coefs(coefs_plus, L, sector='cos', n_grid=1600)
        xs_m, phi_m = phi_from_coefs(coefs_minus, L, sector='sin', n_grid=1600)
        S_exact = compute_gap_exact(chi_fn, lam, phi_p, xs_p, phi_m, xs_m, coeff=-2.0)

        gap_galerkin = ev_minus - ev_plus
        S_with_arch = S_exact + arch_diff

        dt = time.time() - t0
        print(f"  ev_plus = {ev_plus:+.5f}, ev_minus = {ev_minus:+.5f}")
        print(f"  gap_galerkin = {gap_galerkin:+.5f}")
        print(f"  S_exact = {S_exact:+.5e}, arch_diff = {arch_diff:+.5e}")
        print(f"  S_with_arch = {S_with_arch:+.5e}   [dt={dt:.1f}s]")

        results.append({
            "chi": chi_name,
            "D": D,
            "N": N,
            "lambda": lam,
            "ev_plus": float(ev_plus),
            "ev_minus": float(ev_minus),
            "gap_galerkin": float(gap_galerkin),
            "arch_plus": float(arch_plus),
            "arch_minus": float(arch_minus),
            "arch_diff": float(arch_diff),
            "S_exact": float(S_exact),
            "S_with_arch": float(S_with_arch),
            "walltime_s": float(dt),
        })

    return results


def main():
    print("=" * 75)
    print("chi_21 N-KONVERGENZ-ANALYSE (Session 7 Schritt a)")
    print("=" * 75)

    # Empirische Ground Truth
    empirical = load_empirical_gaps()

    all_results = []
    for name, D in TARGETS:
        print(f"\n{'#' * 60}\n# {name} (D={D})\n{'#' * 60}")
        if 20000 in empirical.get(name, {}):
            g_emp, N_used = empirical[name][20000]
            print(f"[emp] gap_emp (N={N_used}) = {g_emp:+.5f}")
        res = analyze_convergence(name, D, lam=20000, N_values=(200, 400, 600))
        all_results.extend(res)

    # Zusammenfassungs-Tabelle
    print("\n" + "=" * 75)
    print("KONVERGENZ-TABELLE")
    print("=" * 75)
    print(f"{'chi':>8} {'N':>4} {'gap_gal':>10} {'S_exact':>12} {'arch_diff':>12} {'S+arch':>12}")
    print("-" * 75)
    for r in all_results:
        print(f"{r['chi']:>8} {r['N']:>4} {r['gap_galerkin']:+10.5f} "
              f"{r['S_exact']:+12.5e} {r['arch_diff']:+12.5e} {r['S_with_arch']:+12.5e}")

    # Empirical reference
    print("\nReferenz (gap_emp):")
    for name, D in TARGETS:
        if 20000 in empirical.get(name, {}):
            g, N = empirical[name][20000]
            print(f"  {name}: gap_emp = {g:+.5f}  (N_used = {N})")

    RES.mkdir(exist_ok=True)
    out_json = RES / "CHI21_N_CONVERGENCE.json"
    with open(out_json, "w") as f:
        json.dump({
            "parameters": {"lambda": 20000},
            "results": all_results,
            "empirical_reference": {
                name: (empirical[name][20000] if 20000 in empirical.get(name, {}) else None)
                for name, _ in TARGETS
            },
        }, f, indent=2, default=str)
    print(f"\n[out] {out_json}")


if __name__ == "__main__":
    main()
