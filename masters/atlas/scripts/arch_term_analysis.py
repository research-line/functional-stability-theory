#!/usr/bin/env python3
# coding: utf-8
"""
Session 7 Woche 1 Option 1: Arch-Term-Integration.

Erweitert die Ground-State-Differenz-Analyse um den expliziten archimedischen
Beitrag. Die Weil-Kern-Matrix zerlegt sich in
    W = W_arch + W_prime
wobei W_arch die digamma-Diagonale ist und W_prime die Prim-Reflexions-Integrale.

Aus dem Eigenvektor-Ansatz folgt:
    ev_pm = <phi_pm, W phi_pm> = arch_pm + prime_pm
    gap_galerkin := ev_minus - ev_plus = (arch_minus - arch_plus) + (prime_minus - prime_plus)
                 =: arch_diff + prime_diff

wobei prime_diff = -2 * sum_{p,m} chi(p)^m * log p / p^(m/2) * Delta_chi(m log p) (S_exact, siehe GSD).

Ziel:
  1. Verify: S_exact_theory + arch_diff_theory ~ gap_galerkin (Konsistenz-Check).
  2. Test: [arch_diff + S_exact] als Predictor vs. gap_emp.
  3. Vergleich sign_ok, R^2 ohne/mit arch-Term.

Erwartung: arch-Term bringt predictions naeher an gap_galerkin (per definition),
und damit evtl. 9/10 sign_ok (wenn gap_galerkin selbst das Vorzeichen trifft).
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
    make_chi_D, CHARS, load_empirical_gaps,
)
from ground_state_difference_analysis import (
    autocorr, autocorr_at, compute_gap_exact,
)


# --------------------------------------------------------------------------
# Arch-Term als separate Quadratische Form
# --------------------------------------------------------------------------
def build_W_arch_only(sector, N, L, kappa=0):
    """Nur die digamma-Diagonale (arch-Teil) von build_W.

    Identische Konvention wie in all10_high_N_server.py / build_W:
      W[n,n] = digamma((2*kappa+1)/4 + i * n * pi / (2L)).real - log(pi)
    """
    bs = 0 if sector == 'cos' else 1
    idx = np.arange(N) + bs
    tau = np.pi * idx / L
    shift = (2 * kappa + 1) / 4.0
    W = np.zeros((N, N))
    diag = np.array([digamma(shift + 1j * t / 2).real - np.log(np.pi) for t in tau])
    # Normierung wie in build_W:
    if sector == 'cos':
        norm = np.where(idx > 0, 1.0 / np.sqrt(L), 1.0 / np.sqrt(2 * L))
    else:
        norm = np.full(N, 1.0 / np.sqrt(L))
    # In build_W wird die Normierung auf die Off-Diagonal via NM = outer(norm,norm) angewandt,
    # aber die digamma-Diagonale NICHT normiert (siehe np.fill_diagonal(W, diag) VOR NM).
    # Also: diagonal bleibt diag, dazu kein NM.
    np.fill_diagonal(W, diag)
    return W


def arch_contrib(coefs, sector, L, kappa=0):
    """Arch-Term <phi, W_arch phi> fuer gegebenen Grundzustand-Eigenvektor."""
    N = len(coefs)
    W_arch = build_W_arch_only(sector, N, L, kappa)
    return float(coefs @ W_arch @ coefs)


def prime_contrib_from_W(coefs, sector, N, L, primes, chi_vals, kappa=0):
    """Prim-Kern-Beitrag durch Differenz (full W) - (arch W) in der Matrix."""
    W_full = build_W(sector, N, L, primes, chi_vals, kappa)
    W_arch = build_W_arch_only(sector, N, L, kappa)
    W_prime = W_full - W_arch
    return float(coefs @ W_prime @ coefs), float(coefs @ W_full @ coefs)


# --------------------------------------------------------------------------
# Hauptroutine
# --------------------------------------------------------------------------
def analyze_with_arch(lam_target=20000, N_Galerkin=200, verbose=True):
    empirical = load_empirical_gaps()
    results = []

    L = math.log(lam_target)
    print(f"\n[analyze] lam = {lam_target}, L = {L:.4f}, N = {N_Galerkin}")

    for name, D in CHARS:
        if lam_target not in empirical.get(name, {}):
            continue
        gap_emp, N_emp = empirical[name][lam_target]
        chi_fn = make_chi_D(D)

        t0 = time.time()
        # Grundzustaende (Full W, inkl. arch + prime)
        ev_plus, coefs_plus = compute_ground_state(chi_fn, lam_target, N=N_Galerkin, sector='cos')
        ev_minus, coefs_minus = compute_ground_state(chi_fn, lam_target, N=N_Galerkin, sector='sin')

        # Arch-Beitrag (Quadratik-Form der Diagonalmatrix)
        arch_plus = arch_contrib(coefs_plus, 'cos', L)
        arch_minus = arch_contrib(coefs_minus, 'sin', L)
        arch_diff = arch_minus - arch_plus

        # Prim-Kontribution via (W_full - W_arch) — Konsistenz-Check
        primes = [p for p in sympy.primerange(2, int(lam_target) + 1) if chi_fn(p) != 0]
        cv = [chi_fn(p) for p in primes]
        prime_plus, full_plus = prime_contrib_from_W(coefs_plus, 'cos', N_Galerkin, L, primes, cv)
        prime_minus, full_minus = prime_contrib_from_W(coefs_minus, 'sin', N_Galerkin, L, primes, cv)
        prime_diff_matrix = prime_minus - prime_plus
        # full_plus sollte = ev_plus sein (Konsistenz)
        assert abs(full_plus - ev_plus) < 1e-6, f"{name}: full_plus={full_plus} vs ev_plus={ev_plus}"
        assert abs(full_minus - ev_minus) < 1e-6

        gap_galerkin = ev_minus - ev_plus  # exakter Galerkin-Gap

        # S_exact via Autokorrelation (aus GSD)
        xs_p, phi_p = phi_from_coefs(coefs_plus, L, sector='cos', n_grid=1600)
        xs_m, phi_m = phi_from_coefs(coefs_minus, L, sector='sin', n_grid=1600)
        S_exact = compute_gap_exact(chi_fn, lam_target, phi_p, xs_p, phi_m, xs_m, coeff=-2.0)

        # Kombinierter Predictor
        S_with_arch = S_exact + arch_diff
        S_with_arch_via_matrix = prime_diff_matrix + arch_diff  # = gap_galerkin (trivial)

        dt = time.time() - t0
        sign_ok_s = np.sign(S_exact) == np.sign(gap_emp)
        sign_ok_sa = np.sign(S_with_arch) == np.sign(gap_emp)
        sign_ok_gal = np.sign(gap_galerkin) == np.sign(gap_emp)
        print(f"[{name:>7}] gap_emp={gap_emp:+8.5f}, gap_gal={gap_galerkin:+.5f}, "
              f"S={S_exact:+.4e}, arch_diff={arch_diff:+.4e}, "
              f"S+arch={S_with_arch:+.4e}, "
              f"prime_diff_mat={prime_diff_matrix:+.4e} "
              f"[{'OK' if sign_ok_sa else 'FAIL'}] [{dt:.1f}s]")

        results.append({
            "chi": name,
            "D": D,
            "lambda": lam_target,
            "gap_emp": float(gap_emp),
            "N_emp": int(N_emp),
            "N_Galerkin": int(N_Galerkin),
            "gap_galerkin": float(gap_galerkin),
            "ev_plus": float(ev_plus),
            "ev_minus": float(ev_minus),
            "arch_plus": float(arch_plus),
            "arch_minus": float(arch_minus),
            "arch_diff": float(arch_diff),
            "prime_plus_matrix": float(prime_plus),
            "prime_minus_matrix": float(prime_minus),
            "prime_diff_matrix": float(prime_diff_matrix),
            "S_exact": float(S_exact),
            "S_with_arch": float(S_with_arch),
            "consistency_prime_gap": abs(prime_diff_matrix - S_exact),
            "sign_ok_S": bool(sign_ok_s),
            "sign_ok_S_arch": bool(sign_ok_sa),
            "sign_ok_galerkin": bool(sign_ok_gal),
        })

    return results


def correlation_stats(values, refs):
    x = np.array(values)
    y = np.array(refs)
    if np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return float("nan"), float("nan"), 0
    r = float(np.corrcoef(x, y)[0, 1])
    sign_ok = int(np.sum(np.sign(x) == np.sign(y)))
    return r, r * r, sign_ok


def main():
    print("=" * 75)
    print("ARCH-TERM INTEGRATION ANALYSIS (Session 7 Woche 1 Option 1)")
    print("=" * 75)

    LAM = 20000
    N_G = 200
    results = analyze_with_arch(lam_target=LAM, N_Galerkin=N_G, verbose=True)

    emp = [r["gap_emp"] for r in results]
    print("\n" + "=" * 75)
    print("KORRELATIONS-ANALYSE")
    print("=" * 75)

    for field in ["S_exact", "arch_diff", "S_with_arch", "gap_galerkin",
                  "prime_diff_matrix"]:
        vals = [r[field] for r in results]
        r, r2, sok = correlation_stats(vals, emp)
        print(f"  {field:>22}:  R = {r:+.4f},  R^2 = {r2:.4f},  sign_ok = {sok}/{len(emp)}")

    # Konsistenz-Check: prime_diff_matrix sollte = S_exact sein
    print("\n" + "=" * 75)
    print("KONSISTENZ: prime_diff_matrix - S_exact (sollte ~ 0)")
    print("=" * 75)
    for r in results:
        print(f"  {r['chi']:>8}: prime_diff_mat={r['prime_diff_matrix']:+.5e}, "
              f"S_exact={r['S_exact']:+.5e}, diff={r['consistency_prime_gap']:.3e}")

    # JSON + MD
    RES.mkdir(exist_ok=True)
    with open(RES / "ARCH_TERM_ANALYSIS.json", "w") as f:
        json.dump({
            "parameters": {"lambda": LAM, "N_Galerkin": N_G},
            "results": results,
        }, f, indent=2)
    print(f"\n[out] {RES / 'ARCH_TERM_ANALYSIS.json'}")

    # Markdown Summary
    s_arch_r, s_arch_r2, s_arch_sok = correlation_stats(
        [r["S_with_arch"] for r in results], emp)
    s_r, s_r2, s_sok = correlation_stats([r["S_exact"] for r in results], emp)
    gal_r, gal_r2, gal_sok = correlation_stats([r["gap_galerkin"] for r in results], emp)

    md = []
    md.append("# Arch-Term-Integration — Session 7 Woche 1 Option 1")
    md.append("")
    md.append(f"**Datum:** {time.strftime('%Y-%m-%d')}")
    md.append(f"**Parameter:** lambda = {LAM}, N_Galerkin = {N_G}")
    md.append("")
    md.append("## Zerlegung")
    md.append("")
    md.append("gap_galerkin := ev_minus - ev_plus = arch_diff + prime_diff")
    md.append("- arch_diff = <phi_minus, W_arch_minus phi_minus> - <phi_plus, W_arch_plus phi_plus>")
    md.append("- prime_diff = S_exact (aus Autokorrelations-Formel, coeff=-2)")
    md.append("")
    md.append("## Predictor-Vergleich")
    md.append("")
    md.append("| Predictor | sign_ok | R | R^2 |")
    md.append("|-----------|---------|---|-----|")
    md.append(f"| S_exact (Prim nur, aus GSD) | {s_sok}/10 | {s_r:+.4f} | {s_r2:.4f} |")
    md.append(f"| S_with_arch (Prim + Arch) | {s_arch_sok}/10 | {s_arch_r:+.4f} | {s_arch_r2:.4f} |")
    md.append(f"| gap_galerkin (voll, exakt bei N={N_G}) | {gal_sok}/10 | {gal_r:+.4f} | {gal_r2:.4f} |")
    md.append("")
    md.append("## Tabelle")
    md.append("")
    md.append("| chi | D | gap_emp | gap_gal | S_exact | arch_diff | S+arch | sign+arch |")
    md.append("|-----|---|---------|---------|---------|-----------|--------|-----------|")
    for r in results:
        md.append(f"| {r['chi']} | {r['D']} | {r['gap_emp']:+.5f} | {r['gap_galerkin']:+.5f} | "
                  f"{r['S_exact']:+.4e} | {r['arch_diff']:+.4e} | {r['S_with_arch']:+.4e} | "
                  f"{'OK' if r['sign_ok_S_arch'] else 'FAIL'} |")
    md.append("")
    md.append("## Interpretation")
    md.append("")
    if s_arch_sok > s_sok:
        md.append(f"Arch-Term verbessert sign_ok von {s_sok}/10 auf {s_arch_sok}/10.")
    elif s_arch_sok == s_sok:
        md.append(f"Arch-Term aendert sign_ok nicht ({s_sok}/10). Vorzeichen-Information bereits in S_exact.")
    else:
        md.append(f"Arch-Term verschlechtert sign_ok auf {s_arch_sok}/10 (von {s_sok}/10).")
    if abs(s_arch_r2 - gal_r2) < 0.02:
        md.append("S_with_arch rekonstruiert gap_galerkin (erwartet: beide decomponieren denselben Eigenwert).")
    md.append("")

    with open(RES / "ARCH_TERM_ANALYSIS.md", "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print(f"[out] {RES / 'ARCH_TERM_ANALYSIS.md'}")


if __name__ == "__main__":
    main()
