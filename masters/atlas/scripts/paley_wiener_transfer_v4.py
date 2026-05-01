#!/usr/bin/env python3
# coding: utf-8
"""
Session 11 Paley-Wiener Transfer Validation v4.

Uses the rigorous Thm 5.1 formula from ANALYTIC_KERNEL_V2.md:
    gap_chi^{prim} = -2 sum_{p,m} chi(p)^m log(p)/p^{m/2}
                     * [(phi+ * phi+)(m log p) + (phi- * phi-)(m log p)]

With Hermite trial functions (character-independent, archimedean ansatz):
    phi+ = (pi sigma^2)^{-1/4} exp(-x^2/(2 sigma^2))
    phi- = sqrt(2/(sigma^3 sqrt(pi))) * x * exp(-x^2/(2 sigma^2))

The convolutions are computed analytically:
    (phi+ * phi+)(u) = exp(-u^2 / (4 sigma^2))
    (phi- * phi-)(u) = (u^2/(2 sigma^2) - 1) exp(-u^2 / (4 sigma^2))

Sum:
    b(u) = (phi+ * phi+)(u) + (phi- * phi-)(u) = (u^2 / (2 sigma^2)) exp(-u^2 / (4 sigma^2))

The character-signal lives in the sign-weighted prime-power sum
    S_chi = sum_{p,m} chi(p)^m (log p) / p^{m/2} * (m log p)^2 / (2 sigma^2)
            * exp(-(m log p)^2 / (4 sigma^2))

Sweep sigma; compare sign(gap^{prim}) := -2 * S_chi with empirical gap.

Critical test: chi_33 must give gap^{prim} < 0 (i.e. S_chi > 0).
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
    ("chi_5", 5), ("chi_8", 8), ("chi_12", 12), ("chi_13", 13),
    ("chi_17", 17), ("chi_21", 21), ("chi_24", 24), ("chi_29", 29),
    ("chi_33", 33), ("chi_60", 60),
]


def load_empirical_gaps():
    path = RES / "ARCH_TERM_N600_SERVER.json"
    with path.open(encoding="utf-8") as f:
        data = json.load(f)
    gaps = {}
    for entry in data["results"]:
        gaps[entry["chi"]] = entry["gap_galerkin"]
    return gaps


def collect_prime_terms(chi_fn, D, log_pm_max, max_m=8):
    terms = []
    p_max = int(math.exp(log_pm_max)) + 1
    for p in sympy.primerange(2, p_max + 1):
        if D % p == 0:
            continue
        cp = chi_fn(p)
        if cp == 0:
            continue
        log_p = math.log(p)
        for m in range(1, max_m + 1):
            log_pm = m * log_p
            if log_pm > log_pm_max:
                break
            c_pm = cp ** m
            weight = log_p / (p ** (m / 2.0))
            terms.append((c_pm, log_pm, weight))
    return terms


def bump_hermite(u, sigma):
    """(phi+ * phi+)(u) + (phi- * phi-)(u) with Hermite trial functions."""
    # = (u^2 / (2 sigma^2)) * exp(-u^2/(4 sigma^2))
    return (u * u / (2.0 * sigma * sigma)) * np.exp(-(u * u) / (4.0 * sigma * sigma))


def compute_gap_prim(terms, sigma):
    """gap_prim = -2 sum chi(p)^m (log p / p^{m/2}) * bump(m log p)."""
    c_arr = np.array([t[0] for t in terms], dtype=float)
    lp_arr = np.array([t[1] for t in terms], dtype=float)
    w_arr = np.array([t[2] for t in terms], dtype=float)
    b = bump_hermite(lp_arr, sigma)
    s = float(np.sum(c_arr * w_arr * b))
    return -2.0 * s


def pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mx = x.mean(); my = y.mean()
    sx = x.std(); sy = y.std()
    if sx == 0 or sy == 0:
        return 0.0
    return float(((x - mx) * (y - my)).mean() / (sx * sy))


def aggregate(rows):
    gaps = np.array([r["gap_emp"] for r in rows])
    vals = np.array([r["gap_prim"] for r in rows])
    sign_ok = sum(int(np.sign(r["gap_prim"]) == np.sign(r["gap_emp"]))
                  for r in rows if r["gap_emp"] != 0)
    R = pearson(vals, gaps)
    idx = [i for i, r in enumerate(rows) if r["chi"] != "chi_21"]
    gaps_n21 = gaps[idx]; vals_n21 = vals[idx]
    sign_ok_n21 = sum(int(np.sign(rows[i]["gap_prim"]) == np.sign(rows[i]["gap_emp"]))
                      for i in idx if rows[i]["gap_emp"] != 0)
    R_n21 = pearson(vals_n21, gaps_n21)
    return {
        "sign_ok_all": sign_ok, "R_all": R, "R2_all": R * R,
        "sign_ok_n21": sign_ok_n21, "R_n21": R_n21, "R2_n21": R_n21 * R_n21,
    }


def main():
    lam = 20000
    L = math.log(lam)
    print(f"[setup] lambda={lam}, L={L:.4f}")
    emp = load_empirical_gaps()

    log_pm_max = 9.0
    print(f"[setup] log_pm_max = {log_pm_max} (p^m <= {math.exp(log_pm_max):.0f})")

    t0 = time.time()
    terms_per_chi = {}
    for name, D in CHARS:
        chi = make_chi_D(D)
        terms = collect_prime_terms(chi, D, log_pm_max, max_m=8)
        terms_per_chi[name] = terms
    print(f"[prime collection] {time.time()-t0:.2f} s")

    # Sweep sigma from 0.5 to 10
    sigmas = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0]

    all_results = []
    for sigma in sigmas:
        rows = []
        for name, D in CHARS:
            terms = terms_per_chi[name]
            gap_prim = compute_gap_prim(terms, sigma)
            rows.append({
                "chi": name, "D": D,
                "gap_prim": gap_prim,
                "gap_emp": emp[name],
            })
        agg = aggregate(rows)
        print(f"[sigma={sigma:5.2f}] "
              f"all: {agg['sign_ok_all']}/10 R2={agg['R2_all']:.3f} | "
              f"n21: {agg['sign_ok_n21']}/9 R2={agg['R2_n21']:.3f}")
        all_results.append({"sigma": sigma, "rows": rows, "aggregate": agg})

    # Best by (sign_ok_all, R2_all)
    best = max(all_results, key=lambda r: (
        r["aggregate"]["sign_ok_all"], r["aggregate"]["R2_all"]))
    print(f"\n[best sigma] {best['sigma']}: "
          f"{best['aggregate']['sign_ok_all']}/10, "
          f"R2={best['aggregate']['R2_all']:.3f}")
    print(f"\n[per-character at best sigma]")
    print(f"  {'chi':<10}{'D':>5}{'gap_prim':>14}{'gap_emp':>14}{'ratio':>10}{'sign':>6}")
    for r in best["rows"]:
        ratio = r["gap_prim"] / r["gap_emp"] if r["gap_emp"] != 0 else float("inf")
        ok = "OK" if np.sign(r["gap_prim"]) == np.sign(r["gap_emp"]) else "FAIL"
        print(f"  {r['chi']:<10}{r['D']:>5}{r['gap_prim']:>14.4f}"
              f"{r['gap_emp']:>14.5f}{ratio:>10.2f}{ok:>6}")

    # Save
    out_json = RES / "PALEY_WIENER_TRANSFER_V4.json"
    with out_json.open("w", encoding="utf-8") as f:
        json.dump({
            "parameters": {"lambda": lam, "L": L, "log_pm_max": log_pm_max, "sigmas": sigmas},
            "sigma_sweep": all_results,
            "best_sigma": best["sigma"],
        }, f, indent=2)
    print(f"[write] {out_json}")

    md_path = RES / "PALEY_WIENER_TRANSFER_V4.md"
    with md_path.open("w", encoding="utf-8") as f:
        f.write("# Paley-Wiener Transfer Validation v4\n\n")
        f.write("**Date:** 2026-04-17 (Session 11)\n")
        f.write("**Script:** `_scripts/paley_wiener_transfer_v4.py`\n")
        f.write("**Approach:** Direct Thm 5.1 formula (ANALYTIC_KERNEL_V2) with Hermite trial functions.\n\n")
        f.write("## Formula\n\n")
        f.write("gap_chi^{prim} = -2 sum_{p,m} chi(p)^m (log p / p^{m/2}) * b(m log p)\n\n")
        f.write("with b(u) = (phi+ * phi+)(u) + (phi- * phi-)(u) "
                "= (u^2/(2 sigma^2)) exp(-u^2/(4 sigma^2))\n\n")
        f.write("## Sigma sweep\n\n")
        f.write("| sigma | sign(all) | R²(all) | sign(n21) | R²(n21) |\n")
        f.write("|---:|:---:|---:|:---:|---:|\n")
        for r in all_results:
            agg = r["aggregate"]
            f.write(f"| {r['sigma']:.2f} | {agg['sign_ok_all']}/10 | "
                    f"{agg['R2_all']:.3f} | {agg['sign_ok_n21']}/9 | "
                    f"{agg['R2_n21']:.3f} |\n")
        f.write(f"\n**Best sigma**: {best['sigma']}: "
                f"{best['aggregate']['sign_ok_all']}/10, "
                f"R²={best['aggregate']['R2_all']:.3f}\n\n")
        f.write(f"## Per-character detail at best sigma\n\n")
        f.write("| chi | D | gap_prim | gap_emp | ratio | sign |\n")
        f.write("|---|---|---:|---:|---:|:---:|\n")
        for r in best["rows"]:
            ratio = r["gap_prim"] / r["gap_emp"] if r["gap_emp"] != 0 else float("inf")
            ok = "OK" if np.sign(r["gap_prim"]) == np.sign(r["gap_emp"]) else "**FAIL**"
            f.write(f"| {r['chi']} | {r['D']} | {r['gap_prim']:+.4f} | "
                    f"{r['gap_emp']:+.5f} | {ratio:+.2f} | {ok} |\n")
    print(f"[write] {md_path}")


if __name__ == "__main__":
    main()
