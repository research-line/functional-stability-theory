#!/usr/bin/env python3
# coding: utf-8
"""
Session 8 - Paley-Wiener Transfer Validation Script.

Tests the perturbative sector-difference formula from
PALEY_WIENER_DIRICHLET_TRANSFER.md §5 against the empirical N=600 gaps.

Core formula (first-order in eps_L = L^{-1/2}):

    gap_chi ~ eps_L * A_chi(L),
    A_chi(L) = sum_{(p,m),(p',m')} chi(p)^m chi(p')^{m'}
               * log(p)*log(p') / sqrt(p^m * (p')^{m'})
               * Phi(m*log p - m'*log p'),

where Phi(u) = (2*pi*sigma_Phi^2)^{-1/2} * exp(-u^2 / (4*sigma_Phi^2)) is the
archimedean-ground-state autocorrelation Gaussian of width sigma_Phi,
explicitly sigma_Phi^2 = sqrt(L/4.21) (from -psi''(1/4)/2 ~ 4.21).

The diagonal (p,m) == (p',m') part gives a character-UNIVERSAL term
(since |chi(p)^m|^2 = 1 for p nmid q). The character-specific signature
lives in the OFF-DIAGONAL contribution

    A_chi^off(L) = A_chi(L) - A_chi^diag(L).

Prediction (Conjecture C2): sign(A_chi^off) agrees with sign(gap_chi)
for 9 out of 10 test characters, in particular sign(A_chi33^off) = -1.

Inputs
------
  _results/ARCH_TERM_N600_SERVER.json     empirical gaps at lambda=20000, N=600

Outputs
-------
  _results/PALEY_WIENER_TRANSFER_VALIDATION.json   raw data
  _results/PALEY_WIENER_TRANSFER_VALIDATION.md     human report

The script is self-contained: it does not import Galerkin machinery. Only
the empirical gap values (N=600) are read from disk.
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


# --------------------------------------------------------------------------
# Kronecker-Symbol (Standard, identisch zu anderen DIRICHLET_L-Scripts)
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
# Archimedean autocorrelation kernel Phi(u)
# --------------------------------------------------------------------------
def sigma_Phi_squared(L):
    """sigma_Phi^2 = sqrt(L / K_arch''(0)/2), with -psi''(1/4)/2 ~= 4.21.

    Reference: ANALYTIC_GROUNDSTATE.md §2.3 and §4.1 (Gaussian width of
    archimedean ground state at Fourier minimum xi=0).
    """
    Karch_pp0 = 4.21  # = -psi''(1/4)/2, numerical value
    return math.sqrt(L / Karch_pp0)


def Phi_autocorr(u, sigma2):
    """Archimedean ground-state autocorrelation (Gaussian).

    Phi(u) = (2*pi*sigma2)^{-1/2} * exp(-u^2 / (4*sigma2))
    This is the self-convolution of a Gaussian of width sqrt(sigma2).
    Width of the autocorrelation: sqrt(2*sigma2).
    """
    return math.exp(-(u * u) / (4.0 * sigma2)) / math.sqrt(2.0 * math.pi * sigma2)


# --------------------------------------------------------------------------
# Prime-pair enumeration
# --------------------------------------------------------------------------
def collect_prime_terms(chi_fn, D, L, max_m=8, phi_width_cutoff=6.0):
    """Collect all (p, m) such that p nmid D, m*log(p) is within the
    effective support of Phi (|u| <= phi_width_cutoff * sqrt(sigma_Phi^2)).

    Returns list of tuples (c_pm, log_p_m, weight), where
        c_pm    = chi(p)^m                      (complex, but real here)
        log_p_m = m * log(p)
        weight  = log(p) / p^(m/2)
    """
    sigma2 = sigma_Phi_squared(L)
    u_max = phi_width_cutoff * math.sqrt(sigma2)
    # |log(p^m)| <= u_max means p^m <= exp(u_max)
    pm_max = math.exp(u_max)
    terms = []
    p_max_m1 = int(pm_max) + 1
    for p in sympy.primerange(2, p_max_m1 + 1):
        if D % p == 0:
            continue
        cp = chi_fn(p)
        if cp == 0:
            continue
        log_p = math.log(p)
        for m in range(1, max_m + 1):
            log_pm = m * log_p
            if log_pm > u_max:
                break
            c_pm = cp ** m  # character value on prime power
            p_half_m = p ** (m / 2.0)
            weight = log_p / p_half_m
            terms.append((c_pm, log_pm, weight))
    return terms


# --------------------------------------------------------------------------
# A_chi computation: diagonal and off-diagonal
# --------------------------------------------------------------------------
def compute_A_chi(chi_fn, D, L, max_m=8, phi_width_cutoff=6.0):
    """Compute the perturbation functional A_chi(L).

    Returns dict with diagonal, off-diagonal, and total components.
    """
    terms = collect_prime_terms(chi_fn, D, L, max_m=max_m,
                                phi_width_cutoff=phi_width_cutoff)
    sigma2 = sigma_Phi_squared(L)

    # Vectorized computation
    n = len(terms)
    if n == 0:
        return {"n_terms": 0, "A_diag": 0.0, "A_off": 0.0, "A_total": 0.0}

    c_arr = np.array([t[0] for t in terms], dtype=float)
    log_arr = np.array([t[1] for t in terms], dtype=float)
    w_arr = np.array([t[2] for t in terms], dtype=float)

    # Pairwise Phi(u_i - u_j): full n x n matrix
    diff = log_arr[:, None] - log_arr[None, :]
    K = np.exp(-(diff * diff) / (4.0 * sigma2)) / math.sqrt(2.0 * math.pi * sigma2)
    C = np.outer(c_arr, c_arr)
    W = np.outer(w_arr, w_arr)

    M = C * W * K  # pairwise contribution matrix

    A_total = M.sum()
    A_diag = np.diag(M).sum()
    A_off = A_total - A_diag

    return {
        "n_terms": n,
        "A_diag": float(A_diag),
        "A_off": float(A_off),
        "A_total": float(A_total),
        "ratio_off_over_total": float(A_off / A_total) if A_total != 0 else 0.0,
    }


# --------------------------------------------------------------------------
# Load empirical N=600 gaps
# --------------------------------------------------------------------------
def load_empirical_gaps():
    path = RES / "ARCH_TERM_N600_SERVER.json"
    with path.open(encoding="utf-8") as f:
        data = json.load(f)
    gaps = {}
    for entry in data["results"]:
        gaps[entry["chi"]] = {
            "D": entry["D"],
            "gap_galerkin": entry["gap_galerkin"],
            "arch_diff": entry["arch_diff"],
            "S_exact": entry["S_exact"],
            "S_with_arch": entry["S_with_arch"],
        }
    return gaps


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    lam = 20000
    L = math.log(lam)
    sigma2 = sigma_Phi_squared(L)
    eps_L = L ** (-0.5)
    print(f"[setup] lambda = {lam}, L = {L:.4f}")
    print(f"[setup] sigma_Phi^2 = {sigma2:.4f} (width ~= {math.sqrt(sigma2):.3f})")
    print(f"[setup] eps_L = L^{{-1/2}} = {eps_L:.4f}")

    emp = load_empirical_gaps()

    rows = []
    t0 = time.time()
    for name, D in CHARS:
        chi = make_chi_D(D)
        t_start = time.time()
        res = compute_A_chi(chi, D, L, max_m=8, phi_width_cutoff=6.0)
        walltime = time.time() - t_start
        gap = emp[name]["gap_galerkin"]
        sign_emp = 1 if gap > 0 else (-1 if gap < 0 else 0)

        # Predictions from the three components
        A_total = res["A_total"]
        A_off = res["A_off"]
        A_diag = res["A_diag"]

        sign_total_ok = int(np.sign(A_total) == sign_emp)
        sign_off_ok = int(np.sign(A_off) == sign_emp)

        # Predicted gap value (off-diagonal only, times eps_L)
        pred_off = eps_L * A_off
        pred_total = eps_L * A_total

        row = {
            "chi": name,
            "D": D,
            "n_terms": res["n_terms"],
            "A_diag": A_diag,
            "A_off": A_off,
            "A_total": A_total,
            "ratio_off_over_total": res["ratio_off_over_total"],
            "pred_off_scaled": pred_off,
            "pred_total_scaled": pred_total,
            "gap_empirical": gap,
            "sign_off_ok": sign_off_ok,
            "sign_total_ok": sign_total_ok,
            "walltime_s": walltime,
        }
        rows.append(row)
        print(f"  [{name}] D={D}, n_terms={res['n_terms']}, "
              f"A_diag={A_diag:+.5f}, A_off={A_off:+.5f}, A_total={A_total:+.5f}, "
              f"gap={gap:+.5f}, sign_off={sign_off_ok}, sign_total={sign_total_ok}, "
              f"{walltime:.2f}s")

    total_walltime = time.time() - t0

    # Aggregate
    n = len(rows)
    sign_off_sum = sum(r["sign_off_ok"] for r in rows)
    sign_total_sum = sum(r["sign_total_ok"] for r in rows)

    # Correlation (Pearson) of A_off vs gap
    gaps = np.array([r["gap_empirical"] for r in rows])
    A_offs = np.array([r["A_off"] for r in rows])
    A_totals = np.array([r["A_total"] for r in rows])

    def pearson(x, y):
        mx = x.mean(); my = y.mean()
        sx = x.std(); sy = y.std()
        if sx == 0 or sy == 0:
            return 0.0
        return float(((x - mx) * (y - my)).mean() / (sx * sy))

    R_off = pearson(A_offs, gaps)
    R_total = pearson(A_totals, gaps)

    # Without chi_21 (known N-oscillator)
    idx_no21 = [i for i, r in enumerate(rows) if r["chi"] != "chi_21"]
    gaps_n21 = gaps[idx_no21]
    A_offs_n21 = A_offs[idx_no21]
    A_totals_n21 = A_totals[idx_no21]
    R_off_n21 = pearson(A_offs_n21, gaps_n21)
    R_total_n21 = pearson(A_totals_n21, gaps_n21)
    sign_off_n21 = sum(rows[i]["sign_off_ok"] for i in idx_no21)
    sign_total_n21 = sum(rows[i]["sign_total_ok"] for i in idx_no21)

    summary = {
        "parameters": {
            "lambda": lam,
            "L": L,
            "sigma_Phi_squared": sigma2,
            "eps_L": eps_L,
            "max_m": 8,
            "phi_width_cutoff": 6.0,
            "total_walltime_s": total_walltime,
        },
        "rows": rows,
        "aggregate": {
            "all_10": {
                "sign_off_ok": sign_off_sum,
                "sign_total_ok": sign_total_sum,
                "R_off": R_off,
                "R_off_squared": R_off ** 2,
                "R_total": R_total,
                "R_total_squared": R_total ** 2,
            },
            "without_chi_21": {
                "sign_off_ok": sign_off_n21,
                "sign_total_ok": sign_total_n21,
                "R_off": R_off_n21,
                "R_off_squared": R_off_n21 ** 2,
                "R_total": R_total_n21,
                "R_total_squared": R_total_n21 ** 2,
            },
        },
    }

    # Write JSON
    out_json = RES / "PALEY_WIENER_TRANSFER_VALIDATION.json"
    with out_json.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"[write] {out_json}")

    # Write MD report
    md_path = RES / "PALEY_WIENER_TRANSFER_VALIDATION.md"
    with md_path.open("w", encoding="utf-8") as f:
        f.write("# Paley-Wiener Transfer Validation (Session 8)\n\n")
        f.write(f"**Date:** 2026-04-16\n")
        f.write(f"**Script:** `_scripts/paley_wiener_transfer_validation.py`\n")
        f.write(f"**Source:** `PALEY_WIENER_DIRICHLET_TRANSFER.md` §6.3\n\n")
        f.write(f"## Parameters\n\n")
        f.write(f"- lambda = {lam}, L = {L:.4f}\n")
        f.write(f"- sigma_Phi^2 = {sigma2:.4f}\n")
        f.write(f"- eps_L = L^(-1/2) = {eps_L:.4f}\n")
        f.write(f"- max_m = 8, phi_width_cutoff = 6.0\n")
        f.write(f"- total walltime: {total_walltime:.2f} s\n\n")
        f.write("## Results per character\n\n")
        f.write("| chi | D | n_terms | A_diag | A_off | A_total | gap_emp | sign(A_off) | sign(A_total) |\n")
        f.write("|---|---|---|---:|---:|---:|---:|:---:|:---:|\n")
        for r in rows:
            ok_off = "OK" if r["sign_off_ok"] else "**FAIL**"
            ok_tot = "OK" if r["sign_total_ok"] else "**FAIL**"
            f.write(f"| {r['chi']} | {r['D']} | {r['n_terms']} | "
                    f"{r['A_diag']:+.4f} | {r['A_off']:+.4f} | {r['A_total']:+.4f} | "
                    f"{r['gap_empirical']:+.5f} | {ok_off} | {ok_tot} |\n")
        f.write("\n## Aggregate\n\n")
        f.write("### All 10 characters\n")
        f.write(f"- sign(A_off) ok: **{sign_off_sum}/10**\n")
        f.write(f"- sign(A_total) ok: **{sign_total_sum}/10**\n")
        f.write(f"- Pearson R(A_off, gap): **{R_off:+.4f}**, R^2 = {R_off**2:.4f}\n")
        f.write(f"- Pearson R(A_total, gap): **{R_total:+.4f}**, R^2 = {R_total**2:.4f}\n\n")
        f.write("### Without chi_21 (known N-oscillator)\n")
        f.write(f"- sign(A_off) ok: **{sign_off_n21}/9**\n")
        f.write(f"- sign(A_total) ok: **{sign_total_n21}/9**\n")
        f.write(f"- Pearson R(A_off, gap): **{R_off_n21:+.4f}**, R^2 = {R_off_n21**2:.4f}\n")
        f.write(f"- Pearson R(A_total, gap): **{R_total_n21:+.4f}**, R^2 = {R_total_n21**2:.4f}\n\n")
        f.write("## Critical test\n\n")
        chi33 = next(r for r in rows if r["chi"] == "chi_33")
        crit = "PASSED" if chi33["sign_off_ok"] else "FAILED"
        f.write(f"**chi_33 prediction:** A_off = {chi33['A_off']:+.5f}, "
                f"gap_emp = {chi33['gap_empirical']:+.5f} -> **{crit}**\n\n")
        f.write("## Success criterion\n\n")
        f.write(f"Conjecture C2 (PALEY_WIENER_DIRICHLET_TRANSFER §6.2):\n")
        f.write(f"  sign(A_off) agreement >= 9/10 AND chi_33 negative.\n\n")
        success = sign_off_sum >= 9 and chi33["sign_off_ok"]
        f.write(f"**Result: {'CONFIRMED' if success else 'NOT CONFIRMED'}** "
                f"(sign_off = {sign_off_sum}/10, chi_33 sign = {'OK' if chi33['sign_off_ok'] else 'FAIL'})\n")
    print(f"[write] {md_path}")
    print(f"[done] total walltime {total_walltime:.2f} s")
    print(f"[summary] sign_off = {sign_off_sum}/10, chi_33 sign ok = {chi33['sign_off_ok']}")


if __name__ == "__main__":
    main()
