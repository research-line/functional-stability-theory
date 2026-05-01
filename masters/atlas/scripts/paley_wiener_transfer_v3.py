#!/usr/bin/env python3
# coding: utf-8
"""
Session 11 Paley-Wiener Transfer Validation v3 with FULL character-specific
Weil-multiplier resolvent kernel.

Motivation:
v1 (Gaussian kernel):     4/10 signs, R²=0.03 — character-invariant kernel
v2 (pole-cos heuristic):  6/10 signs, R²=0.31 (best) — ad-hoc weights
v3 (this):                Full multiplier resolvent 1/(K_chi_hat(xi) - mu_chi).

The v3 kernel is built from the exact character-specific Weil multiplier
    K_chi_hat(xi) = h_arch(xi) + 2 * sum_{p, m} w_{p,m} cos(m*log(p) * xi)
with w_{p,m} = chi(p)^m * log(p) / p^{m/2} (real characters).

Ground-state energy mu_chi is estimated as min_xi K_chi_hat(xi).
Resolvent R_chi(xi) = 1/(K_chi_hat(xi) - mu_chi + i*epsilon) is Fourier-inverted
to a position-space kernel K_L^chi(u) that includes character-specific
pole structure AT ZEROS OF THE MULTIPLIER (which, via Weil formula,
track the character L-zeros).

Then:
    A_chi^{off} = sum_{p != p', m, m'} chi(p)^m chi(p')^{m'}
                  * log(p) log(p') / p^{m/2}(p')^{m'/2}
                  * Re K_L^chi(m*log p - m'*log p')

Critical test: chi_33 must come out negative.
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import sympy
from scipy.special import digamma as scipy_digamma

sys.stdout.reconfigure(line_buffering=True)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "_results"


# --------------------------------------------------------------------------
# Kronecker symbol
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
# Data loaders
# --------------------------------------------------------------------------
def load_empirical_gaps():
    path = RES / "ARCH_TERM_N600_SERVER.json"
    with path.open(encoding="utf-8") as f:
        data = json.load(f)
    gaps = {}
    for entry in data["results"]:
        gaps[entry["chi"]] = entry["gap_galerkin"]
    return gaps


# --------------------------------------------------------------------------
# Archimedean multiplier h_arch(xi) = log(q/pi) + Re psi(1/4 + i xi/2)
# --------------------------------------------------------------------------
def h_arch(xi_grid, q):
    # scipy digamma is complex; take real part
    z = 0.25 + 0.5j * xi_grid
    return np.log(q / math.pi) + np.real(scipy_digamma(z))


# --------------------------------------------------------------------------
# Prime-collection
# --------------------------------------------------------------------------
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


# --------------------------------------------------------------------------
# v3 kernel: multiplier resolvent
# --------------------------------------------------------------------------
def build_multiplier(xi_grid, terms, q):
    """K_chi_hat(xi) = h_arch(xi) + 2 sum w_{p,m} cos(m log p * xi)."""
    K_hat = h_arch(xi_grid, q)
    for c_pm, log_pm, weight in terms:
        K_hat = K_hat + 2.0 * c_pm * weight * np.cos(log_pm * xi_grid)
    return K_hat


def build_kernel_v3(xi_grid, K_hat, epsilon, u_eval):
    """K_L^chi(u) = Fourier inverse of 1/(K_hat - mu + i*epsilon) at u in u_eval.

    Returns a real-valued array K_L[u] of same length as u_eval.
    """
    mu = float(K_hat.min())
    # Regularized resolvent (kept complex; we take Re of the integral)
    denom = (K_hat - mu) + 1j * epsilon
    R = 1.0 / denom
    # Inverse Fourier transform via direct summation
    # K_L(u) = int R(xi) exp(i xi u) dxi / (2 pi)
    # Discretize: sum_k R_k exp(i xi_k u) * dxi / (2 pi)
    dxi = xi_grid[1] - xi_grid[0]
    # For each u_eval, K_L(u) = dxi / (2 pi) * sum_k R_k exp(i xi_k u)
    # Use matmul
    # shape: u_eval (nU,) , xi_grid (nXi,)
    phase = np.exp(1j * np.outer(u_eval, xi_grid))  # (nU, nXi)
    K_u = phase @ R  # (nU,)
    K_u = K_u * (dxi / (2.0 * math.pi))
    return np.real(K_u), mu


# --------------------------------------------------------------------------
# Compute A_chi from a kernel function K_L(u)
# --------------------------------------------------------------------------
def compute_A_from_kernel(terms, K_L_values, u_eval):
    """A = sum_{i,j} chi_i chi_j w_i w_j K_L(log_pm_i - log_pm_j).

    We interpolate K_L onto the pairwise differences.
    """
    c_arr = np.array([t[0] for t in terms], dtype=float)
    lp_arr = np.array([t[1] for t in terms], dtype=float)
    w_arr = np.array([t[2] for t in terms], dtype=float)
    # Pairwise difference u_ij = lp_i - lp_j
    u_matrix = lp_arr[:, None] - lp_arr[None, :]
    # Interpolate K_L onto u_matrix
    K_flat = np.interp(u_matrix.ravel(), u_eval, K_L_values)
    K_matrix = K_flat.reshape(u_matrix.shape)
    C = np.outer(c_arr, c_arr)
    W = np.outer(w_arr, w_arr)
    M = C * W * K_matrix
    A_total = float(M.sum())
    A_diag = float(np.diag(M).sum())
    A_off = A_total - A_diag
    return {"n_terms": len(terms), "A_diag": A_diag, "A_off": A_off,
            "A_total": A_total}


def pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mx = x.mean(); my = y.mean()
    sx = x.std(); sy = y.std()
    if sx == 0 or sy == 0:
        return 0.0
    return float(((x - mx) * (y - my)).mean() / (sx * sy))


def aggregate(rows, field):
    gaps = np.array([r["gap_emp"] for r in rows])
    vals = np.array([r[field] for r in rows])
    sign_ok = sum(int(np.sign(r[field]) == np.sign(r["gap_emp"]))
                  for r in rows if r["gap_emp"] != 0)
    R = pearson(vals, gaps)
    idx = [i for i, r in enumerate(rows) if r["chi"] != "chi_21"]
    gaps_n21 = gaps[idx]; vals_n21 = vals[idx]
    sign_ok_n21 = sum(int(np.sign(rows[i][field]) == np.sign(rows[i]["gap_emp"]))
                      for i in idx if rows[i]["gap_emp"] != 0)
    R_n21 = pearson(vals_n21, gaps_n21)
    return {
        "sign_ok_all": sign_ok, "R_all": R, "R2_all": R * R,
        "sign_ok_n21": sign_ok_n21, "R_n21": R_n21, "R2_n21": R_n21 * R_n21,
    }


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    lam = 20000
    L = math.log(lam)
    eps_L = 1.0 / math.sqrt(L)
    print(f"[setup] lambda={lam}, L={L:.4f}, eps_L={eps_L:.4f}")

    emp = load_empirical_gaps()

    log_pm_max = 9.0
    print(f"[setup] log_pm_max = {log_pm_max} (p^m <= {math.exp(log_pm_max):.0f})")

    # xi-grid for multiplier evaluation.
    # Range: up to first few zeros (~20 is enough for chi with gamma_1 ~ 2-6).
    xi_max = 30.0
    n_xi = 8192  # sufficient for FFT-like accuracy
    xi_grid = np.linspace(-xi_max, xi_max, n_xi)
    print(f"[setup] xi-grid: [-{xi_max}, {xi_max}], n_xi={n_xi}, dxi={xi_grid[1]-xi_grid[0]:.5f}")

    # u-grid for position kernel (covers all possible log-p-differences)
    u_max = log_pm_max * 2.0  # ~18
    n_u = 2048
    u_eval = np.linspace(-u_max, u_max, n_u)
    print(f"[setup] u-grid: [-{u_max}, {u_max}], n_u={n_u}")

    # Regularization
    # Range swept to find best value
    epsilons = [0.1, 0.3, 0.5, 1.0, 2.0]

    # Precompute prime terms per character
    t0 = time.time()
    terms_per_chi = {}
    for name, D in CHARS:
        chi = make_chi_D(D)
        terms = collect_prime_terms(chi, D, log_pm_max, max_m=8)
        terms_per_chi[name] = terms
        print(f"  [{name}] D={D}, n_terms={len(terms)}")
    print(f"[prime collection] {time.time()-t0:.2f} s")

    # Results container
    all_results = []

    for epsilon in epsilons:
        rows = []
        t1 = time.time()
        for name, D in CHARS:
            terms = terms_per_chi[name]
            # Build character-specific multiplier
            K_hat = build_multiplier(xi_grid, terms, D)
            # Kernel via resolvent
            K_L_values, mu_chi = build_kernel_v3(xi_grid, K_hat, epsilon, u_eval)
            # Compute A
            A = compute_A_from_kernel(terms, K_L_values, u_eval)
            rows.append({
                "chi": name, "D": D,
                "mu_chi": mu_chi,
                "A_diag": A["A_diag"], "A_off": A["A_off"],
                "A_total": A["A_total"],
                "gap_emp": emp[name],
            })
        dt = time.time() - t1
        a_off = aggregate(rows, "A_off")
        a_total = aggregate(rows, "A_total")
        print(f"[epsilon={epsilon:.2f}] off: {a_off['sign_ok_all']}/10 "
              f"R2={a_off['R2_all']:.3f} | "
              f"total: {a_total['sign_ok_all']}/10 R2={a_total['R2_all']:.3f} "
              f"| n21: off {a_off['sign_ok_n21']}/9 R2={a_off['R2_n21']:.3f} "
              f"({dt:.2f}s)")

        all_results.append({
            "epsilon": epsilon,
            "rows": rows,
            "aggregate_off": a_off,
            "aggregate_total": a_total,
        })

    # Find best
    best = max(all_results, key=lambda r: (
        r["aggregate_off"]["sign_ok_all"], r["aggregate_off"]["R2_all"]))

    print(f"\n[best epsilon] {best['epsilon']}: off {best['aggregate_off']['sign_ok_all']}/10, "
          f"R2={best['aggregate_off']['R2_all']:.3f}")

    # Print per-character detail for best
    print(f"\n[best epsilon details]")
    print(f"  {'chi':<10}{'D':>5}{'mu_chi':>12}{'A_off':>14}{'A_total':>14}{'gap_emp':>12}{'sign':>6}")
    for r in best["rows"]:
        ok = "OK" if np.sign(r["A_off"]) == np.sign(r["gap_emp"]) else "FAIL"
        print(f"  {r['chi']:<10}{r['D']:>5}{r['mu_chi']:>12.4f}"
              f"{r['A_off']:>14.4f}{r['A_total']:>14.4f}{r['gap_emp']:>12.5f}{ok:>6}")

    # Save
    out_json = RES / "PALEY_WIENER_TRANSFER_V3.json"
    with out_json.open("w", encoding="utf-8") as f:
        json.dump({
            "parameters": {
                "lambda": lam, "L": L, "eps_L": eps_L,
                "log_pm_max": log_pm_max,
                "xi_max": xi_max, "n_xi": n_xi,
                "u_max": u_max, "n_u": n_u,
                "epsilons": epsilons,
            },
            "epsilon_sweep": all_results,
            "best_epsilon": best["epsilon"],
        }, f, indent=2)
    print(f"[write] {out_json}")

    # Markdown summary
    md_path = RES / "PALEY_WIENER_TRANSFER_V3.md"
    with md_path.open("w", encoding="utf-8") as f:
        f.write("# Paley-Wiener Transfer Validation v3\n\n")
        f.write("**Date:** 2026-04-17 (Session 11)\n")
        f.write("**Script:** `_scripts/paley_wiener_transfer_v3.py`\n")
        f.write("**Approach:** Full character-specific multiplier resolvent kernel.\n\n")
        f.write("## Motivation\n\n")
        f.write("v1 (Gauss kernel): 4/10, R²=0.03. v2 (ad-hoc pole-cos): 6/10, R²=0.31. "
                "v3 uses the **physically correct** Weil multiplier "
                "$\\widehat K_\\chi(\\xi) = h_\\mathrm{arch}(\\xi) + 2\\sum w_{p,m}\\cos(m\\log p\\,\\xi)$ "
                "and takes the character-specific resolvent "
                "$R_\\chi(\\xi) = 1/(\\widehat K_\\chi(\\xi) - \\mu_\\chi + i\\epsilon)$ "
                "whose Fourier inverse is the position kernel.\n\n")
        f.write("## Parameters\n\n")
        f.write(f"- lambda={lam}, L={L:.4f}\n")
        f.write(f"- log_pm_max={log_pm_max}\n")
        f.write(f"- xi-grid: [-{xi_max}, {xi_max}], n_xi={n_xi}\n")
        f.write(f"- u-grid: [-{u_max}, {u_max}], n_u={n_u}\n")
        f.write(f"- epsilon sweep: {epsilons}\n\n")
        f.write("## Epsilon sweep results\n\n")
        f.write("| epsilon | sign(A_off) | R²(A_off) | sign(A_total) | R²(A_total) | sign_n21(A_off) | R²_n21(A_off) |\n")
        f.write("|---:|:---:|---:|:---:|---:|:---:|---:|\n")
        for r in all_results:
            a_off = r["aggregate_off"]
            a_tot = r["aggregate_total"]
            f.write(f"| {r['epsilon']:.2f} | {a_off['sign_ok_all']}/10 | "
                    f"{a_off['R2_all']:.3f} | {a_tot['sign_ok_all']}/10 | "
                    f"{a_tot['R2_all']:.3f} | {a_off['sign_ok_n21']}/9 | "
                    f"{a_off['R2_n21']:.3f} |\n")
        f.write(f"\n**Best epsilon**: {best['epsilon']}: "
                f"off {best['aggregate_off']['sign_ok_all']}/10, "
                f"R²={best['aggregate_off']['R2_all']:.3f}\n\n")
        f.write(f"## Per-character details (epsilon={best['epsilon']})\n\n")
        f.write("| chi | D | mu_chi | A_diag | A_off | A_total | gap_emp | sign |\n")
        f.write("|---|---|---:|---:|---:|---:|---:|:---:|\n")
        for r in best["rows"]:
            ok = "OK" if np.sign(r["A_off"]) == np.sign(r["gap_emp"]) else "**FAIL**"
            f.write(f"| {r['chi']} | {r['D']} | {r['mu_chi']:+.4f} | "
                    f"{r['A_diag']:+.4f} | {r['A_off']:+.4f} | "
                    f"{r['A_total']:+.4f} | {r['gap_emp']:+.5f} | {ok} |\n")
    print(f"[write] {md_path}")


if __name__ == "__main__":
    main()
