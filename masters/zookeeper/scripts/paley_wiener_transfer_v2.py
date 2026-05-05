#!/usr/bin/env python3
# coding: utf-8
"""
Session 8 Iteration 2 - Paley-Wiener Transfer Validation with POLE-AWARE kernel.

Previous result (Gaussian-only kernel, v1): sign_off = 4/10, R^2 = 0.26.
Diagnosis: archimedean-centred Gaussian of width sigma_Phi ~ 1.24 concentrates
Fourier support at xi=0, while character-specific zeros gamma_chi^(k) >= 3
contribute only via exponentially suppressed Fourier tails
(e^{-gamma^2 sigma^2} ~ e^{-21} for gamma=3, sigma=1.24).

This iteration tests two modifications:

  (A) SCAN over sigma width: the archimedean width is a leading-order
      estimate; the true perturbed ground state is broader due to pole
      contributions. Scan sigma in {1.2, 2.0, 3.0, 4.0, 5.0}.

  (B) POLE-AWARE kernel with explicit zero-pole oscillations:
      K_L(u; chi) = (1/sqrt(4*pi*sigma^2)) exp(-u^2/(4 sigma^2))
                    + alpha * sum_k weight_k * cos(gamma_chi^(k) * u)
      with weight_k = 1/gamma_k (Hadamard-style) or = exp(-gamma_k^2 sigma^2)
      (Fourier-suppressed). Alpha and weighting are both tested.

For each variant, we compute A_chi^{off,total} and compare signs with
the empirical gap_galerkin at N=600, lambda=20000.

Output: _results/PALEY_WIENER_TRANSFER_V2.{json,md}
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
# Kronecker symbol (same as v1)
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


def load_zeros():
    path = RES / "zeros_all_chars.json"
    with path.open(encoding="utf-8") as f:
        return json.load(f)


# --------------------------------------------------------------------------
# Prime term collection (width-aware)
# --------------------------------------------------------------------------
def collect_prime_terms(chi_fn, D, log_pm_max, max_m=8):
    """Collect (chi^m, m*log p, log p / p^{m/2}) with m*log(p) <= log_pm_max."""
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
# Kernel variants
# --------------------------------------------------------------------------
def K_gauss(u_matrix, sigma):
    """Pure Gaussian kernel (v1 default)."""
    return np.exp(-(u_matrix ** 2) / (4.0 * sigma * sigma)) / math.sqrt(
        4.0 * math.pi * sigma * sigma)


def K_pole_cos(u_matrix, zeros, weighting="hadamard", sigma=None):
    """Pure pole-oscillation kernel: sum_k w_k * cos(gamma_k * u)."""
    out = np.zeros_like(u_matrix)
    for gamma in zeros:
        if weighting == "hadamard":
            w = 1.0 / gamma
        elif weighting == "hadamard2":
            w = 1.0 / (gamma * gamma)
        elif weighting == "fourier":
            # exp(-gamma^2 * sigma^2) factor (Fourier-suppressed)
            w = math.exp(-gamma * gamma * sigma * sigma)
        else:
            w = 1.0
        out += w * np.cos(gamma * u_matrix)
    return out


def K_mixed(u_matrix, sigma, zeros, alpha, weighting="hadamard"):
    """Mixed: Gaussian + alpha * pole oscillation."""
    return K_gauss(u_matrix, sigma) + alpha * K_pole_cos(
        u_matrix, zeros, weighting=weighting, sigma=sigma)


# --------------------------------------------------------------------------
# Compute A_chi for given kernel
# --------------------------------------------------------------------------
def compute_A(terms, K_matrix):
    if not terms:
        return {"n_terms": 0, "A_diag": 0.0, "A_off": 0.0, "A_total": 0.0}
    c_arr = np.array([t[0] for t in terms], dtype=float)
    w_arr = np.array([t[2] for t in terms], dtype=float)
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
    # Without chi_21
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
# Run one variant
# --------------------------------------------------------------------------
def run_variant(label, kernel_fn, terms_per_chi, emp, zeros_dict):
    """kernel_fn(u_matrix, chi_name) -> K_matrix."""
    rows = []
    for name, D in CHARS:
        terms = terms_per_chi[name]
        log_arr = np.array([t[1] for t in terms], dtype=float)
        u = log_arr[:, None] - log_arr[None, :]
        K = kernel_fn(u, name)
        A = compute_A(terms, K)
        rows.append({
            "chi": name, "D": D,
            "A_diag": A["A_diag"], "A_off": A["A_off"], "A_total": A["A_total"],
            "gap_emp": emp[name],
        })
    return {
        "label": label,
        "rows": rows,
        "aggregate_off": aggregate(rows, "A_off"),
        "aggregate_total": aggregate(rows, "A_total"),
    }


# --------------------------------------------------------------------------
# Main: sweep variants
# --------------------------------------------------------------------------
def main():
    lam = 20000
    L = math.log(lam)
    eps_L = 1.0 / math.sqrt(L)
    print(f"[setup] lambda={lam}, L={L:.4f}, eps_L={eps_L:.4f}")

    emp = load_empirical_gaps()
    zeros_dict = load_zeros()

    # Kernel support cutoff: primes with log(p^m) <= log_pm_max are included.
    # Practical choice: 9.0 (p^m <= e^9 ~ 8100, ~1000 primes). This is wide
    # enough to capture the dominant contributions (since weights decay as
    # 1/p^{m/2}) and small enough for pairwise matrix to fit in RAM:
    # 1000^2 = 1e6 floats = 8 MB per variant.
    log_pm_max = 9.0
    print(f"[setup] log_pm_max = {log_pm_max} (p^m <= {math.exp(log_pm_max):.0f})")

    # Precompute prime terms once per character (same for all kernel variants)
    t0 = time.time()
    terms_per_chi = {}
    for name, D in CHARS:
        chi = make_chi_D(D)
        terms = collect_prime_terms(chi, D, log_pm_max, max_m=8)
        terms_per_chi[name] = terms
        print(f"  [{name}] D={D}, n_terms={len(terms)}")
    print(f"[prime collection] {time.time()-t0:.2f} s")

    # --- Variant sweep ---
    variants = []

    # (A) Pure Gaussian with varying sigma
    for sigma in [1.2, 2.0, 3.0, 4.0, 5.0]:
        label = f"gauss_sigma{sigma}"
        def kf(u, name, s=sigma):
            return K_gauss(u, s)
        variants.append((label, kf))

    # (B) Pure pole-cos with different weighting
    for weighting in ["hadamard", "hadamard2", "equal"]:
        label = f"pole_{weighting}"
        def kf(u, name, w=weighting):
            return K_pole_cos(u, zeros_dict[name], weighting=w)
        variants.append((label, kf))

    # (C) Mixed Gaussian + pole
    for sigma in [2.0, 3.0]:
        for alpha in [0.1, 0.5, 1.0, 3.0]:
            for weighting in ["hadamard", "hadamard2"]:
                label = f"mix_sigma{sigma}_alpha{alpha}_{weighting}"
                def kf(u, name, s=sigma, a=alpha, w=weighting):
                    return K_mixed(u, s, zeros_dict[name], a, weighting=w)
                variants.append((label, kf))

    # Run all
    results = []
    t0 = time.time()
    for label, kf in variants:
        t1 = time.time()
        res = run_variant(label, kf, terms_per_chi, emp, zeros_dict)
        dt = time.time() - t1
        a_off = res["aggregate_off"]
        a_tot = res["aggregate_total"]
        print(f"[variant] {label:40s}  "
              f"off: {a_off['sign_ok_all']:2d}/10 R2={a_off['R2_all']:.3f}  "
              f"total: {a_tot['sign_ok_all']:2d}/10 R2={a_tot['R2_all']:.3f}  "
              f"{dt:.2f}s")
        results.append(res)
    total = time.time() - t0
    print(f"[all variants done] {total:.2f} s, {len(variants)} variants")

    # Find best variant
    best_off = max(results, key=lambda r: (
        r["aggregate_off"]["sign_ok_all"], r["aggregate_off"]["R2_all"]))
    best_total = max(results, key=lambda r: (
        r["aggregate_total"]["sign_ok_all"], r["aggregate_total"]["R2_all"]))

    print(f"\n[best A_off] {best_off['label']} -> "
          f"{best_off['aggregate_off']['sign_ok_all']}/10, "
          f"R2={best_off['aggregate_off']['R2_all']:.3f}")
    print(f"[best A_total] {best_total['label']} -> "
          f"{best_total['aggregate_total']['sign_ok_all']}/10, "
          f"R2={best_total['aggregate_total']['R2_all']:.3f}")

    # Save
    out_json = RES / "PALEY_WIENER_TRANSFER_V2.json"
    with out_json.open("w", encoding="utf-8") as f:
        json.dump({
            "parameters": {"lambda": lam, "L": L, "eps_L": eps_L,
                           "log_pm_max": log_pm_max, "total_time_s": total},
            "variants": results,
            "best_off": best_off["label"],
            "best_total": best_total["label"],
        }, f, indent=2)
    print(f"[write] {out_json}")

    # MD summary
    md_path = RES / "PALEY_WIENER_TRANSFER_V2.md"
    with md_path.open("w", encoding="utf-8") as f:
        f.write("# Paley-Wiener Transfer Validation v2 (Session 8 Iteration 2)\n\n")
        f.write(f"**Date:** 2026-04-16\n")
        f.write(f"**Script:** `_scripts/paley_wiener_transfer_v2.py`\n")
        f.write(f"**Approach:** Scan over kernel variants (Gaussian, pole, mixed)\n\n")
        f.write(f"## Parameters\n")
        f.write(f"- lambda={lam}, L={L:.4f}, eps_L={eps_L:.4f}\n")
        f.write(f"- log_pm_max={log_pm_max} (primes cover up to exp({log_pm_max}))\n")
        f.write(f"- Total walltime: {total:.2f} s, {len(variants)} variants\n\n")
        f.write(f"## Variant scan results\n\n")
        f.write("| Variant | sign(A_off) | R2(A_off) | sign(A_total) | R2(A_total) |\n")
        f.write("|---|:---:|---:|:---:|---:|\n")
        for res in results:
            a_off = res["aggregate_off"]
            a_tot = res["aggregate_total"]
            f.write(f"| `{res['label']}` | "
                    f"{a_off['sign_ok_all']}/10 | {a_off['R2_all']:.3f} | "
                    f"{a_tot['sign_ok_all']}/10 | {a_tot['R2_all']:.3f} |\n")
        f.write("\n## Best variants\n\n")
        f.write(f"- **Best A_off:** `{best_off['label']}` - "
                f"{best_off['aggregate_off']['sign_ok_all']}/10, "
                f"R²={best_off['aggregate_off']['R2_all']:.3f}\n")
        f.write(f"- **Best A_total:** `{best_total['label']}` - "
                f"{best_total['aggregate_total']['sign_ok_all']}/10, "
                f"R²={best_total['aggregate_total']['R2_all']:.3f}\n\n")
        # Detailed best-off table
        f.write(f"## Per-character results for best A_off: `{best_off['label']}`\n\n")
        f.write("| chi | D | A_diag | A_off | A_total | gap_emp | sign_off_ok |\n")
        f.write("|---|---|---:|---:|---:|---:|:---:|\n")
        for r in best_off["rows"]:
            ok = "OK" if np.sign(r["A_off"]) == np.sign(r["gap_emp"]) else "**FAIL**"
            f.write(f"| {r['chi']} | {r['D']} | {r['A_diag']:+.4f} | "
                    f"{r['A_off']:+.4f} | {r['A_total']:+.4f} | "
                    f"{r['gap_emp']:+.5f} | {ok} |\n")
    print(f"[write] {md_path}")


if __name__ == "__main__":
    main()
