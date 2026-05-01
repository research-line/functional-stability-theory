#!/usr/bin/env python3
# coding: utf-8
"""
Dedekind NE-B analog test for SGE (Semigroup-Group Equivalence) hypothesis.

Tests whether the prime-ideal shift operators of a number field K yield
a universal commuting operator (HP-BL = YES) or only the scalar centraliser
(HP-BL = NO, analogous to Riemann NE-B).

Under the SGE hypothesis (NE-B paper §9, SEMIGROUP_GROUP_EQUIVALENCE.md):
  Zeta_K has a semigroup transfer algebra (prime ideals with norm composition)
  and therefore HP-BL should be NO (only scalar centraliser).

Test setup:
  1. Enumerate prime ideals of K = Q(sqrt(-5)) up to norm bound N_p_max
     (split / inert / ramified behaviour classified).
  2. Compute normalised shifts r_p = 2 * log(N_p) / L for a chosen L.
  3. Build NE-B shift matrices D_N(r) as in RH_II Def. (sec:NE-B), namely
        D_{nm}(r) = <phi_n, (S_{rL} + S_{-rL}) phi_m> - <psi_n, (...)  psi_m>,
     with phi_n(t) = cos(n*pi*t/L)/sqrt(L) (even),
          psi_n(t) = sin((n+1)*pi*t/L)/sqrt(L) (odd).
     Using L = 1 WLOG (L-independence lemma).
  4. Solve for the common symmetric centraliser:
        find M = M^T with [M, D_N(r_p)] = 0 for all selected p.
     Return dim of centraliser via nullspace of the corresponding linear
     system on the sym(N) = N(N+1)/2 degrees of freedom.

Outputs:
  masters/zeta-zoo/results/DEDEKIND_NE_B_TEST.json
  masters/zeta-zoo/results/DEDEKIND_NE_B_TEST.md

Three comparisons:
  (A) Q-primes (control, known NE-B = scalar centraliser)
  (B) Q(sqrt(-5))-prime-ideal norms (SGE test)
  (C) Uniform random shifts (control, should also yield scalar)
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
PROJECT_ROOT = REPO_ROOT / "masters" / "zeta-zoo"
RES = PROJECT_ROOT / "results"
RES.mkdir(exist_ok=True)


# ----------------------------------------------------------------------
# Prime ideal enumeration for K = Q(sqrt(d))
# ----------------------------------------------------------------------
def legendre_neg5(p):
    """Compute (d/p) for d = -5 via Kronecker symbol at odd p > 5."""
    if p == 2:
        return 0  # ramified at 2 for disc = -20
    if p == 5:
        return 0  # ramified
    # (-5/p) = (-1/p)(5/p)
    # (-1/p) = 1 if p ≡ 1 mod 4 else -1
    # (5/p): quadratic reciprocity, (5/p) = (p/5) since 5 ≡ 1 mod 4
    # (p/5) depends on p mod 5: 1,4 are QR; 2,3 are NR
    m1 = 1 if p % 4 == 1 else -1
    pm5 = p % 5
    m5 = 1 if pm5 in (1, 4) else -1
    return m1 * m5


def prime_ideals_Q_sqrt_neg5(norm_max):
    """Enumerate prime ideals of K = Q(sqrt(-5)) with norm <= norm_max.

    Returns list of (norm, description).
    Split:     p -> two prime ideals of norm p
    Inert:     p -> one prime ideal of norm p^2
    Ramified:  p -> one prime ideal of norm p (ramified at 2, 5)
    """
    # Sieve primes up to norm_max (for inert we need p <= sqrt(norm_max))
    from sympy import primerange
    ideals = []
    for p in primerange(2, norm_max + 1):
        chi = legendre_neg5(p)
        if chi == 0:  # ramified
            ideals.append((p, f"ram_{p}"))
        elif chi == 1:  # split
            ideals.append((p, f"split_{p}_a"))
            ideals.append((p, f"split_{p}_b"))
        else:  # inert, norm = p^2
            if p * p <= norm_max:
                ideals.append((p * p, f"inert_{p}"))
    return sorted(ideals, key=lambda x: x[0])


def rational_primes(p_max):
    """Primes of Q with p <= p_max, as (norm, desc) tuples."""
    from sympy import primerange
    return [(p, f"Q_{p}") for p in primerange(2, p_max + 1)]


# ----------------------------------------------------------------------
# NE-B shift matrices D_N(r) via numerical quadrature
# ----------------------------------------------------------------------
def phi(n, t):
    """Even basis phi_n(t) = cos(n*pi*t) on [-1, 1], L = 1."""
    if n == 0:
        return 1.0 / math.sqrt(2.0)
    return math.cos(n * math.pi * t)


def psi(n, t):
    """Odd basis psi_n(t) = sin((n+1)*pi*t)."""
    return math.sin((n + 1) * math.pi * t)


def shift_overlap(basis_fn, n, m, r):
    """Compute <f_n, (S_r + S_{-r}) f_m> on [-1, 1] (L=1).

    (S_r f)(t) = f(t - r) if |t-r| <= 1 else 0.
    """
    def integrand(t):
        val = 0.0
        if -1.0 <= t - r <= 1.0:
            val += basis_fn(m, t - r)
        if -1.0 <= t + r <= 1.0:
            val += basis_fn(m, t + r)
        return basis_fn(n, t) * val

    # Quad with breakpoints to avoid discontinuities
    breakpoints = sorted({-1.0, 1.0, -1 - r, -1 + r, 1 - r, 1 + r})
    breakpoints = [b for b in breakpoints if -1.0 <= b <= 1.0]
    if breakpoints[0] > -1.0:
        breakpoints.insert(0, -1.0)
    if breakpoints[-1] < 1.0:
        breakpoints.append(1.0)

    total = 0.0
    for i in range(len(breakpoints) - 1):
        val, _ = quad(integrand, breakpoints[i], breakpoints[i + 1],
                      limit=100, epsabs=1e-12, epsrel=1e-10)
        total += val
    return total


def build_D(N, r):
    """Build N x N NE-B matrix D_N(r) = <phi, shift phi> - <psi, shift psi>."""
    D = np.zeros((N, N))
    for n in range(N):
        for m in range(n, N):
            cos_part = shift_overlap(phi, n, m, r)
            sin_part = shift_overlap(psi, n, m, r)
            D[n, m] = cos_part - sin_part
            if n != m:
                D[m, n] = D[n, m]
    return D


# ----------------------------------------------------------------------
# Common centraliser of a family of symmetric matrices
# ----------------------------------------------------------------------
def centraliser_dim(matrices, symmetric=True, tol=1e-8):
    """Compute dim of {M : [M, D_j] = 0 for all j, M = M^T}.

    We enumerate the N(N+1)/2 degrees of freedom of a symmetric matrix,
    and require [M, D_j] = 0, which is an N^2 linear system per D_j.
    Stack all constraints and find nullspace dimension.
    """
    N = matrices[0].shape[0]
    # Build constraints. Parametrise M by upper triangle.
    # Variable index: (a, b) with a <= b maps to k = a*N - a*(a-1)//2 + (b-a)
    n_dof = N * (N + 1) // 2

    def var_index(a, b):
        if a > b:
            a, b = b, a
        return a * N - (a * (a - 1)) // 2 + (b - a)

    def sym_matrix(x):
        M = np.zeros((N, N))
        for a in range(N):
            for b in range(a, N):
                k = var_index(a, b)
                M[a, b] = x[k]
                M[b, a] = x[k]
        return M

    # Build constraint matrix A: [A x] flattened = [M, D_j] flattened
    rows = []
    for D in matrices:
        for i in range(N):
            for j in range(N):
                # (M D - D M)_{ij} = sum_k M_{ik} D_{kj} - D_{ik} M_{kj}
                row = np.zeros(n_dof)
                for k in range(N):
                    row[var_index(i, k)] += D[k, j]
                    row[var_index(k, j)] -= D[i, k]
                rows.append(row)

    A = np.array(rows)
    # Rank via SVD
    u, s, vt = np.linalg.svd(A, full_matrices=False)
    rank = int(np.sum(s > tol * s[0])) if len(s) > 0 and s[0] > 0 else 0
    dim_ker = n_dof - rank
    return dim_ker, s


# ----------------------------------------------------------------------
# Main test
# ----------------------------------------------------------------------
def run_test(label, norms, L, N_list):
    """Given a list of norms, convert to shifts r = 2 log(norm) / L,
    restrict to r in (0, 2), build D_N(r) for each, compute centraliser dim.
    """
    shifts = []
    for nrm in norms:
        r = 2.0 * math.log(nrm) / L
        if 0 < r < 2:
            shifts.append(r)
    shifts = sorted(set(shifts))  # unique
    print(f"\n=== {label} ===")
    print(f"  N norms: {len(norms)}, shifts in (0,2): {len(shifts)}")
    if len(shifts) >= 3:
        print(f"  sample shifts: {[f'{s:.4f}' for s in shifts[:5]]}"
              f" ... {[f'{s:.4f}' for s in shifts[-3:]]}"
              if len(shifts) > 8 else f"  all shifts: {[f'{s:.4f}' for s in shifts]}")

    results = {}
    for N in N_list:
        print(f"  --- N = {N} ---")
        matrices = [build_D(N, r) for r in shifts]
        dim, sigma = centraliser_dim(matrices)
        s_min = sigma[-1] if len(sigma) > 0 else None
        s_max = sigma[0] if len(sigma) > 0 else None
        results[N] = {
            "dim_centraliser": dim,
            "n_shifts": len(shifts),
            "s_min": float(s_min) if s_min is not None else None,
            "s_max": float(s_max) if s_max is not None else None,
        }
        print(f"    dim(Z) = {dim}  (expected: 1 if NE-B-analog holds)")
        print(f"    singular spectrum range: [{s_min:.2e}, {s_max:.2e}]")
    return {"shifts": shifts, "L": L, "results": results}


def main():
    # Parameters
    L = 10.0   # log(lambda) for lambda ~= 22000
    # Bound so that many primes have r = 2 log(norm)/L in (0, 2),
    # i.e. norm in (1, e^L). For L = 10: norm < e^10 ≈ 22026.
    norm_max_Q = 500       # primes of Q up to 500 (well within e^10)
    norm_max_K = 500       # prime ideals of Q(sqrt -5) with norm up to 500
    N_list = [5, 7, 10]

    print(f"[setup] L = {L}, norm_max_Q = {norm_max_Q}, norm_max_K = {norm_max_K}")
    print(f"[setup] N_list = {N_list}")

    # Control A: Q primes
    Q_primes = [n for n, _ in rational_primes(norm_max_Q)]
    res_Q = run_test("Control A: Q primes (NE-B classical)", Q_primes, L, N_list)

    # Test B: Q(sqrt -5) prime ideals
    K_ideals = prime_ideals_Q_sqrt_neg5(norm_max_K)
    K_norms = [n for n, _ in K_ideals]
    res_K = run_test("Test B: Q(sqrt -5) prime ideals (SGE prediction: NO)",
                     K_norms, L, N_list)

    # Control C: random shifts
    rng = np.random.default_rng(42)
    n_random = max(50, len(K_norms))  # comparable density
    random_shifts_raw = rng.uniform(0.01, 1.99, size=n_random).tolist()
    random_norms = [math.exp(r * L / 2.0) for r in random_shifts_raw]
    res_rand = run_test("Control C: random shifts", random_norms, L, N_list)

    # ---- Summary ----
    summary = {
        "parameters": {
            "L": L,
            "norm_max_Q": norm_max_Q,
            "norm_max_K": norm_max_K,
            "N_list": N_list,
        },
        "Q_primes": {"n_primes": len(Q_primes), **res_Q},
        "K_ideals": {
            "n_ideals": len(K_ideals),
            "n_ramified": sum(1 for _, d in K_ideals if d.startswith("ram")),
            "n_split": sum(1 for _, d in K_ideals if d.startswith("split")),
            "n_inert": sum(1 for _, d in K_ideals if d.startswith("inert")),
            **res_K,
        },
        "random": {"n_shifts": n_random, **res_rand},
    }

    out_json = RES / "DEDEKIND_NE_B_TEST.json"
    with out_json.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"\n[write] {out_json}")

    # MD report
    md = RES / "DEDEKIND_NE_B_TEST.md"
    with md.open("w", encoding="utf-8") as f:
        f.write("# Dedekind NE-B Analog Test (SGE Probe)\n\n")
        f.write("**Datum:** 2026-04-16\n")
        f.write("**Skript:** `scripts/zeta-zoo/dedekind_ne_b_test.py`\n")
        f.write("**Motivation:** Test der SGE-Hypothese an Q(sqrt(-5)) als erste "
                "Nicht-Q-Familie. SGE sagt HP-BL = NO voraus (Semigruppe der Primideale).\n\n")
        f.write("## Parameter\n\n")
        f.write(f"- L = {L}, norm_max_Q = {norm_max_Q}, norm_max_K = {norm_max_K}\n")
        f.write(f"- N_list = {N_list}\n\n")
        f.write(f"## Q(sqrt(-5)) Primideal-Struktur\n\n")
        f.write(f"- Total Primideale: {len(K_ideals)}\n")
        f.write(f"- Ramified: {summary['K_ideals']['n_ramified']} (p in {{2, 5}})\n")
        f.write(f"- Split: {summary['K_ideals']['n_split']} (p ≡ 1,3,7,9 mod 20)\n")
        f.write(f"- Inert (norm = p^2): {summary['K_ideals']['n_inert']} (p ≡ 11,13,17,19 mod 20)\n\n")
        f.write("## Ergebnisse: dim(Zentralisator)\n\n")
        f.write("| N | Q-primes (Control A) | Q(sqrt -5)-Ideale (Test B) | Random (Control C) |\n")
        f.write("|:-:|:-:|:-:|:-:|\n")
        for N in N_list:
            f.write(f"| {N} | {res_Q['results'][N]['dim_centraliser']} "
                    f"| {res_K['results'][N]['dim_centraliser']} "
                    f"| {res_rand['results'][N]['dim_centraliser']} |\n")
        f.write("\n**Erwartung unter SGE:** dim = 1 (nur skalar) in allen drei Spalten "
                "für genügend dichte Shift-Sets.\n\n")
        f.write("## Interpretation\n\n")
        dim_Q = res_Q['results'][N_list[-1]]['dim_centraliser']
        dim_K = res_K['results'][N_list[-1]]['dim_centraliser']
        if dim_Q == 1 and dim_K == 1:
            f.write("**Ergebnis: SGE-konsistent.** Der Zentralisator beider Familien ist "
                    "nur skalar (dim = 1) bei höchstem getesteten N. Q(sqrt(-5)) zeigt "
                    "dieselbe NE-B-Struktur wie Q, konsistent mit der Vorhersage "
                    "HP-BL(zeta_K) = NO unter der SGE-Hypothese (Semigruppen-Indexmenge).\n\n")
            f.write("**Erste Bestätigung der SGE-Hypothese ausserhalb der drei "
                    "kanonischen Bausteine (Riemann, Selberg, Prime-Hub).**\n")
        elif dim_K > 1:
            f.write(f"**WARNUNG: Unerwartetes Ergebnis.** dim(Z) = {dim_K} für "
                    f"Q(sqrt(-5)) bei N = {N_list[-1]}. Das wuerde SGE widerlegen "
                    "oder den Test-Setup in Frage stellen. Detailliertes Debug "
                    "erforderlich.\n")
        else:
            f.write(f"**Unvollständiges Ergebnis.** Q-primes geben dim = {dim_Q}, "
                    f"Q(sqrt(-5))-Ideale geben dim = {dim_K}. Überprüfung "
                    "der Shift-Dichte und Numerik nötig.\n")


if __name__ == "__main__":
    main()
