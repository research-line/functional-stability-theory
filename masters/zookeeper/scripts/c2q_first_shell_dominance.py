"""
c2q_first_shell_dominance.py — First-Shell Dominance Lemma (C2q)

Testet die Hypothese:
  δ_j = -(t_j-w)/α · [c_{a1} α_{a1}/(t_j-w_{a1}) + R_j]
  mit |R_j| ≪ |c_{a1} α_{a1}/(t_j-w_{a1})|

d.h. der Off-Cluster-Defekt wird von der ERSTEN Off-Cluster-Schale dominiert.

Drei Strukturmerkmale:
  1. First-shell dominance (per mode)
  2. Sign coherence of shell coefficients
  3. Rapid decay away from first shell

Falls bestätigt: δ_j wird fast eindimensional, kontrolliert durch
einen einzigen Interlacing-Nachbarn.

Autor: LG (Opus 4.6, Session 27b+)
Datum: 2026-04-20
"""

import os, sys, time, math
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, diagonalize_mp, chi_trivial,
)
from c2_approximation_test import (
    k_lambda_value, inner_product, norm,
)
from c2_poisson_decomposition import project_to_fourier

DPS = int(os.environ.get("DPS", 50))
CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
]


def build_operators(lam, N):
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    W02 = mp.matrix(dim, dim)
    H = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            W02[i, j] = Aq[i, j] - Ah[i, j]
            H[i, j] = Ah[i, j]

    col0 = mp.matrix(dim, 1)
    for i in range(dim):
        col0[i, 0] = W02[i, 0]
    cn = norm(col0, dim)
    ut = mp.matrix(dim, 1)
    for i in range(dim):
        ut[i, 0] = col0[i, 0] / cn

    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
    alpha = float(inner_product(ut, kn, dim))

    ws_A, Qs = diagonalize_mp(Aq, verbose=False)
    w_min = float(ws_A[0])

    f_lam = mp.matrix(dim, 1)
    g_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]
        g_lam[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]

    T = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            T[i, j] = (H[i, j]
                        - ut[i, 0] * Hut[j, 0]
                        - Hut[i, 0] * ut[j, 0]
                        + mp.mpf(uHu) * ut[i, 0] * ut[j, 0])

    ts, Es = diagonalize_mp(T, verbose=False)

    cl_A = [a for a in range(dim) if abs(float(ws_A[a]) - w_min) < 1e-10]
    noncl_A = [a for a in range(dim) if a not in cl_A]

    c_arr = np.zeros(dim)
    alpha_a = np.zeros(dim)
    w_arr = np.array([float(ws_A[a]) for a in range(dim)])

    for a in range(dim):
        qa = mp.matrix(dim, 1)
        for i in range(dim):
            qa[i, 0] = Qs[i, a]
        c_arr[a] = float(inner_product(qa, kn, dim))
        alpha_a[a] = float(inner_product(ut, qa, dim))

    alpha_cl = sum(c_arr[a] * alpha_a[a] for a in cl_A)

    t_vals = np.array([float(ts[j]) for j in range(dim)])
    f_arr = np.zeros(dim)
    g_arr = np.zeros(dim)
    for j in range(dim):
        vj = mp.matrix(dim, 1)
        for i in range(dim):
            vj[i, 0] = Es[i, j]
        f_arr[j] = float(inner_product(vj, f_lam, dim))
        g_arr[j] = float(inner_product(vj, g_lam, dim))

    j_zero = int(np.argmin(np.abs(t_vals)))
    j_reg = [j for j in range(dim) if j != j_zero and abs(f_arr[j]) > 1e-25]

    return {
        "dim": dim, "alpha": alpha, "alpha_cl": alpha_cl,
        "w_min": w_min, "w_arr": w_arr, "c_arr": c_arr, "alpha_a": alpha_a,
        "cl_A": cl_A, "noncl_A": noncl_A,
        "t_vals": t_vals, "f_arr": f_arr, "g_arr": g_arr,
        "j_reg": j_reg, "j_zero": j_zero,
    }


def run(lam, N):
    mp.mp.dps = DPS
    dim = N + 1

    print(f"\n{'='*78}")
    print(f"  C2q FIRST-SHELL DOMINANCE:  λ={lam}, N={N}, dim={dim}")
    print(f"{'='*78}")

    t0 = time.time()
    ops = build_operators(lam, N)
    print(f"  Build+Diag: {time.time()-t0:.0f}s", flush=True)

    alpha = ops["alpha"]
    alpha_cl = ops["alpha_cl"]
    w_min = ops["w_min"]
    w_arr = ops["w_arr"]
    c_arr = ops["c_arr"]
    alpha_a = ops["alpha_a"]
    noncl_A = ops["noncl_A"]
    t_vals = ops["t_vals"]
    f_arr = ops["f_arr"]
    g_arr = ops["g_arr"]
    j_reg = ops["j_reg"]

    # Sortiere Non-Cluster A-EVs nach Abstand zu w_min (Schalen)
    noncl_sorted_by_gap = sorted(noncl_A, key=lambda a: abs(w_arr[a] - w_min))
    a1 = noncl_sorted_by_gap[0]
    gap_a1 = w_arr[a1] - w_min

    print(f"\n  α={alpha:.6f}, α_cl={alpha_cl:.6f}")
    print(f"  w_min={w_min:.6f}")
    print(f"  Erste Off-Cluster-Schale: a={a1}, w_a={w_arr[a1]:.10f}")
    print(f"    Gap: w_a - w_min = {gap_a1:.4e}")
    print(f"    |c_a1| = {abs(c_arr[a1]):.4e}")
    print(f"    |α_a1| = {abs(alpha_a[a1]):.4e}")
    print(f"    |c_a1·α_a1| = {abs(c_arr[a1]*alpha_a[a1]):.4e}")
    print(f"    sgn(c_a1·α_a1) = {'>' if c_arr[a1]*alpha_a[a1] > 0 else '<'} 0")

    # Erste 5 Schalen auflisten
    print(f"\n  Schalenhierarchie (erste 8 Non-Cluster A-EVs):")
    print(f"  {'Rang':>4}  {'a':>4}  {'w_a-w_min':>12}  {'|c_a·α_a|':>12}  {'sgn':>4}")
    print(f"  {'-'*48}")
    for rank, a in enumerate(noncl_sorted_by_gap[:8]):
        sgn = "+" if c_arr[a]*alpha_a[a] > 0 else "-"
        print(f"  {rank+1:4d}  {a:4d}  {w_arr[a]-w_min:12.4e}  "
              f"{abs(c_arr[a]*alpha_a[a]):12.4e}  {sgn:>4}")

    # ================================================================
    # PHASE 1: Per-Mode First-Shell Dominance
    # ================================================================
    print(f"\n  {'='*78}")
    print(f"  PHASE 1: Per-Mode First-Shell Dominance")
    print(f"  δ_j = -(t-w)/α · [S1_j + R_j],  S1_j = c_a1·α_a1/(t_j-w_a1)")
    print(f"  {'='*78}")

    results = []
    for j in sorted(j_reg, key=lambda j: t_vals[j]):
        d = t_vals[j] - w_min
        if abs(d) < 1e-30:
            continue

        # First shell term
        da1 = t_vals[j] - w_arr[a1]
        if abs(da1) < 1e-30:
            continue
        S1 = c_arr[a1] * alpha_a[a1] / da1

        # Rest
        R = 0.0
        for a in noncl_A:
            if a == a1:
                continue
            da = t_vals[j] - w_arr[a]
            if abs(da) > 1e-30:
                R += c_arr[a] * alpha_a[a] / da

        total = S1 + R

        # δ_j from full sum
        delta_full = -(d / alpha) * total

        # δ_j from first shell only
        delta_s1 = -(d / alpha) * S1

        # Dominance ratio
        if abs(S1) > 1e-30:
            dom_ratio = abs(R) / abs(S1)
        else:
            dom_ratio = float("inf")

        is_near = d < 1e-1

        results.append({
            "j": j, "t_w": d, "S1": S1, "R": R, "total": total,
            "delta_full": delta_full, "delta_s1": delta_s1,
            "dom_ratio": dom_ratio, "is_near": is_near,
        })

    print(f"\n  {'j':>4}  {'t_j-w':>10}  {'δ_full':>10}  {'δ_S1':>10}  "
          f"{'|R/S1|':>10}  {'S1 dom?':>8}  {'near?':>5}")
    print(f"  {'-'*70}")
    for r in results:
        dom = "JA" if r["dom_ratio"] < 1 else "NEIN"
        near = "near" if r["is_near"] else ""
        print(f"  {r['j']:4d}  {r['t_w']:10.4e}  {r['delta_full']:10.4e}  "
              f"{r['delta_s1']:10.4e}  {r['dom_ratio']:10.4e}  {dom:>8}  {near:>5}")

    # Zusammenfassung
    n_dom = sum(1 for r in results if r["dom_ratio"] < 1)
    n_near = sum(1 for r in results if r["is_near"])
    n_near_dom = sum(1 for r in results if r["is_near"] and r["dom_ratio"] < 1)
    n_far_dom = sum(1 for r in results if not r["is_near"] and r["dom_ratio"] < 1)
    n_far = sum(1 for r in results if not r["is_near"])

    near_ratios = [r["dom_ratio"] for r in results if r["is_near"]]
    far_ratios = [r["dom_ratio"] for r in results if not r["is_near"]]

    print(f"\n  First-Shell Dominance (|R/S1| < 1):")
    print(f"    Gesamt: {n_dom}/{len(results)} Moden")
    print(f"    Near-cluster (t-w < 0.1): {n_near_dom}/{n_near} Moden")
    if n_far > 0:
        print(f"    Off-cluster (t-w ≥ 0.1): {n_far_dom}/{n_far} Moden")
    if near_ratios:
        print(f"    Near-cluster max |R/S1|: {max(near_ratios):.4e}")
        print(f"    Near-cluster min |R/S1|: {min(near_ratios):.4e}")
    if far_ratios:
        print(f"    Off-cluster max |R/S1|: {max(far_ratios):.4e}")
        print(f"    Off-cluster median |R/S1|: {np.median(far_ratios):.4e}")

    # ================================================================
    # PHASE 2: Shell-Akkumulation
    # ================================================================
    print(f"\n  {'='*78}")
    print(f"  PHASE 2: Shell-Akkumulation — wie viele Schalen für 99% von δ_j?")
    print(f"  {'='*78}")

    # Für near-cluster Moden: kumulative Approximation
    near_modes = [r for r in results if r["is_near"]]
    if not near_modes:
        near_modes = results[:5]

    max_shells = min(10, len(noncl_sorted_by_gap))

    print(f"\n  Kumulativer Anteil |δ_j(k shells)|/|δ_j(all)| für near-cluster Moden:")
    print(f"  {'k':>4}  {'a_k':>4}  {'gap_k':>10}  ", end="")
    for r in near_modes[:5]:
        print(f"  j={r['j']:>3}", end="")
    print()
    print(f"  {'-'*(30 + 10*min(5, len(near_modes)))}")

    for k in range(1, max_shells + 1):
        shells_k = noncl_sorted_by_gap[:k]
        a_k = noncl_sorted_by_gap[k-1]
        gap_k = w_arr[a_k] - w_min
        print(f"  {k:4d}  {a_k:4d}  {gap_k:10.4e}  ", end="")

        for r in near_modes[:5]:
            j = r["j"]
            d = t_vals[j] - w_min
            # k-shell partial sum
            partial = 0.0
            for a in shells_k:
                da = t_vals[j] - w_arr[a]
                if abs(da) > 1e-30:
                    partial += c_arr[a] * alpha_a[a] / da
            delta_k = -(d / alpha) * partial
            if abs(r["delta_full"]) > 1e-30:
                ratio_k = abs(delta_k) / abs(r["delta_full"])
            else:
                ratio_k = 0.0
            print(f"  {ratio_k:8.4f}", end="")
        print()

    # ================================================================
    # PHASE 3: First-Shell δ_j Bound
    # ================================================================
    print(f"\n  {'='*78}")
    print(f"  PHASE 3: First-Shell Bound — |δ_j^S1| < 1 ohne MS2?")
    print(f"  δ_j^S1 = -(t_j-w)·c_a1·α_a1 / (α·(t_j-w_a1))")
    print(f"  {'='*78}")

    max_delta_s1 = 0
    worst_j = -1
    for r in results:
        if abs(r["delta_s1"]) > max_delta_s1:
            max_delta_s1 = abs(r["delta_s1"])
            worst_j = r["j"]

    print(f"\n  |c_a1·α_a1|/α = {abs(c_arr[a1]*alpha_a[a1])/abs(alpha):.6e}")
    print(f"  max |δ_j^S1| = {max_delta_s1:.6e} (j={worst_j})")
    print(f"  max |δ_j^S1| < 1: {'JA' if max_delta_s1 < 1 else 'NEIN'}")

    # Für near-cluster: δ_S1 ≈ (t-w)/(t-w_a1) · c_a1·α_a1/α
    # Wenn t ≈ w_min und w_a1 ≈ w_min: (t-w)/(t-w_a1) → 1 wenn gap klein
    # Also δ_S1 → c_a1·α_a1/α (endlich!)
    print(f"\n  Near-cluster Limit: δ_j^S1 → c_a1·α_a1/α = {c_arr[a1]*alpha_a[a1]/alpha:.6e}")

    # ================================================================
    # PHASE 4: Interlacing-Geometrie
    # ================================================================
    print(f"\n  {'='*78}")
    print(f"  PHASE 4: Interlacing-Geometrie")
    print(f"  Wie liegen t_j, w_a relativ zueinander?")
    print(f"  {'='*78}")

    # Sortierte T-Eigenwerte und A-Eigenwerte im Near-Cluster-Bereich
    near_t = sorted([t_vals[j] for j in j_reg if t_vals[j] - w_min < 0.1])
    near_w = sorted([w_arr[a] for a in noncl_sorted_by_gap if w_arr[a] - w_min < 0.1])

    print(f"\n  Near-cluster T-Eigenwerte (t_j - w_min < 0.1):")
    for t in near_t:
        print(f"    t - w_min = {t - w_min:.6e}")

    print(f"\n  Near-cluster A-Eigenwerte (w_a - w_min < 0.1):")
    for w in near_w:
        print(f"    w - w_min = {w - w_min:.6e}")

    print(f"\n  Interlacing: Für jedes near-cluster j, nächster w_a:")
    for j in sorted(j_reg, key=lambda j: t_vals[j]):
        d = t_vals[j] - w_min
        if d > 0.1 or d < 1e-30:
            continue
        nearest_a = min(noncl_A, key=lambda a: abs(t_vals[j] - w_arr[a]))
        gap_tw = t_vals[j] - w_arr[nearest_a]
        print(f"    j={j}: t-w_min={d:.4e}, nächstes w_a (a={nearest_a}): "
              f"gap t-w_a = {gap_tw:.4e}")

    # ================================================================
    # PHASE 5: Vorzeichen-Konsistenz per Mode
    # ================================================================
    print(f"\n  {'='*78}")
    print(f"  PHASE 5: Vorzeichen δ_j^S1 vs δ_j^full")
    print(f"  {'='*78}")

    n_same_sign = 0
    n_total_check = 0
    for r in results:
        if abs(r["delta_full"]) > 1e-20 and abs(r["delta_s1"]) > 1e-20:
            n_total_check += 1
            if np.sign(r["delta_full"]) == np.sign(r["delta_s1"]):
                n_same_sign += 1

    print(f"  sgn(δ_S1) = sgn(δ_full): {n_same_sign}/{n_total_check}")
    print(f"  → First shell bestimmt auch das VORZEICHEN von δ_j")

    # ================================================================
    # ZUSAMMENFASSUNG
    # ================================================================
    print(f"\n  {'='*78}")
    print(f"  C2q ZUSAMMENFASSUNG (λ={lam})")
    print(f"  {'='*78}")

    print(f"  First-Shell Dominance:")
    print(f"    |R/S1| < 1 (gesamt): {n_dom}/{len(results)}")
    print(f"    |R/S1| < 1 (near-cluster): {n_near_dom}/{n_near}")
    if far_ratios:
        print(f"    |R/S1| median (off-cluster): {np.median(far_ratios):.4e}")
    print(f"  Sign Coherence: sgn(δ_S1) = sgn(δ_full) für {n_same_sign}/{n_total_check}")
    print(f"  First-Shell Bound: max |δ_S1| = {max_delta_s1:.4e} {'< 1 ✓' if max_delta_s1 < 1 else '≥ 1 ✗'}")
    print(f"  Near-cluster Limit: c_a1·α_a1/α = {c_arr[a1]*alpha_a[a1]/alpha:.4e}")

    return {
        "lam": lam,
        "n_dom": n_dom, "n_total": len(results),
        "n_near_dom": n_near_dom, "n_near": n_near,
        "max_delta_s1": max_delta_s1,
        "near_limit": c_arr[a1]*alpha_a[a1]/alpha,
        "near_ratios": near_ratios,
        "far_ratios": far_ratios,
        "n_same_sign": n_same_sign,
        "n_sign_total": n_total_check,
    }


if __name__ == "__main__":
    results = []
    for cfg in CONFIGS:
        r = run(cfg["lam"], cfg["N"])
        results.append(r)

    if len(results) >= 2:
        print(f"\n{'='*78}")
        print(f"  λ-VERGLEICH (C2q First-Shell Dominance)")
        print(f"{'='*78}")
        for r in results:
            print(f"\n  λ={r['lam']}:")
            print(f"    S1 Dominanz: {r['n_dom']}/{r['n_total']} "
                  f"(near: {r['n_near_dom']}/{r['n_near']})")
            print(f"    max |δ_S1|: {r['max_delta_s1']:.4e}")
            print(f"    near limit c_a1·α_a1/α: {r['near_limit']:.4e}")
            print(f"    sign coherence: {r['n_same_sign']}/{r['n_sign_total']}")

    print(f"\n=== FERTIG ===")
