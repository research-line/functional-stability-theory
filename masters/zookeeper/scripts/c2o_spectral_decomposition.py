"""
c2o_spectral_decomposition.py — Spektrale Zerlegung von δ_j für C2o.1

Verifiziert die Schlüsselidentität:
  g_j^⊥ = -f_j Σ_{a∉Cl} c_a α_a / (t_j - w_a)

Und berechnet den Cauchy-Schwarz-Bound:
  |δ_j| ≤ ||k^⊥||/α · (t_j-w_min) · √(Σ_{a∉Cl} α_a²/(t_j-w_a)²)

Kernfrage: Reicht ||k^⊥|| ≤ 1 (trivial) für |δ_j| < 1,
oder brauchen wir ||k^⊥|| ~ λ^{-1.84} (= MS2)?

Falls ||k^⊥|| ≤ 1 reicht: C2o.1 folgt OHNE MS2!
Falls nicht: C2o.1 ⟺ MS2 (zirkulär).

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


def run(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    print(f"\n{'='*78}")
    print(f"  C2o SPEKTRALE ZERLEGUNG:  λ={lam}, N={N}, dim={dim}")
    print(f"{'='*78}")

    t0 = time.time()
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

    # A = QW, T = PHP
    ws_A, Qs = diagonalize_mp(Aq, verbose=False)
    w_min = float(ws_A[0])

    # f = P₀ Hũ
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
    print(f"  Build+Diag: {time.time()-t0:.0f}s", flush=True)

    # === A-Eigenbasis: c_a, α_a ===
    cl_A = [a for a in range(dim) if abs(float(ws_A[a]) - w_min) < 1e-10]
    noncl_A = [a for a in range(dim) if a not in cl_A]

    c_arr = np.zeros(dim)  # c_a = ⟨q_a, k⟩
    alpha_a = np.zeros(dim)  # α_a = ⟨ũ, q_a⟩
    w_arr = np.array([float(ws_A[a]) for a in range(dim)])

    for a in range(dim):
        qa = mp.matrix(dim, 1)
        for i in range(dim):
            qa[i, 0] = Qs[i, a]
        c_arr[a] = float(inner_product(qa, kn, dim))
        alpha_a[a] = float(inner_product(ut, qa, dim))

    alpha_cl = sum(c_arr[a] * alpha_a[a] for a in cl_A)
    k_perp_sq = sum(c_arr[a]**2 for a in noncl_A)
    k_perp_norm = math.sqrt(k_perp_sq)

    print(f"  α={alpha:.6f}, α_cl={alpha_cl:.6f}, α_cl/α={alpha_cl/alpha:.8f}")
    print(f"  ||k^⊥||={k_perp_norm:.6e}, ||k^⊥||²={k_perp_sq:.6e}")
    print(f"  Cluster: {len(cl_A)}, Non-Cluster: {len(noncl_A)}")

    # === T-Eigenbasis: f_j, g_j ===
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

    # === Identität verifizieren ===
    print(f"\n  {'='*78}")
    print(f"  IDENTITÄT: g_j^⊥ = -f_j Σ_{{a∉Cl}} c_a α_a / (t_j - w_a)")
    print(f"  {'='*78}")

    # Berechne g_j via Identität
    g_cl_arr = np.zeros(dim)
    g_perp_direct = np.zeros(dim)  # direkt: g - g^cl
    g_perp_identity = np.zeros(dim)  # via Identität

    for j in j_reg:
        d = t_vals[j] - w_min
        if abs(d) < 1e-30:
            continue
        g_cl_arr[j] = -alpha_cl * f_arr[j] / d

        g_perp_direct[j] = g_arr[j] - g_cl_arr[j]

        # Identität: g_j^⊥ = -f_j Σ_{a∉cl} c_a α_a / (t_j - w_a)
        s = 0.0
        for a in noncl_A:
            wa = w_arr[a]
            denom = t_vals[j] - wa
            if abs(denom) > 1e-30:
                s += c_arr[a] * alpha_a[a] / denom
        g_perp_identity[j] = -f_arr[j] * s

    # Vergleich
    print(f"\n  {'j':>4}  {'g^⊥ (direkt)':>14}  {'g^⊥ (Ident.)':>14}  {'Rel. Fehler':>14}")
    print(f"  {'-'*50}")
    max_rel_err = 0
    for j in j_reg[:20]:
        if abs(g_perp_direct[j]) < 1e-30:
            continue
        rel = abs(g_perp_direct[j] - g_perp_identity[j]) / abs(g_perp_direct[j])
        max_rel_err = max(max_rel_err, rel)
        print(f"  {j:4d}  {g_perp_direct[j]:14.6e}  {g_perp_identity[j]:14.6e}  {rel:14.6e}")

    print(f"\n  Max rel. Fehler: {max_rel_err:.6e}")
    print(f"  IDENTITÄT {'VERIFIZIERT' if max_rel_err < 1e-6 else 'FEHLER'}!")

    # === δ_j Zerlegung: welche a∉Cl tragen bei? ===
    print(f"\n  {'='*78}")
    print(f"  δ_j ZERLEGUNG: Welche Non-Cluster A-EVs tragen bei?")
    print(f"  {'='*78}")

    # Für jeden j∈J_reg: berechne δ_j und die Top-Beiträge
    delta_arr = np.zeros(dim)
    for j in j_reg:
        d = t_vals[j] - w_min
        if abs(d) < 1e-30:
            continue
        delta_arr[j] = d * g_perp_direct[j] / (alpha * f_arr[j])

    # Gesamt-Beitrag jedes a∉Cl zu Σ|δ_j|²
    contrib_a = np.zeros(dim)
    for a in noncl_A:
        for j in j_reg:
            d = t_vals[j] - w_min
            da = t_vals[j] - w_arr[a]
            if abs(d) < 1e-30 or abs(da) < 1e-30:
                continue
            term = (c_arr[a] * alpha_a[a] / alpha) * (d / da)
            contrib_a[a] += term**2

    # Top Non-Cluster A-EVs nach Beitrag
    noncl_sorted = sorted(noncl_A, key=lambda a: -contrib_a[a])
    print(f"\n  Top Non-Cluster A-EVs nach Beitrag zu Σ|δ_j|²:")
    print(f"  {'a':>4}  {'w_a':>10}  {'w_a-w_min':>10}  {'|c_a|':>10}  {'|α_a|':>10}  {'|c_a·α_a|':>10}  {'Beitrag':>12}")
    print(f"  {'-'*80}")
    for rank in range(min(15, len(noncl_sorted))):
        a = noncl_sorted[rank]
        print(f"  {a:4d}  {w_arr[a]:10.4f}  {w_arr[a]-w_min:10.4e}  {abs(c_arr[a]):10.4e}  "
              f"{abs(alpha_a[a]):10.4e}  {abs(c_arr[a]*alpha_a[a]):10.4e}  {contrib_a[a]:12.4e}")

    # === Cauchy-Schwarz Bound ===
    print(f"\n  {'='*78}")
    print(f"  CAUCHY-SCHWARZ BOUND: |δ_j| ≤ ||k^⊥||/α · (t_j-w) · √Q_j")
    print(f"  mit Q_j = Σ_{{a∉Cl}} α_a²/(t_j-w_a)²")
    print(f"  {'='*78}")

    # Berechne Q_j für jede Mode
    Q_arr = np.zeros(dim)
    for j in j_reg:
        for a in noncl_A:
            da = t_vals[j] - w_arr[a]
            if abs(da) > 1e-30:
                Q_arr[j] += alpha_a[a]**2 / da**2

    # Cauchy-Schwarz Bound
    cs_bound = np.zeros(dim)
    cs_bound_trivial = np.zeros(dim)  # mit ||k^⊥|| = 1
    for j in j_reg:
        d = t_vals[j] - w_min
        cs_bound[j] = (k_perp_norm / abs(alpha)) * d * math.sqrt(Q_arr[j])
        cs_bound_trivial[j] = (1.0 / abs(alpha)) * d * math.sqrt(Q_arr[j])

    print(f"\n  {'j':>4}  {'t_j-w':>10}  {'|δ_j|':>10}  {'CS(k⊥)':>10}  {'CS(1)':>10}  "
          f"{'CS(k⊥)<1':>9}  {'CS(1)<1':>7}  {'√Q_j':>10}")
    print(f"  {'-'*85}")

    n_cs_ok = 0
    n_cs_trivial_ok = 0
    for j in sorted(j_reg, key=lambda j: t_vals[j]):
        d = t_vals[j] - w_min
        if abs(d) < 1e-30:
            continue
        cs_ok = "JA" if cs_bound[j] < 1 else "NEIN"
        cs_triv = "JA" if cs_bound_trivial[j] < 1 else "NEIN"
        if cs_bound[j] < 1:
            n_cs_ok += 1
        if cs_bound_trivial[j] < 1:
            n_cs_trivial_ok += 1
        print(f"  {j:4d}  {d:10.4e}  {abs(delta_arr[j]):10.4e}  {cs_bound[j]:10.4e}  "
              f"{cs_bound_trivial[j]:10.4e}  {cs_ok:>9}  {cs_triv:>7}  {math.sqrt(Q_arr[j]):10.4e}")

    print(f"\n  === ZUSAMMENFASSUNG ===")
    print(f"  ||k^⊥|| = {k_perp_norm:.6e}")
    print(f"  CS-Bound mit ||k^⊥||: {n_cs_ok}/{len(j_reg)} Moden erfüllt")
    print(f"  CS-Bound mit ||k^⊥||=1 (trivial): {n_cs_trivial_ok}/{len(j_reg)} Moden erfüllt")

    if n_cs_trivial_ok == len(j_reg):
        print(f"\n  *** DURCHBRUCH: CS-Bound mit ||k^⊥||≤1 REICHT! ***")
        print(f"  *** C2o.1 folgt OHNE MS2! ***")
    elif n_cs_ok == len(j_reg):
        print(f"\n  CS-Bound mit tatsächlichem ||k^⊥|| reicht.")
        print(f"  Braucht ||k^⊥|| < {1/max(cs_bound_trivial[j]/max(cs_bound[j],1e-30) for j in j_reg if cs_bound[j]>0):.2e}")
    else:
        print(f"\n  CS-Bound reicht NICHT — braucht schärferen Bound.")

    # Tightness: wie eng ist CS?
    tightness = []
    for j in j_reg:
        if cs_bound[j] > 0:
            tightness.append(abs(delta_arr[j]) / cs_bound[j])
    if tightness:
        print(f"\n  Tightness |δ_j|/CS-Bound:")
        print(f"    max  = {max(tightness):.6e}")
        print(f"    mean = {np.mean(tightness):.6e}")
        print(f"    min  = {min(tightness):.6e}")

    # === Schärferer Bound: Sign-Coherence ausnutzen ===
    print(f"\n  {'='*78}")
    print(f"  SIGN-COHERENCE: c_a α_a Vorzeichen-Analyse")
    print(f"  {'='*78}")

    pos_ca = sum(1 for a in noncl_A if c_arr[a] * alpha_a[a] > 0)
    neg_ca = len(noncl_A) - pos_ca
    print(f"  c_a · α_a > 0: {pos_ca}/{len(noncl_A)}")
    print(f"  c_a · α_a < 0: {neg_ca}/{len(noncl_A)}")
    sum_pos = sum(c_arr[a] * alpha_a[a] for a in noncl_A if c_arr[a] * alpha_a[a] > 0)
    sum_neg = sum(c_arr[a] * alpha_a[a] for a in noncl_A if c_arr[a] * alpha_a[a] < 0)
    print(f"  Σ (c_a α_a, pos) = {sum_pos:.6e}")
    print(f"  Σ (c_a α_a, neg) = {sum_neg:.6e}")
    print(f"  Σ c_a α_a (total) = {sum_pos + sum_neg:.6e}")
    print(f"  Cancellation: {abs(sum_pos+sum_neg) / (abs(sum_pos)+abs(sum_neg)):.6e}" if abs(sum_pos)+abs(sum_neg) > 0 else "")

    return {
        "lam": lam,
        "k_perp": k_perp_norm,
        "max_delta": max(abs(delta_arr[j]) for j in j_reg),
        "max_cs_bound": max(cs_bound[j] for j in j_reg),
        "max_cs_trivial": max(cs_bound_trivial[j] for j in j_reg),
        "n_cs_ok": n_cs_ok,
        "n_cs_trivial_ok": n_cs_trivial_ok,
        "n_reg": len(j_reg),
        "identity_err": max_rel_err,
    }


if __name__ == "__main__":
    results = []
    for cfg in CONFIGS:
        r = run(cfg["lam"], cfg["N"])
        results.append(r)

    if len(results) >= 2:
        print(f"\n{'='*78}")
        print(f"  λ-VERGLEICH (Spektrale Zerlegung)")
        print(f"{'='*78}")
        for name, key in [
            ("||k^⊥||", "k_perp"),
            ("max |δ_j|", "max_delta"),
            ("max CS(k⊥)", "max_cs_bound"),
            ("max CS(1)", "max_cs_trivial"),
            ("CS(k⊥) OK", "n_cs_ok"),
            ("CS(1) OK", "n_cs_trivial_ok"),
            ("Identität err", "identity_err"),
        ]:
            vals = [r[key] for r in results]
            print(f"  {name:>20}  ", end="")
            for v in vals:
                print(f"{v:14.4e}  ", end="")
            if len(vals) >= 2 and isinstance(vals[0], float) and abs(vals[0]) > 0 and abs(vals[-1]) > 0:
                lams = [r["lam"] for r in results]
                try:
                    exp = math.log(abs(vals[-1]) / abs(vals[0])) / math.log(lams[-1] / lams[0])
                    print(f"  {exp:8.2f}")
                except (ValueError, ZeroDivisionError):
                    print()
            else:
                print()

    print(f"\n=== FERTIG ===")
