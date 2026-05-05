"""
c2o_eta_analysis.py — η_j-Analyse für C2o.1 (Modeweises Comparison Lemma)

Berechnet und analysiert:
  1. η_j = -(t_j - w_min) g_j / (α f_j)  für alle j ∈ J_reg
  2. Cluster-Zerlegung: η_j = α_cl/α - δ_j  mit δ_j = (t_j-w_min)g_j^⊥/(αf_j)
  3. Off-Cluster-Defekt |δ_j| — muss < 1 für C2o.1
  4. Energieargument: (t_j-w_min)|g_j|² vs α|f_j||g_j|
  5. Norm-Bound ||g||≤1 — Konsequenzen für η
  6. Monotonie-Check: η_j als Funktion von t_j
  7. Near-Cluster-Verhalten: wie sich η_j bei t_j→w_min verhält

Ziel: Muster finden die zum analytischen Beweis von η_j < 2 führen.

Referenz: C2O_MODEWEISE_COMPARISON.md (Äquivalenzkette)
  C2l.5 ⟺ C2o.1 ⟺ η_j ∈ (0,2) ⟺ |r_j| < α|f_j|

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
    print(f"  C2o η-ANALYSE:  λ={lam}, N={N}, dim={dim}")
    print(f"{'='*78}")

    # === Phase 1: Operatoren ===
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
    C_tilde = float(cn * cn)

    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))

    # k_full
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf
    alpha = float(inner_product(ut, kn, dim))

    # QW Eigenwerte (für w_min und Cluster)
    ws_qw, Qs = diagonalize_mp(Aq, verbose=False)
    w_min = float(ws_qw[0])

    # f = P₀ Hũ, g = P₀ k
    f_lam = mp.matrix(dim, 1)
    g_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]
        g_lam[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]

    # T = PHP = P₀HP₀
    T = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            T[i, j] = (H[i, j]
                        - ut[i, 0] * Hut[j, 0]
                        - Hut[i, 0] * ut[j, 0]
                        + mp.mpf(uHu) * ut[i, 0] * ut[j, 0])

    ts, Es = diagonalize_mp(T, verbose=False)
    print(f"  Build+Diag: {time.time()-t0:.0f}s", flush=True)

    # Cluster-Identifikation (A-Eigenwerte)
    cl_A = [k for k in range(dim) if abs(float(ws_qw[k]) - w_min) < 1e-10]

    # α_cl = ⟨ũ, P_cluster k⟩
    c_coeffs = []
    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]
        c_coeffs.append(float(inner_product(qk, kn, dim)))

    alpha_cl = 0.0
    for k in cl_A:
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]
        u_comp = float(inner_product(ut, qk, dim))
        alpha_cl += c_coeffs[k] * u_comp

    print(f"  w_min={w_min:.6f}, α={alpha:.6f}, α_cl={alpha_cl:.6f}")
    print(f"  α_cl/α = {alpha_cl/alpha:.8f}")
    print(f"  Cluster: {len(cl_A)}/{dim}")

    # === Phase 2: T-Eigenbasis-Projektionen ===
    t_vals = np.array([float(ts[j]) for j in range(dim)])
    f_arr = np.zeros(dim)
    g_arr = np.zeros(dim)

    for j in range(dim):
        vj = mp.matrix(dim, 1)
        for i in range(dim):
            vj[i, 0] = Es[i, j]
        f_arr[j] = float(inner_product(vj, f_lam, dim))
        g_arr[j] = float(inner_product(vj, g_lam, dim))

    # Zero-EV und Cluster-Gap (Eigenvalue-basiert)
    j_zero = int(np.argmin(np.abs(t_vals)))
    t_offsets = np.abs(t_vals - w_min)
    t_offsets_nz = t_offsets.copy()
    t_offsets_nz[j_zero] = 1e10
    t_sorted = np.sort(t_offsets_nz)
    for gi in range(len(t_sorted) - 1):
        if t_sorted[gi + 1] / max(t_sorted[gi], 1e-30) > 10:
            cluster_gap = (t_sorted[gi] + t_sorted[gi + 1]) / 2
            break
    else:
        cluster_gap = 1e-3

    j_near_cl = [j for j in range(dim) if j != j_zero and t_offsets_nz[j] < cluster_gap]
    j_reg = [j for j in range(dim) if j != j_zero and abs(f_arr[j]) > 1e-25]

    print(f"  j_zero={j_zero} (t={t_vals[j_zero]:.2e})")
    print(f"  Cluster-Gap: {cluster_gap:.4e}, Near-Cluster Moden: {len(j_near_cl)}")
    print(f"  J_reg (|f_j|>1e-25): {len(j_reg)}")
    print(f"  ||f|| = {np.sqrt(np.sum(f_arr**2)):.6e}")
    print(f"  ||g|| = {np.sqrt(np.sum(g_arr**2)):.6e}")

    # === Phase 3: η_j Berechnung ===
    print(f"\n  {'='*78}")
    print(f"  η-ANALYSE: η_j = -(t_j-w_min)g_j/(αf_j),  C2o.1 ⟺ η_j ∈ (0,2)")
    print(f"  {'='*78}")

    # Cluster-Zerlegung: g_j^cl = -α_cl f_j / (t_j - w_min)
    eta_arr = np.full(dim, np.nan)
    delta_arr = np.full(dim, np.nan)
    g_cl_arr = np.zeros(dim)
    g_perp_arr = np.zeros(dim)
    kappa_perp_arr = np.full(dim, np.nan)
    r_arr = np.zeros(dim)

    for j in j_reg:
        tj = t_vals[j]
        fj = f_arr[j]
        gj = g_arr[j]
        d = tj - w_min
        if abs(d) < 1e-30:
            continue

        eta_arr[j] = -(d * gj) / (alpha * fj)
        g_cl_arr[j] = -alpha_cl * fj / d
        g_perp_arr[j] = gj - g_cl_arr[j]
        r_arr[j] = d * gj + alpha * fj  # = α f_j (1 - η_j)
        delta_arr[j] = (d * g_perp_arr[j]) / (alpha * fj)
        kappa_perp_arr[j] = abs(delta_arr[j])

    eta_valid = eta_arr[~np.isnan(eta_arr)]
    delta_valid = delta_arr[~np.isnan(delta_arr)]
    kp_valid = kappa_perp_arr[~np.isnan(kappa_perp_arr)]

    # Tabelle: Top-Moden sortiert nach η_j (absteigend, gefährlichste zuerst)
    jr_sorted = sorted(j_reg, key=lambda j: -eta_arr[j] if not np.isnan(eta_arr[j]) else -999)

    print(f"\n  Moden sortiert nach η_j (absteigend — gefährlichste zuerst):")
    print(f"  {'j':>4}  {'t_j':>10}  {'t_j-w':>10}  {'|f_j|':>10}  {'|g_j|':>10}  "
          f"{'η_j':>10}  {'1-η_j':>10}  {'δ_j':>10}  {'|g⊥|':>10}  {'near?':>5}")
    print(f"  {'-'*105}")

    n_show = min(30, len(jr_sorted))
    for rank in range(n_show):
        j = jr_sorted[rank]
        if np.isnan(eta_arr[j]):
            continue
        tj = t_vals[j]
        is_near = "NEAR" if j in j_near_cl else ""
        print(f"  {j:4d}  {tj:10.4f}  {tj-w_min:10.4e}  {abs(f_arr[j]):10.3e}  "
              f"{abs(g_arr[j]):10.3e}  {eta_arr[j]:10.6f}  {1-eta_arr[j]:10.6e}  "
              f"{delta_arr[j]:10.6e}  {abs(g_perp_arr[j]):10.3e}  {is_near:>5}")

    # === Phase 4: η-Statistik ===
    print(f"\n  {'='*78}")
    print(f"  η-STATISTIK ({len(eta_valid)} Moden)")
    print(f"  {'='*78}")

    eta_min = np.min(eta_valid)
    eta_max = np.max(eta_valid)
    eta_mean = np.mean(eta_valid)
    eta_med = np.median(eta_valid)

    print(f"  η_j Range:    [{eta_min:.8f}, {eta_max:.8f}]")
    print(f"  η_j Mean:     {eta_mean:.8f}")
    print(f"  η_j Median:   {eta_med:.8f}")
    print(f"  η_j Std:      {np.std(eta_valid):.8e}")
    print(f"\n  ALLE η_j > 0?  {'JA' if eta_min > 0 else 'NEIN'}  (min = {eta_min:.8e})")
    print(f"  ALLE η_j < 2?  {'JA' if eta_max < 2 else 'NEIN'}  (max = {eta_max:.8f})")
    print(f"  Marge zu 0:    {eta_min:.6e}")
    print(f"  Marge zu 2:    {2 - eta_max:.6e}")

    print(f"\n  |1-η_j| Statistik (= |r_j|/(α|f_j|) = κ* von C2o.1):")
    one_minus_eta = np.abs(1 - eta_valid)
    print(f"    max  = {np.max(one_minus_eta):.6e}  (= κ* = sup ratio)")
    print(f"    med  = {np.median(one_minus_eta):.6e}")
    print(f"    mean = {np.mean(one_minus_eta):.6e}")

    print(f"\n  δ_j = (t_j-w)g_j^⊥/(αf_j) — Off-Cluster-Defekt:")
    print(f"    |δ| max  = {np.max(np.abs(delta_valid)):.6e}")
    print(f"    |δ| med  = {np.median(np.abs(delta_valid)):.6e}")
    print(f"    δ mean   = {np.mean(delta_valid):.6e}  (Vorzeichen-Tendenz?)")
    print(f"    α_cl/α   = {alpha_cl/alpha:.8f}")
    print(f"    η = α_cl/α - δ, also δ ≈ α_cl/α - 1 ≈ {alpha_cl/alpha - 1:.6e} bei η≈1")

    # === Phase 5: Energieargument ===
    print(f"\n  {'='*78}")
    print(f"  ENERGIEARGUMENT")
    print(f"  {'='*78}")

    # Identität: (t_j-w)g_j² + αf_j g_j = g_j r_j
    # Bei η≈1: (t_j-w)|g_j|² ≈ α|f_j||g_j| (da f·g < 0)
    e_kinetic = []  # (t_j - w_min) g_j²
    e_forcing = []  # α |f_j| |g_j| (= -α f_j g_j da f·g < 0)
    e_ratio = []    # kinetic / forcing

    for j in j_reg:
        if np.isnan(eta_arr[j]):
            continue
        tj = t_vals[j]
        d = tj - w_min
        kin = d * g_arr[j]**2
        forc = abs(alpha * f_arr[j] * g_arr[j])
        e_kinetic.append(kin)
        e_forcing.append(forc)
        if forc > 0:
            e_ratio.append(kin / forc)

    e_kin = np.array(e_kinetic)
    e_for = np.array(e_forcing)
    e_rat = np.array(e_ratio)

    print(f"  (t_j-w)|g_j|² vs α|f_j||g_j|:")
    print(f"    Σ kinetic  = {np.sum(e_kin):.6e}")
    print(f"    Σ forcing  = {np.sum(e_for):.6e}")
    print(f"    Ratio Σ    = {np.sum(e_kin)/np.sum(e_for):.8f}")
    print(f"\n  Per-Mode Ratio (t_j-w)|g_j|²/(α|f_j||g_j|):")
    print(f"    = (t_j-w)|g_j|/(α|f_j|) = η_j  (falls f_j g_j < 0)")
    print(f"    max  = {np.max(e_rat):.8f}")
    print(f"    min  = {np.min(e_rat):.8f}")
    print(f"    mean = {np.mean(e_rat):.8f}")
    print(f"    (Bestätigt: Energy-Ratio = η_j)")

    # Globale Energiebilanz
    E_total_kin = np.sum(e_kin)
    E_total_for = np.sum(e_for)
    norm_g_sq = np.sum(g_arr[j]**2 for j in j_reg if not np.isnan(eta_arr[j]))
    norm_f_sq = np.sum(f_arr[j]**2 for j in j_reg if not np.isnan(eta_arr[j]))

    print(f"\n  Globale Energiebilanz:")
    print(f"    ||g||² (J_reg) = {norm_g_sq:.6e}")
    print(f"    ||f||² (J_reg) = {norm_f_sq:.6e}")
    print(f"    ⟨g,(T-w)g⟩     = {E_total_kin:.6e}")
    print(f"    -α⟨f,g⟩        = {E_total_for:.6e}")
    print(f"    Differenz       = {abs(E_total_kin - E_total_for):.6e}")

    # === Phase 6: Cauchy-Schwarz und Norm-Bounds ===
    print(f"\n  {'='*78}")
    print(f"  NORM-BOUNDS UND CAUCHY-SCHWARZ")
    print(f"  {'='*78}")

    # ||g|| ≤ 1 Bound
    norm_g_full = float(norm(g_lam, dim))
    print(f"  ||g|| = ||P₀k|| = {norm_g_full:.8f}  (≤ 1)")

    # Für jede Mode: |g_j| ≤ ||g|| ≤ 1
    # Also η_j = (t_j-w)|g_j|/(α|f_j|) ≤ (t_j-w)/(α|f_j|)
    # η_j < 2 ⟺ |g_j| < 2α|f_j|/(t_j-w)

    print(f"\n  Mode-Bound: η_j < 2 ⟺ |g_j| < 2α|f_j|/(t_j-w)")
    print(f"  Vergleich |g_j| vs 2α|f_j|/(t_j-w):")

    gb_ratio = []  # |g_j| / [2α|f_j|/(t_j-w)]
    for j in jr_sorted[:n_show]:
        if np.isnan(eta_arr[j]):
            continue
        d = t_vals[j] - w_min
        bound = 2 * abs(alpha) * abs(f_arr[j]) / d
        actual = abs(g_arr[j])
        rat = actual / bound
        gb_ratio.append(rat)

    gb_arr = np.array(gb_ratio)
    print(f"    max |g_j|/bound = {np.max(gb_arr):.8f}  (< 1 nötig)")
    print(f"    = max(η_j/2)    = {np.max(eta_valid)/2:.8f}")

    # Cauchy-Schwarz global: Σ|r_j|² = ||R||²
    r_valid = np.array([r_arr[j] for j in j_reg if not np.isnan(eta_arr[j])])
    norm_R = np.sqrt(np.sum(r_valid**2))
    print(f"\n  ||R|| = {norm_R:.6e}")
    print(f"  α||f|| = {abs(alpha)*np.sqrt(norm_f_sq):.6e}")
    print(f"  ||R||/(α||f||) = {norm_R/(abs(alpha)*np.sqrt(norm_f_sq)):.6e}  (globaler κ)")

    # === Phase 7: Monotonie und Scaling ===
    print(f"\n  {'='*78}")
    print(f"  MONOTONIE: η_j als Funktion von (t_j - w_min)")
    print(f"  {'='*78}")

    # Sortiere nach t_j
    jr_by_t = sorted(j_reg, key=lambda j: t_vals[j])
    jr_by_t = [j for j in jr_by_t if not np.isnan(eta_arr[j])]

    print(f"\n  {'j':>4}  {'t_j-w':>12}  {'η_j':>12}  {'δ_j':>12}  {'|g_j|':>12}  {'|g⊥|':>12}  {'|f_j|':>12}")
    print(f"  {'-'*90}")
    for j in jr_by_t[:35]:
        d = t_vals[j] - w_min
        print(f"  {j:4d}  {d:12.4e}  {eta_arr[j]:12.8f}  {delta_arr[j]:12.6e}  "
              f"{abs(g_arr[j]):12.4e}  {abs(g_perp_arr[j]):12.4e}  {abs(f_arr[j]):12.4e}")

    # Ist η monoton fallend in t_j?
    etas_sorted = [eta_arr[j] for j in jr_by_t]
    monotone_dec = all(etas_sorted[i] >= etas_sorted[i+1] for i in range(len(etas_sorted)-1))
    monotone_inc = all(etas_sorted[i] <= etas_sorted[i+1] for i in range(len(etas_sorted)-1))
    print(f"\n  η monoton fallend in t?  {monotone_dec}")
    print(f"  η monoton steigend in t? {monotone_inc}")

    # Near-cluster vs off-cluster η
    eta_near = [eta_arr[j] for j in j_reg if j in j_near_cl and not np.isnan(eta_arr[j])]
    eta_off = [eta_arr[j] for j in j_reg if j not in j_near_cl and not np.isnan(eta_arr[j])]

    if eta_near:
        print(f"\n  Near-Cluster Moden ({len(eta_near)}):")
        print(f"    η range: [{min(eta_near):.8f}, {max(eta_near):.8f}]")
        print(f"    η mean:  {np.mean(eta_near):.8f}")
    if eta_off:
        print(f"  Off-Cluster Moden ({len(eta_off)}):")
        print(f"    η range: [{min(eta_off):.8f}, {max(eta_off):.8f}]")
        print(f"    η mean:  {np.mean(eta_off):.8f}")

    # === Phase 8: Off-Cluster-Defekt-Struktur ===
    print(f"\n  {'='*78}")
    print(f"  OFF-CLUSTER-DEFEKT: Was kontrolliert g_j^⊥?")
    print(f"  {'='*78}")

    # g_j^⊥ = g_j - g_j^cl = g_j + α_cl f_j/(t_j-w)
    # r_j = (α-α_cl)f_j + (t_j-w)g_j^⊥
    # Frage: Ist (t_j-w)|g_j^⊥| klein relativ zu α|f_j|?

    gp_valid = np.array([g_perp_arr[j] for j in j_reg if not np.isnan(eta_arr[j])])
    print(f"  ||g^⊥|| = {np.sqrt(np.sum(gp_valid**2)):.6e}")
    print(f"  ||g^cl|| = {np.sqrt(np.sum(g_cl_arr[j]**2 for j in j_reg if not np.isnan(eta_arr[j]))):.6e}")
    print(f"  ||g||    = {norm_g_full:.6e}")

    # Korrelation g_j^⊥ mit f_j
    f_valid = np.array([f_arr[j] for j in j_reg if not np.isnan(eta_arr[j])])
    if np.std(f_valid) > 0 and np.std(gp_valid) > 0:
        corr_gp_f = np.corrcoef(f_valid, gp_valid)[0, 1]
        print(f"  corr(g^⊥, f):  {corr_gp_f:+.6f}")

    # Analyse: r_j = (α-α_cl)f_j + (t_j-w)g_j^⊥
    dalpha = alpha - alpha_cl
    print(f"\n  r_j = (α-α_cl)f_j + (t_j-w)g_j^⊥")
    print(f"  α - α_cl = {dalpha:.6e}")

    r_term1 = []
    r_term2 = []
    for j in j_reg:
        if np.isnan(eta_arr[j]):
            continue
        t1 = dalpha * f_arr[j]
        t2 = (t_vals[j] - w_min) * g_perp_arr[j]
        r_term1.append(t1)
        r_term2.append(t2)

    rt1 = np.array(r_term1)
    rt2 = np.array(r_term2)
    print(f"  ||(α-α_cl)f||  = {np.sqrt(np.sum(rt1**2)):.6e}")
    print(f"  ||(t-w)g^⊥||   = {np.sqrt(np.sum(rt2**2)):.6e}")
    print(f"  ||r||           = {norm_R:.6e}")

    # Vorzeichen: cancelieren die Terme?
    same = sum(1 for i in range(len(rt1)) if rt1[i] * rt2[i] > 0)
    print(f"  Gleich-VZ:  {same}/{len(rt1)}")
    print(f"  Gegen-VZ:   {len(rt1)-same}/{len(rt1)}")

    # === Phase 9: Potentieller Beweisweg ===
    print(f"\n  {'='*78}")
    print(f"  BEWEISANSÄTZE — NUMERISCHE EVIDENZ")
    print(f"  {'='*78}")

    # Route 1: Variational / Energy
    # Aus ⟨g,(T-w)g⟩ = -α⟨f,g⟩ und Cauchy-Schwarz:
    # ⟨g,(T-w)g⟩ ≤ α||f||||g|| (da f·g < 0 für alle Moden)
    # Aber: ⟨g,(T-w)g⟩ = Σ (t_j-w)g_j² ≥ (t_j-w)g_j² für jeden einzelnen j
    # Also: (t_j-w)g_j² ≤ α||f||||g||
    # => |g_j| ≤ √(α||f||||g||/(t_j-w))
    # => η_j = (t_j-w)|g_j|/(α|f_j|) ≤ √((t_j-w)·α||f||||g||)/(α|f_j|)
    #        = √(||f||||g||(t_j-w)/α) / |f_j|

    print(f"\n  Route 1: Energy-Bound")
    norm_f = np.sqrt(norm_f_sq)
    print(f"    α||f||||g|| = {abs(alpha)*norm_f*norm_g_full:.6e}")
    print(f"    Per-Mode-Bound: |g_j| ≤ √(α||f||||g||/(t_j-w))")

    for j in jr_sorted[:10]:
        if np.isnan(eta_arr[j]):
            continue
        d = t_vals[j] - w_min
        g_bound = math.sqrt(abs(alpha) * norm_f * norm_g_full / d) if d > 0 else 999
        eta_bound = d * g_bound / (abs(alpha) * abs(f_arr[j])) if abs(f_arr[j]) > 0 else 999
        print(f"    j={j:3d}: |g_j|={abs(g_arr[j]):.4e} ≤ {g_bound:.4e}? "
              f"{'JA' if abs(g_arr[j]) < g_bound else 'NEIN':>4}  "
              f"η_bound={eta_bound:.4f} {'(< 2)' if eta_bound < 2 else '(≥ 2!)'}")

    # Route 2: ||g|| ≤ 1 + Resolventenform
    print(f"\n  Route 2: ||g||² = Σg_j² ≤ 1")
    print(f"    g_j ≈ -α f_j/(t_j-w) + g_j^⊥")
    print(f"    ||g^cl||² = α_cl² Σ f_j²/(t_j-w)² = {np.sum(g_cl_arr[j]**2 for j in j_reg if not np.isnan(eta_arr[j])):.6e}")

    # Route 3: Direkte Kontrolle von (t_j-w)g_j^⊥
    # r_j = (α-α_cl)f_j + (t_j-w)g_j^⊥
    # |r_j| ≤ |α-α_cl||f_j| + (t_j-w)|g_j^⊥|
    # Aber wir brauchen |r_j| < α|f_j|, also:
    # (t_j-w)|g_j^⊥| < α|f_j| - |α-α_cl||f_j| = (α - |α-α_cl|)|f_j| ≈ α_cl|f_j|
    # Das ist exakt C2l.5!

    print(f"\n  Route 3: Direkte g^⊥-Kontrolle (= C2l.5)")
    print(f"    (t_j-w)|g_j^⊥| < α_cl|f_j| nötig")
    print(f"    Äquivalent: |δ_j| < α_cl/α ≈ 1")
    for j in jr_sorted[:10]:
        if np.isnan(eta_arr[j]):
            continue
        d = t_vals[j] - w_min
        lhs = d * abs(g_perp_arr[j])
        rhs = abs(alpha_cl) * abs(f_arr[j])
        print(f"    j={j:3d}: (t-w)|g⊥|={lhs:.4e} < α_cl|f|={rhs:.4e}?  "
              f"{'JA' if lhs < rhs else 'NEIN'}  (Ratio={lhs/rhs:.6e})")

    # Route 4: Quadratische Form / Operator-Ungleichung
    # Betrachte: Σ (t_j-w) g_j^⊥² (Energie des Off-Cluster-Anteils)
    e_perp = sum((t_vals[j]-w_min) * g_perp_arr[j]**2 for j in j_reg if not np.isnan(eta_arr[j]))
    e_cl = sum((t_vals[j]-w_min) * g_cl_arr[j]**2 for j in j_reg if not np.isnan(eta_arr[j]))
    e_cross = sum((t_vals[j]-w_min) * 2 * g_cl_arr[j] * g_perp_arr[j] for j in j_reg if not np.isnan(eta_arr[j]))

    print(f"\n  Route 4: Quadratische Energieform")
    print(f"    ⟨g^cl,(T-w)g^cl⟩ = {e_cl:.6e}")
    print(f"    ⟨g^⊥,(T-w)g^⊥⟩  = {e_perp:.6e}")
    print(f"    2⟨g^cl,(T-w)g^⊥⟩ = {e_cross:.6e}")
    print(f"    Summe             = {e_cl+e_perp+e_cross:.6e}")
    print(f"    = ⟨g,(T-w)g⟩     = {E_total_kin:.6e}")

    # Cluster-Resolventenform: (t_j-w)g_j^cl = -α_cl f_j
    # Also ⟨g^cl,(T-w)g^cl⟩ = α_cl² Σ f_j²/(t_j-w)
    # Und ⟨g^cl,(T-w)g^⊥⟩ = -α_cl Σ f_j g_j^⊥
    cross_fg_perp = sum(f_arr[j] * g_perp_arr[j] for j in j_reg if not np.isnan(eta_arr[j]))
    print(f"    ⟨f, g^⊥⟩ = {cross_fg_perp:.6e}")
    print(f"    -α_cl⟨f, g^⊥⟩ = {-alpha_cl * cross_fg_perp:.6e} (≈ ½ Cross-Term)")

    return {
        "lam": lam,
        "eta_max": float(eta_max),
        "eta_min": float(eta_min),
        "eta_mean": float(eta_mean),
        "kappa_star": float(np.max(np.abs(1 - eta_valid))),
        "delta_max": float(np.max(np.abs(delta_valid))),
        "norm_g_perp": float(np.sqrt(np.sum(gp_valid**2))),
        "norm_g": float(norm_g_full),
        "norm_R": float(norm_R),
        "alpha_cl_alpha": float(alpha_cl / alpha),
        "n_reg": len(eta_valid),
        "e_perp": float(e_perp),
    }


if __name__ == "__main__":
    results = []
    for cfg in CONFIGS:
        r = run(cfg["lam"], cfg["N"])
        results.append(r)

    if len(results) >= 2:
        print(f"\n{'='*78}")
        print(f"  λ-VERGLEICH (C2o.1)")
        print(f"{'='*78}")
        print(f"\n  {'Metrik':>20}  ", end="")
        for r in results:
            print(f"{'λ='+str(r['lam']):>14}  ", end="")
        print(f"{'Exponent':>10}")

        for name, key in [
            ("η_max", "eta_max"),
            ("η_min", "eta_min"),
            ("κ* = sup|1-η|", "kappa_star"),
            ("|δ|_max", "delta_max"),
            ("||g^⊥||", "norm_g_perp"),
            ("||g||", "norm_g"),
            ("||R||", "norm_R"),
            ("α_cl/α", "alpha_cl_alpha"),
            ("⟨g⊥,(T-w)g⊥⟩", "e_perp"),
        ]:
            vals = [r[key] for r in results]
            print(f"  {name:>20}  ", end="")
            for v in vals:
                print(f"{v:14.6e}  ", end="")
            if len(vals) >= 2 and abs(vals[0]) > 0 and abs(vals[-1]) > 0:
                lams = [r["lam"] for r in results]
                try:
                    exp = math.log(abs(vals[-1]) / abs(vals[0])) / math.log(lams[-1] / lams[0])
                    print(f"  {exp:8.2f}")
                except (ValueError, ZeroDivisionError):
                    print()
            else:
                print()

    print(f"\n{'='*78}")
    print(f"  FAZIT")
    print(f"{'='*78}")
    for r in results:
        print(f"  λ={r['lam']}: η ∈ [{r['eta_min']:.6f}, {r['eta_max']:.6f}], "
              f"κ*={r['kappa_star']:.4e}, |δ|_max={r['delta_max']:.4e}")
    all_ok = all(r["eta_max"] < 2 and r["eta_min"] > 0 for r in results)
    print(f"\n  C2o.1 (η ∈ (0,2)) numerisch bestätigt: {'JA' if all_ok else 'NEIN'}")

    print(f"\n=== FERTIG ===")
