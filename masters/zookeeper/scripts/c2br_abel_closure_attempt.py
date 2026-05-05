"""
c2br_abel_closure_attempt.py

Versuch, die Abel-Schließung (C2aj: |δ_m| < 1) nicht-zirkulär zu beweisen.

Strategie: Zerlege δ_j in strukturelle Komponenten und identifiziere,
welcher Anteil OHNE MS2-Annahme kontrollierbar ist.

Zerlegung (aus C2w/C2y):
  δ_j = -(d_j/α) · Σ_{off} c_a α_a / (t_j - w_a)
      = -(d_j/α) · Σ_{off} [(μ/α)α_a² + α_a h_a] / [(w_a-w_min)(t_j-w_a)]
      = -(μ/α)·(d_j)·Σ α_a²/[(w_a-w_min)(t_j-w_a)]  [Term I: kanonisch]
        -(d_j/α)·Σ α_a h_a/[(w_a-w_min)(t_j-w_a)]     [Term II: h-Korrektur]

Term I nutzt C2v (Rang-1-Nullidentität):
  Σ_{off} α_a²/(w_a - t_j) = A_cl/d_j  (exakt)

Schlüsselfrage: Kann man |δ_j| < 1 beweisen OHNE μ/α oder ||k_perp|| klein vorauszusetzen?

Neue Idee: Die GLÄTTE von k_λ erzwingt Abfall der Koeffizienten c_a mit wachsendem
Eigenwertindex — das ist eine Eigenschaft von k (Poisson-Struktur), nicht von MS2.
Wenn der Abfall schnell genug ist, folgt die Abel-Schranke NICHT-ZIRKULÄR.
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import mpmath as mp
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, chi_trivial,
)
from c2_approximation_test import (
    h_educated_guess, k_lambda_value, norm,
)
from c2_poisson_decomposition import (
    h_hat_analytical, project_to_fourier,
)

DPS = int(os.environ.get("DPS", 25))


def to_np(M, dim):
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = float(M[i, j])
    return A


def compute_abel_structure(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)
    L = float(L_mp)

    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    A_np = to_np(Aq, dim)
    H_np = to_np(Ah, dim)
    W_np = A_np - H_np

    col0 = W_np[:, 0]
    cn = np.linalg.norm(col0)
    ut = col0 / cn

    Ct = float(ut @ W_np @ ut)

    ws_A, vs_A = np.linalg.eigh(A_np)
    w0 = ws_A[0]
    g0 = vs_A[:, 0]
    if g0[0] < 0:
        g0 = -g0

    hs_T, vs_T = np.linalg.eigh(H_np)

    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)
    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
    alpha = float(ut @ kn)

    # P_perp
    P_perp = np.eye(dim) - np.outer(ut, ut)

    # T-Eigensystem (= PHP restricted to P_perp-space)
    T_np = P_perp @ H_np @ P_perp
    ts, vs_T_full = np.linalg.eigh(T_np)

    # A-Eigenwerte und -Vektoren
    # Cluster: Eigenwerte nahe w0
    cluster_tol = 1e-6
    cluster_mask = np.abs(ws_A - w0) < cluster_tol
    n_cl = np.sum(cluster_mask)
    off_mask = ~cluster_mask
    off_indices = np.where(off_mask)[0]
    cl_indices = np.where(cluster_mask)[0]

    # T-Eigenwerte (off-cluster: die mit Eigenwert != 0 aus T_np)
    # T hat einen Nulleigenwert in ut-Richtung
    t_sort = np.sort(ts)
    # Identifiziere off-zero T-Eigenwerte
    t_nonzero_mask = np.abs(ts) > 1e-10
    t_off_indices = np.where(t_nonzero_mask)[0]
    t_off_vals = ts[t_off_indices]
    t_off_vecs = vs_T_full[:, t_off_indices]

    # Sortiere T-Eigenwerte
    t_sort_idx = np.argsort(t_off_vals)
    t_off_vals = t_off_vals[t_sort_idx]
    t_off_vecs = t_off_vecs[:, t_sort_idx]

    # Off-cluster A-Eigenwerte
    w_off = ws_A[off_indices]
    v_off = vs_A[:, off_indices]
    off_sort = np.argsort(w_off)
    w_off = w_off[off_sort]
    v_off = v_off[:, off_sort]

    n_off = len(w_off)
    n_T = len(t_off_vals)

    # alpha_a = <q_a, ut> fuer off-cluster
    alpha_a = np.array([float(v_off[:, a] @ ut) for a in range(n_off)])

    # c_a = <q_a, k> fuer off-cluster
    c_a = np.array([float(v_off[:, a] @ kn) for a in range(n_off)])

    # f_j = <e_j, f> wo f = H*ut - <ut,H*ut>*ut
    uHu = float(ut @ H_np @ ut)
    f_vec = H_np @ ut - uHu * ut
    f_j = np.array([float(t_off_vecs[:, j] @ f_vec) for j in range(n_T)])

    # A_cl (cluster mass in canonical weights)
    A_cl = sum(float(vs_A[:, i] @ ut)**2 for i in cl_indices)

    # mu/alpha
    mu_alpha = float(ut @ (A_np - w0 * np.eye(dim)) @ kn) / alpha

    # h = P_perp(A-w0)k
    h_vec = P_perp @ (A_np - w0 * np.eye(dim)) @ kn
    h_a = np.array([float(v_off[:, a] @ h_vec) for a in range(n_off)])

    # === DELTA-BERECHNUNG ===
    # delta_j fuer jedes T-Eigenwert t_j
    deltas = np.zeros(n_T)
    delta_term1 = np.zeros(n_T)  # kanonischer Term (proportional zu mu/alpha)
    delta_term2 = np.zeros(n_T)  # h-Korrektur

    for j in range(n_T):
        t_j = t_off_vals[j]
        d_j = t_j - w0  # distance to cluster

        if abs(d_j) < 1e-15:
            continue

        # Full delta
        cauchy_sum = 0.0
        for a in range(n_off):
            denom = t_j - w_off[a]
            if abs(denom) > 1e-15:
                cauchy_sum += c_a[a] * alpha_a[a] / denom

        deltas[j] = 1.0 + (d_j / alpha) * cauchy_sum

        # Decompose: c_a*alpha_a = (mu*alpha_a^2 + alpha_a*h_a) / (w_a - w0)
        # So the Cauchy sum = mu/alpha * Σ alpha_a^2/[(w_a-w0)(t_j-w_a)] + ...
        sum1 = 0.0
        sum2 = 0.0
        for a in range(n_off):
            denom = t_j - w_off[a]
            w_gap = w_off[a] - w0
            if abs(denom) > 1e-15 and abs(w_gap) > 1e-15:
                sum1 += alpha_a[a]**2 / (w_gap * denom)
                sum2 += alpha_a[a] * h_a[a] / (w_gap * denom)

        delta_term1[j] = (d_j / alpha) * (mu_alpha * alpha) * sum1  # = mu * d_j * sum1 / alpha...
        # Actually: term1 = (d_j/alpha) * (mu_alpha*alpha/1) * sum1... let me redo this
        # c_a = (mu*alpha_a + h_a) / (w_a - w0)
        # c_a * alpha_a = (mu*alpha_a^2 + alpha_a*h_a) / (w_a - w0)
        # Cauchy = Σ c_a*alpha_a/(t_j-w_a) = Σ (mu*alpha_a^2 + alpha_a*h_a) / [(w_a-w0)(t_j-w_a)]
        # delta = 1 + (d_j/alpha) * Cauchy
        # term1 part: (d_j/alpha) * mu * Σ alpha_a^2/[(w_a-w0)(t_j-w_a)]
        # But mu = alpha * mu_alpha, so term1 = d_j * mu_alpha * Σ alpha_a^2/[(w_a-w0)(t_j-w_a)]
        delta_term1[j] = d_j * mu_alpha * sum1
        delta_term2[j] = (d_j / alpha) * sum2

    # === ABEL PARTIAL SUMS ===
    # Off-cluster Abel partial sums of c_a * alpha_a (gewichtet, geordnet nach w_a)
    weights = c_a * alpha_a
    B_partial = np.cumsum(weights)  # B_m = Σ_{a<=m} c_a α_a
    C_partial = np.sum(weights) - B_partial  # C_{m+1} = Σ_{a>m} c_a α_a

    # Canonical partial sums alpha_a^2
    can_weights = alpha_a**2
    B_can = np.cumsum(can_weights)
    C_can = np.sum(can_weights) - B_can

    # === ABFALL VON c_a (Glätte-Hypothese) ===
    # Erwartung: |c_a| fällt mit wachsendem Eigenwert w_a
    # Weil k_λ glatt ist, sollten hohe Eigenmoden wenig Overlap haben

    # === ABEL BOUND BERECHNUNG ===
    # Für jeden T-Eigenwert t_j: finde den "gefährlichen" Abel-Kanal
    abel_bounds = np.zeros(n_T)
    abel_left = np.zeros(n_T)
    abel_right = np.zeros(n_T)

    for j in range(n_T):
        t_j = t_off_vals[j]
        d_j = t_j - w0
        if abs(d_j) < 1e-15:
            continue

        # Finde Position von t_j relativ zu off-cluster A-Eigenwerten
        # t_j liegt zwischen w_off[m] und w_off[m+1] (Interlacing)
        pos = np.searchsorted(w_off, t_j)

        # Linker Kanal: B_{pos-1} / (t_j - w_{pos-1}) * d_j
        if pos > 0:
            gap_L = t_j - w_off[pos - 1]
            if gap_L > 1e-15:
                abel_left[j] = abs(B_partial[pos - 1]) * d_j / (alpha * gap_L)

        # Rechter Kanal: C_{pos} / (w_{pos} - t_j) * d_j
        if pos < n_off:
            gap_R = w_off[pos] - t_j
            if gap_R > 1e-15:
                abel_right[j] = abs(C_partial[pos - 1] if pos > 0 else np.sum(weights)) * d_j / (alpha * gap_R)
                # Korrektur: C_{pos} = Σ_{a>=pos} c_a α_a
                C_from_pos = np.sum(weights[pos:])
                abel_right[j] = abs(C_from_pos) * d_j / (alpha * gap_R)

        abel_bounds[j] = max(abel_left[j], abel_right[j])

    return {
        'lam': lam, 'N': N, 'dim': dim,
        'n_off': n_off, 'n_T': n_T, 'n_cl': n_cl,
        'w0': w0, 'alpha': alpha, 'mu_alpha': mu_alpha, 'A_cl': A_cl,
        'w_off': w_off, 't_off': t_off_vals,
        'c_a': c_a, 'alpha_a': alpha_a, 'h_a': h_a,
        'deltas': deltas, 'delta_term1': delta_term1, 'delta_term2': delta_term2,
        'B_partial': B_partial, 'C_partial': C_partial,
        'B_can': B_can, 'C_can': C_can,
        'abel_bounds': abel_bounds, 'abel_left': abel_left, 'abel_right': abel_right,
        'f_j': f_j,
        'h_norm': np.linalg.norm(h_vec),
        'k_perp_sq': np.sum(c_a**2),
    }


def main():
    print("C2br: Abel-Schließungsversuch")
    print(f"DPS={DPS}")

    for lam in [3.0, 5.0]:
        N = 20
        print(f"\n{'='*90}")
        print(f"  lam={lam}, N={N}")
        print(f"{'='*90}")

        t0 = time.time()
        d = compute_abel_structure(lam, N)
        print(f"  Berechnung: {time.time()-t0:.0f}s", flush=True)

        print(f"\n  Grunddaten:")
        print(f"  n_cl={d['n_cl']}, n_off={d['n_off']}, n_T={d['n_T']}")
        print(f"  w0={d['w0']:.10e}, alpha={d['alpha']:.6f}")
        print(f"  mu/alpha={d['mu_alpha']:.6e}, A_cl={d['A_cl']:.6f}")
        print(f"  ||h||={d['h_norm']:.6e}, ||k_perp||^2={d['k_perp_sq']:.6e}")

        # Abfall von c_a
        print(f"\n  === ABFALL DER OFF-CLUSTER-KOEFFIZIENTEN ===")
        print(f"  {'a':>4} {'w_a-w0':>14} {'|c_a|':>12} {'|alpha_a|':>12} {'|c_a*alpha_a|':>14} {'|c_a|/|alpha_a|':>16}")
        n_show = min(10, d['n_off'])
        for a in range(n_show):
            ca = d['c_a'][a]
            aa = d['alpha_a'][a]
            print(f"  {a:4d} {d['w_off'][a]-d['w0']:+14.6e} {abs(ca):12.4e} {abs(aa):12.4e} "
                  f"{abs(ca*aa):14.6e} {abs(ca)/abs(aa) if abs(aa)>1e-30 else 0:16.6e}")
        if d['n_off'] > 10:
            print(f"  ...")
            for a in range(d['n_off']-3, d['n_off']):
                ca = d['c_a'][a]
                aa = d['alpha_a'][a]
                print(f"  {a:4d} {d['w_off'][a]-d['w0']:+14.6e} {abs(ca):12.4e} {abs(aa):12.4e} "
                      f"{abs(ca*aa):14.6e} {abs(ca)/abs(aa) if abs(aa)>1e-30 else 0:16.6e}")

        # Delta-Zerlegung
        print(f"\n  === DELTA-ZERLEGUNG (erste 10 T-Moden) ===")
        print(f"  {'j':>4} {'delta_j':>14} {'Term1(mu/a)':>14} {'Term2(h)':>14} {'|d|<1?':>8}")
        n_show_t = min(10, d['n_T'])
        for j in range(n_show_t):
            ok = "OK" if abs(d['deltas'][j]) < 1 else "FAIL"
            print(f"  {j:4d} {d['deltas'][j]:+14.6e} {d['delta_term1'][j]:+14.6e} "
                  f"{d['delta_term2'][j]:+14.6e} {ok:>8}")

        # Wieviel kommt von Term1 vs Term2?
        valid = np.abs(d['deltas']) > 1e-15
        if np.any(valid):
            ratio_t1 = np.abs(d['delta_term1'][valid]) / (np.abs(d['deltas'][valid]) + 1e-30)
            ratio_t2 = np.abs(d['delta_term2'][valid]) / (np.abs(d['deltas'][valid]) + 1e-30)
            print(f"\n  Term1-Anteil (median): {np.median(ratio_t1):.4f}")
            print(f"  Term2-Anteil (median): {np.median(ratio_t2):.4f}")
            print(f"  Term1 dominiert: {np.sum(ratio_t1 > ratio_t2)}/{np.sum(valid)}")

        # Abel-Bounds
        print(f"\n  === ABEL-BOUNDS ===")
        valid_abel = d['abel_bounds'] > 1e-15
        if np.any(valid_abel):
            max_abel = np.max(d['abel_bounds'][valid_abel])
            n_pass = np.sum(d['abel_bounds'][valid_abel] < 1)
            n_total = np.sum(valid_abel)
            print(f"  max(Abel bound) = {max_abel:.6e}")
            print(f"  Pass (<1): {n_pass}/{n_total}")

        # Kanonische vs. tatsächliche Partialsummen
        print(f"\n  === PARTIALSUMMEN-VERGLEICH ===")
        print(f"  {'m':>4} {'B_m (actual)':>14} {'B_m (can)':>14} {'Ratio':>12} {'C_m (actual)':>14} {'C_m (can)':>14}")
        for m in range(min(8, d['n_off'])):
            B_act = d['B_partial'][m]
            B_can = d['B_can'][m]
            C_act = np.sum(d['c_a'][m+1:] * d['alpha_a'][m+1:]) if m < d['n_off']-1 else 0
            C_can_m = d['C_can'][m]
            ratio = B_act / B_can if abs(B_can) > 1e-30 else 0
            print(f"  {m:4d} {B_act:+14.6e} {B_can:+14.6e} {ratio:+12.4e} {C_act:+14.6e} {C_can_m:+14.6e}")

        # Schlüsselfrage: Ist |c_a/alpha_a| beschränkt ohne MS2?
        print(f"\n  === SCHLÜSSELFRAGE: |c_a/alpha_a| BESCHRÄNKT? ===")
        ratios_ca = np.abs(d['c_a']) / (np.abs(d['alpha_a']) + 1e-30)
        valid_r = np.abs(d['alpha_a']) > 1e-10
        if np.any(valid_r):
            print(f"  max |c_a/alpha_a| = {np.max(ratios_ca[valid_r]):.6e}")
            print(f"  median            = {np.median(ratios_ca[valid_r]):.6e}")
            print(f"  Erwartung (mu/alpha * alpha_a / (w_a-w0)... ): {abs(d['mu_alpha'])*np.max(np.abs(d['alpha_a'])/(np.abs(d['w_off']-d['w0'])+1e-30)):.6e}")

        # Kontraktionstest: Wenn ||k_perp||=1 (trivial), was gibt |delta|?
        # Trivial-Bound: |delta_j| <= (d_j/alpha) * ||k_perp|| * ||alpha_off||_2 / min_gap
        # Das ist die C2y.3-Version
        d_vals = d['t_off'][:d['n_T']] - d['w0']
        alpha_off_norm = np.sqrt(np.sum(d['alpha_a']**2))
        print(f"\n  === KONTRAKTIONSTEST ===")
        print(f"  ||alpha_off||_2 = {alpha_off_norm:.6e}")
        print(f"  ||k_perp|| (actual) = {np.sqrt(d['k_perp_sq']):.6e}")
        # Trivialer Bound mit ||k_perp||=1:
        # |delta_j| <= d_j * sqrt(sum alpha_a^2/(t_j-w_a)^2) / alpha [C2y.3 mit ||R||=||k_perp||*max_gap]
        # Einfacher: |c_a| <= ||k_perp|| = 1 (trivial), also |c_a*alpha_a| <= |alpha_a|
        # => Abel bound mit |c_a*alpha_a| <= |alpha_a|
        trivial_B = np.cumsum(np.abs(d['alpha_a']))
        max_trivial_abel = 0.0
        for j in range(d['n_T']):
            t_j = d['t_off'][j]
            d_j = t_j - d['w0']
            if abs(d_j) < 1e-15:
                continue
            pos = np.searchsorted(d['w_off'], t_j)
            if pos > 0:
                gap_L = t_j - d['w_off'][pos-1]
                if gap_L > 1e-15:
                    ab = trivial_B[pos-1] * d_j / (d['alpha'] * gap_L)
                    max_trivial_abel = max(max_trivial_abel, ab)
            if pos < d['n_off']:
                gap_R = d['w_off'][pos] - t_j
                C_triv = np.sum(np.abs(d['alpha_a'][pos:]))
                if gap_R > 1e-15:
                    ab = C_triv * d_j / (d['alpha'] * gap_R)
                    max_trivial_abel = max(max_trivial_abel, ab)

        print(f"  max(Abel bound mit |c_a|<=1): {max_trivial_abel:.6e}")
        print(f"  {'KONTRAKTION MOEGLICH' if max_trivial_abel < 1 else 'KONTRAKTION SCHEITERT (bound > 1)'}")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
