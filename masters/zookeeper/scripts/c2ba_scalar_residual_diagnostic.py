"""
c2ba_scalar_residual_diagnostic.py — Diagnostik für C2ba.2

Zerlegung von μ = ⟨ũ, (A-w_min)k⟩ in Poisson-Bestandteile:

  μ = α(w̄ - w_min) + ⟨f_λ, k_⊥⟩

mit w̄ = C̃ + ⟨ũ,Hũ⟩ (Rayleigh-Quotient von ũ für A).

Frage: Welcher Term dominiert? Kann man μ/α klein machen
ohne ||k_⊥|| (= MS2) zu benutzen?

Aus C2M: (H_L - w_min)K_L = ((H_∞ - w_min)K)|_{[0,L]} + B_L
Also: ⟨ũ, (H-w_min)k⟩ = ⟨ũ, Bulk⟩ + ⟨ũ, B_L⟩
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import mpmath as mp
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


def scalar_residual_diagnostic(lam, N):
    mp.mp.dps = DPS
    dim = N + 1
    print(f"\n{'='*100}", flush=True)
    print(f"C2ba.2 Scalar Residual Diagnostic: λ={lam}, N={N}", flush=True)
    print(f"{'='*100}", flush=True)

    t0 = time.time()
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq, N, parity="even")
    Ah, _ = project_to_parity(Mh, N, parity="even")

    # W02 = A - H, also C̃|ũ⟩⟨ũ| = Aq - Ah
    W02 = mp.matrix(dim, dim)
    H = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            W02[i, j] = Aq[i, j] - Ah[i, j]
            H[i, j] = Ah[i, j]

    # ũ = normalisierte erste Spalte von W02
    col0 = mp.matrix(dim, 1)
    for i in range(dim):
        col0[i, 0] = W02[i, 0]
    cn = norm(col0, dim)
    ut = mp.matrix(dim, 1)
    for i in range(dim):
        ut[i, 0] = col0[i, 0] / cn

    # C̃ = ||W02 col0||² / ||col0||² × ... actually C̃ = cn² / (something)
    # Simpler: C̃ = ⟨ũ, W02 ũ⟩ / ⟨ũ, ũ⟩ = ⟨ũ, (A-H)ũ⟩
    # W02 = C̃ |ũ⟩⟨ũ|, so ⟨ũ, W02 ũ⟩ = C̃
    Ctilde = float(sum(ut[i, 0] * sum(W02[i, j] * ut[j, 0] for j in range(dim))
                       for i in range(dim)))

    # uHu = ⟨ũ, Hũ⟩
    Hut = mp.matrix(dim, 1)
    for i in range(dim):
        Hut[i, 0] = sum(H[i, j] * ut[j, 0] for j in range(dim))
    uHu = float(inner_product(ut, Hut, dim))

    # w̄ = C̃ + uHu = Rayleigh-Quotient von ũ für A
    w_bar = Ctilde + uHu

    # k_λ im Fourier-Raum
    def kfull(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf

    alpha = float(inner_product(ut, kn, dim))

    # k_⊥ = k - α ũ
    k_perp = mp.matrix(dim, 1)
    for i in range(dim):
        k_perp[i, 0] = kn[i, 0] - mp.mpf(alpha) * ut[i, 0]
    k_perp_norm = float(norm(k_perp, dim))

    # A-Diagonalisierung
    ws_A, Qs = diagonalize_mp(Aq, verbose=False)
    w_min = float(ws_A[0])

    print(f"Build: {time.time()-t0:.1f}s", flush=True)

    # --- Kerngrößen ---
    print(f"\n--- Kerngrößen ---", flush=True)
    print(f"  C̃              = {Ctilde:.6f}", flush=True)
    print(f"  ⟨ũ,Hũ⟩         = {uHu:.6f}", flush=True)
    print(f"  w̄ = C̃ + uHu   = {w_bar:.6f}", flush=True)
    print(f"  w_min           = {w_min:.6f}", flush=True)
    print(f"  w̄ - w_min      = {w_bar - w_min:.6e}  (Rayleigh-Gap)", flush=True)
    print(f"  α               = {alpha:.6f}", flush=True)
    print(f"  ||k_⊥||         = {k_perp_norm:.6e}", flush=True)

    # --- μ-Zerlegung ---
    # μ = ⟨ũ, (A - w_min)k⟩
    # Direkt: über A-Eigenbasis
    c_arr = np.zeros(dim)
    alpha_a = np.zeros(dim)
    w_arr = np.array([float(ws_A[a]) for a in range(dim)])
    for a in range(dim):
        qa = mp.matrix(dim, 1)
        for i in range(dim):
            qa[i, 0] = Qs[i, a]
        c_arr[a] = float(inner_product(qa, kn, dim))
        alpha_a[a] = float(inner_product(ut, qa, dim))

    cl_A = [a for a in range(dim) if abs(w_arr[a] - w_min) < 1e-10]
    noncl_A = [a for a in range(dim) if a not in cl_A]

    mu_spectral = sum((w_arr[a] - w_min) * alpha_a[a] * c_arr[a] for a in range(dim))

    # Zerlegung: μ = α(w̄ - w_min) + ⟨f_λ, k_⊥⟩
    term1 = alpha * (w_bar - w_min)

    # f_λ = Hũ - uHu·ũ = P₀Hũ
    f_lam = mp.matrix(dim, 1)
    for i in range(dim):
        f_lam[i, 0] = Hut[i, 0] - mp.mpf(uHu) * ut[i, 0]

    f_lam_norm = float(norm(f_lam, dim))
    term2 = float(inner_product(f_lam, k_perp, dim))

    mu_decomposed = term1 + term2

    print(f"\n--- μ-Zerlegung ---", flush=True)
    print(f"  μ (spektral)          = {mu_spectral:.6e}", flush=True)
    print(f"  μ (zerlegt)           = {mu_decomposed:.6e}", flush=True)
    print(f"  Konsistenz-Check      = {abs(mu_spectral - mu_decomposed):.2e}", flush=True)
    print(f"", flush=True)
    print(f"  Term 1: α(w̄-w_min)   = {term1:.6e}  (Rayleigh-Gap × α)", flush=True)
    print(f"  Term 2: ⟨f_λ, k_⊥⟩   = {term2:.6e}  (Poisson-Alignment)", flush=True)
    print(f"  Anteil Term 1         = {abs(term1)/abs(mu_spectral)*100:.1f}%", flush=True)
    print(f"  Anteil Term 2         = {abs(term2)/abs(mu_spectral)*100:.1f}%", flush=True)
    print(f"  CANCELLATION?         = {'JA' if np.sign(term1) != np.sign(term2) else 'NEIN'} "
          f"(Vorzeichen: {'+' if term1 > 0 else '-'} vs {'+' if term2 > 0 else '-'})", flush=True)
    print(f"  |T1|+|T2|           = {abs(term1)+abs(term2):.6e}", flush=True)
    print(f"  |μ|/(|T1|+|T2|)    = {abs(mu_spectral)/(abs(term1)+abs(term2)):.4f}  "
          f"(1=keine Cancellation, 0=perfekte)", flush=True)

    # --- μ/α ---
    mu_over_alpha = mu_spectral / alpha
    print(f"\n--- μ/α ---", flush=True)
    print(f"  μ/α                   = {mu_over_alpha:.6e}", flush=True)
    print(f"  = (w̄-w_min) + ⟨f_λ,k_⊥⟩/α", flush=True)
    print(f"  Term 1/α: w̄-w_min    = {w_bar - w_min:.6e}", flush=True)
    print(f"  Term 2/α: ⟨f_λ,k_⊥⟩/α = {term2/alpha:.6e}", flush=True)

    # --- Triviale Schranken ---
    print(f"\n--- Schranken ---", flush=True)
    print(f"  Trivial: |μ/α| ≤ (w̄-w_min) + ||f_λ||·||k_⊥||/α", flush=True)
    print(f"           = {w_bar-w_min:.6e} + {f_lam_norm*k_perp_norm/alpha:.6e}", flush=True)
    print(f"           = {(w_bar-w_min) + f_lam_norm*k_perp_norm/alpha:.6e}", flush=True)
    print(f"  ||f_λ||              = {f_lam_norm:.6e}", flush=True)
    print(f"  Ohne MS2-Term:  (w̄-w_min) allein = {w_bar-w_min:.6e}", flush=True)
    print(f"  Frage: Kann (w̄-w_min) allein μ/α erklären?", flush=True)

    # --- Feinere spektrale Analyse: ⟨f_λ, k_⊥⟩ aufgelöst ---
    print(f"\n--- ⟨f_λ, k_⊥⟩ nach A-Eigenmoden ---", flush=True)
    # ⟨f_λ, k_⊥⟩ = Σ_{a∉Cl} f_a · c_a  wobei f_a = (w_a - w̄)α_a
    f_a_vals = np.array([(w_arr[a] - w_bar) * alpha_a[a] for a in range(dim)])
    print(f"  {'a':>3}  {'w_a-w̄':>12}  {'α_a':>10}  {'c_a':>10}  {'f_a·c_a':>12}  {'kum':>12}", flush=True)

    kum = 0
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    for idx, a in enumerate(noncl_sorted[:8]):
        fa = f_a_vals[a]
        contrib = fa * c_arr[a]
        kum += contrib
        print(f"  {a:3d}  {w_arr[a]-w_bar:+12.4e}  {alpha_a[a]:+10.4e}  "
              f"{c_arr[a]:+10.4e}  {contrib:+12.4e}  {kum:+12.4e}", flush=True)
    print(f"  ...", flush=True)
    rest = term2 - kum
    print(f"  Rest (n≥8):  {rest:+12.4e}", flush=True)
    print(f"  Summe:       {term2:+12.4e}", flush=True)

    # --- Secular-Identität ---
    # 1 + C̃ Σ α_a²/(w_a - w_min) = 0
    m_H_wmin = sum(alpha_a[a]**2 / (w_arr[a] - w_min) for a in noncl_A)
    secular_check = 1 + Ctilde * m_H_wmin
    print(f"\n--- Secular-Identität ---", flush=True)
    print(f"  1 + C̃·m_H(w_min) = {secular_check:.6e}  (sollte ≈ 0)", flush=True)
    print(f"  C̃ = -1/m_H(w_min) = {-1/m_H_wmin:.6e}  (vs C̃ = {Ctilde:.6e})", flush=True)

    # --- Rayleigh-Gap physikalisch ---
    print(f"\n--- Rayleigh-Gap: Woher kommt w̄-w_min? ---", flush=True)
    print(f"  w̄ = ⟨ũ, Aũ⟩ = C̃ + ⟨ũ,Hũ⟩", flush=True)
    print(f"  w_min = kleinster EW von A = C̃|ũ⟩⟨ũ| + H", flush=True)
    print(f"  Eigenvektor q₀ zu w_min hat Cluster-Overlap:", flush=True)
    alpha_0 = alpha_a[0]
    print(f"    α₀ = ⟨ũ, q₀⟩ = {alpha_0:.6f}", flush=True)
    print(f"    1 - α₀² = {1 - alpha_0**2:.6e}", flush=True)
    # w̄ - w_min hängt mit der "Nicht-Alignment" von ũ und q₀ zusammen
    # Über Perturbationstheorie: w̄ - w_min ≈ Σ_{a≠0} α_a² (w_a - w_min) / (1) (grob)
    pert_sum = sum(alpha_a[a]**2 * (w_arr[a] - w_min) for a in noncl_A)
    print(f"    Sigma_off alpha_a^2(w_a-w_min) = {pert_sum:.6e}", flush=True)
    print(f"    w̄ - w_min                 = {w_bar - w_min:.6e}", flush=True)
    print(f"    Ratio pert/gap            = {pert_sum/(w_bar-w_min):.4f}", flush=True)

    # --- Zusammenfassung ---
    print(f"\n{'='*80}", flush=True)
    print(f"ZUSAMMENFASSUNG", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"  μ/α = {mu_over_alpha:.6e}", flush=True)
    print(f"  = (w̄-w_min) + ⟨f_λ,k_⊥⟩/α", flush=True)
    print(f"  = {w_bar-w_min:.6e} + {term2/alpha:.6e}", flush=True)
    cancellation_ratio = abs(mu_over_alpha) / (abs(w_bar - w_min) + abs(term2/alpha))
    print(f"  Cancellation-Ratio: {cancellation_ratio:.4f}", flush=True)
    if abs(term1) > abs(term2):
        print(f"  ⇒ RAYLEIGH-GAP DOMINIERT (Term 1 = {abs(term1)/abs(term2):.1f}× Term 2)", flush=True)
    else:
        print(f"  ⇒ POISSON-ALIGNMENT DOMINIERT (Term 2 = {abs(term2)/abs(term1):.1f}× Term 1)", flush=True)

    if np.sign(term1) != np.sign(term2):
        print(f"  ⇒ Beide Terme haben GEGENSÄTZLICHES Vorzeichen → Cancellation", flush=True)
        print(f"  ⇒ μ/α klein durch DESTRUKTIVE INTERFERENZ", flush=True)
    else:
        print(f"  ⇒ Beide Terme haben GLEICHES Vorzeichen → Addition", flush=True)
        print(f"  ⇒ μ/α direkt klein weil beide Terme klein", flush=True)


if __name__ == "__main__":
    scalar_residual_diagnostic(3.0, 30)
    print("\nDone.", flush=True)
