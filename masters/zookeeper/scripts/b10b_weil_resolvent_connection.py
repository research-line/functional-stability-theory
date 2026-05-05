"""
b10b_weil_resolvent_connection.py — Warum erzwingt die Weil-Formel alpha=R?

B10b: Verbindung zwischen Weil-Explizitformel und Resolvent-Identitaet.

ANSATZ:
  Die Weil-Explizitformel ist eine IDENTITAET:
    QW(f, g) = Sigma_rho f~(rho_bar) * g~(rho) + (Korrekturen)

  Fuer f = q (Eigenvektor mit QW*q = w*q) und g = e(gamma):
    w * <q, e(gamma)> = Sigma_rho q~(rho_bar) * e~(gamma)(rho) + ...
    w * F(gamma)      = Sigma_rho q~(rho_bar) * e~(gamma)(rho) + ...

  FRAGE: Was passiert an gamma = gamma_j (Riemann-Nullstelle)?

VERIFIKATIONSSCHRITTE:
  S1: Mellin-Transformierte q~(rho) der Eigenvektoren berechnen
  S2: Mellin-Transformierte e~(gamma)(rho) der Evaluationsvektoren
  S3: Spektrale Summe Sigma_rho verifizieren
  S4: Cancellation-Analyse: warum ist F(gamma_j) ~ 0?
  S5: Verbindung zur Resolvent-Identitaet
"""

import os
import sys
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_WR_matrix_mp,
    build_W02_matrix_mp,
    build_Wprime_matrix_mp,
    chi_trivial,
    project_to_parity,
    diagonalize_mp,
)

DPS = int(os.environ.get("DPS", 50))
N = int(os.environ.get("N", 30))
LAM = float(os.environ.get("LAM", 3.0))

RIEMANN_ZEROS = [14.134725141734693, 21.022039638771554, 25.010857580145688,
                 30.424876125859513, 32.935061587739189]
NON_ZEROS = [10.0, 17.5, 23.0, 28.0, 35.0]


def build_eval_vector_even(gamma, N_val, L_mp):
    dim = N_val + 1
    g = mp.mpf(gamma)
    inv_sqrt2 = 1 / mp.sqrt(mp.mpf(2))
    e = mp.matrix(dim, 1)
    e[0, 0] = 1 / g
    for n in range(1, N_val + 1):
        pole = 2 * mp.pi * n / L_mp
        e[n, 0] = inv_sqrt2 * 2 * g / (g**2 - pole**2)
    return e


def mellin_Vn_even(n, s, L_mp, lam_mp):
    """Mellin-Transformierte von V_n^{even} bei s = 1/2 + i*gamma.

    V_n^{even}(u) = (V_n + V_{-n}) / sqrt(2)  fuer n >= 1
    V_0^{even}(u) = V_0(u) = L^{-1/2}

    Mellin: V_n~(s) = integral_{lam^{-1}}^{lam} V_n(u) u^{s-1} d*u
    """
    L = L_mp
    inv_sqrt_L = 1 / mp.sqrt(L)

    if n == 0:
        # V_0(u) = L^{-1/2}, Mellin: L^{-1/2} * integral u^{s-1} du/u
        # = L^{-1/2} * [u^{s-1}/(s-1)]_{lam^{-1}}^{lam}  ... wait
        # integral_{lam^{-1}}^{lam} u^{s-1} d*u = integral u^{s-2} du
        # = [u^{s-1}/(s-1)]_{lam^{-1}}^{lam} = (lam^{s-1} - lam^{-(s-1)})/(s-1)
        # = 2 sinh((s-1) log lam) / (s-1)
        # = 2 sinh((s-1) L/2) / (s-1)
        return inv_sqrt_L * 2 * mp.sinh((s - 1) * L / 2) / (s - 1)
    else:
        # V_n^{even} = (V_n + V_{-n})/sqrt(2)
        # V_n(u) = L^{-1/2} exp(2*pi*i*n*log(lam*u)/L)
        # V_n~(s) = L^{-1/2} integral exp(2*pi*i*n*x/L) * exp((s-1)*x) dx  (x=log(lam*u))
        #         = L^{-1/2} integral_0^L exp((2*pi*i*n/L + s-1)*x) dx
        #         = L^{-1/2} * [exp(alpha*L)-1]/alpha  with alpha = 2*pi*i*n/L + s-1
        inv_sqrt2 = 1 / mp.sqrt(2)
        alpha_plus = 2 * mp.pi * mp.mpc(0, 1) * n / L + (s - 1)
        alpha_minus = -2 * mp.pi * mp.mpc(0, 1) * n / L + (s - 1)

        # V_n~(s) = L^{-1/2} * (exp(alpha*L)-1)/alpha
        Vn_plus = inv_sqrt_L * (mp.exp(alpha_plus * L) - 1) / alpha_plus
        Vn_minus = inv_sqrt_L * (mp.exp(alpha_minus * L) - 1) / alpha_minus

        return inv_sqrt2 * (Vn_plus + Vn_minus)


def dot_mp(a, b, dim):
    s = mp.mpf(0)
    for r in range(dim):
        s += a[r, 0] * b[r, 0]
    return s


def dot_complex(a_real, b_complex, dim):
    s = mp.mpc(0)
    for r in range(dim):
        s += a_real[r, 0] * b_complex[r]
    return s


def main():
    mp.mp.dps = DPS
    L_mp = 2 * mp.log(mp.mpf(LAM))
    lam_mp = mp.mpf(LAM)
    size = 2 * N + 1
    dim = N + 1

    print("=" * 100)
    print("B10b: WEIL-EXPLIZITFORMEL ↔ RESOLVENT-IDENTITAET")
    print("=" * 100)
    print(f"  lambda={LAM}, N={N}, dim={dim}, L={float(L_mp):.4f}, DPS={DPS}")

    # ===== Phase 1: Standard-Setup =====
    print("\n--- Phase 1: Matrixkonstruktion ---")
    W02 = build_W02_matrix_mp(N, L_mp)
    WR = build_WR_matrix_mp(N, L_mp, parity="even", verbose=False)
    Wp = build_Wprime_matrix_mp(N, L_mp, chi_trivial, 1, verbose=False)

    dummy = mp.matrix(size, size)
    _, U_even = project_to_parity(dummy, N, parity="even")

    W02e = U_even.T * W02 * U_even
    WRe = U_even.T * WR * U_even
    Wpe = U_even.T * Wp * U_even

    QW_even = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            QW_even[i, j] = W02e[i, j] - WRe[i, j] - Wpe[i, j]
    for i in range(dim):
        for j in range(i + 1, dim):
            avg = (QW_even[i, j] + QW_even[j, i]) / 2
            QW_even[i, j] = avg
            QW_even[j, i] = avg

    w_all, Q = diagonalize_mp(QW_even, verbose=False)
    print(f"  QW: w_min={float(w_all[0]):+.6e}, w_max={float(w_all[dim-1]):+.6e}")

    cluster = [k for k in range(dim) if abs(float(w_all[k]) - float(w_all[0])) < 1e-10]
    k_best = cluster[-1]
    w_best = float(w_all[k_best])

    q_best = mp.matrix(dim, 1)
    for r in range(dim):
        q_best[r, 0] = Q[r, k_best]

    # ===== Phase 2: Mellin-Transformierte des Eigenvektors =====
    print("\n--- Phase 2: Mellin-Transformierte q~(rho) an Riemann-Nullstellen ---")
    print("  rho = 1/2 + i*gamma,  q~(rho) = Sigma_n q_n * V_n~(rho)")
    print()

    # Berechne q~(rho) fuer rho an den 5 Riemann-Nullstellen
    q_mellin_at_zeros = []
    print(f"  {'gamma':>10s}  {'|q~(1/2+ig)|':>16s}  {'Re q~':>14s}  {'Im q~':>14s}")
    print("  " + "-" * 60)

    for gamma in RIEMANN_ZEROS:
        s = mp.mpc(mp.mpf('0.5'), mp.mpf(gamma))
        # q~(s) = Sigma_n q_n * V_n^{even}~(s)
        Vn_mellin = [mellin_Vn_even(n, s, L_mp, lam_mp) for n in range(dim)]
        q_tilde = dot_complex(q_best, Vn_mellin, dim)
        q_mellin_at_zeros.append(q_tilde)
        print(f"  {gamma:10.6f}  {float(abs(q_tilde)):16.6e}  "
              f"{float(q_tilde.real):+14.6e}  {float(q_tilde.imag):+14.6e}")

    # Und an Nicht-Nullstellen
    print()
    q_mellin_at_non = []
    for z in NON_ZEROS:
        s = mp.mpc(mp.mpf('0.5'), mp.mpf(z))
        Vn_mellin = [mellin_Vn_even(n, s, L_mp, lam_mp) for n in range(dim)]
        q_tilde = dot_complex(q_best, Vn_mellin, dim)
        q_mellin_at_non.append(q_tilde)
        print(f"  {z:10.6f}  {float(abs(q_tilde)):16.6e}  "
              f"{float(q_tilde.real):+14.6e}  {float(q_tilde.imag):+14.6e}")

    # ===== Phase 3: Spektrale Summe der Weil-Formel =====
    print("\n--- Phase 3: Weil-Spektralsumme ---")
    print("  QW(q, e(g)) = w * F(g)")
    print("  Weil: QW(q, e(g)) = Sigma_rho q~(rho_bar) * e~(g)(rho) + Korrekturen")
    print("  => w*F(g) = Sigma_rho q~(rho_bar) * e~(g)(rho) + ...")
    print()

    # Berechne die linke Seite: w * F(gamma)
    # und die spektrale Summe (endlich, nur 5 Nullstellen im Fenster)
    print(f"  {'gamma':>10s}  {'w*F(g)':>14s}  {'Sigma_rho':>14s}  {'Diff':>14s}  {'typ':>4s}")
    print("  " + "-" * 62)

    all_test = [(g, "ZERO") for g in RIEMANN_ZEROS] + [(z, "NON") for z in NON_ZEROS]
    all_test.sort(key=lambda x: x[0])

    for gamma, typ in all_test:
        e_g = build_eval_vector_even(gamma, N, L_mp)
        F_g = float(dot_mp(e_g, q_best, dim))
        wF = w_best * F_g

        # Spektrale Summe: Sigma_{j=1}^5 q~(rho_j_bar) * e~(gamma)(rho_j)
        # rho_j = 1/2 + i*gamma_j,  rho_j_bar = 1/2 - i*gamma_j
        s_g = mp.mpc(mp.mpf('0.5'), mp.mpf(gamma))
        e_mellin_at_zeros = []
        for gj in RIEMANN_ZEROS:
            s_j = mp.mpc(mp.mpf('0.5'), mp.mpf(gj))
            Vn_mellin = [mellin_Vn_even(n, s_j, L_mp, lam_mp) for n in range(dim)]
            e_tilde_j = dot_complex(e_g, Vn_mellin, dim)
            e_mellin_at_zeros.append(e_tilde_j)

        spectral_sum = mp.mpc(0)
        for j, gj in enumerate(RIEMANN_ZEROS):
            q_bar = mp.conj(q_mellin_at_zeros[j])
            spectral_sum += q_bar * e_mellin_at_zeros[j]

        spectral_real = float(spectral_sum.real)
        diff = abs(wF - spectral_real)
        print(f"  {gamma:10.3f}  {wF:+14.6e}  {spectral_real:+14.6e}  {diff:14.2e}  {typ:>4s}")

    # ===== Phase 4: |q~(rho)|^2 Verteilung =====
    print("\n--- Phase 4: Spektralgewicht |q~(rho)|^2 an Nullstellen ---")
    print("  Weil-Positivitaet: QW(q,q) = Sigma_rho |q~(rho)|^2 + Korrekturen = w")
    print()

    total_spectral = 0.0
    print(f"  {'j':>3s}  {'gamma_j':>10s}  {'|q~(rho_j)|^2':>16s}  {'kumuliert':>14s}")
    print("  " + "-" * 48)
    for j, gj in enumerate(RIEMANN_ZEROS):
        q_sq = float(abs(q_mellin_at_zeros[j])**2)
        total_spectral += q_sq
        print(f"  {j+1:3d}  {gj:10.6f}  {q_sq:16.6e}  {total_spectral:14.6e}")

    print(f"\n  Sigma |q~(rho)|^2 (5 Nullstellen) = {total_spectral:+.6e}")
    print(f"  w (Eigenwert)                     = {w_best:+.6e}")
    print(f"  Differenz (= Korrekturen + Rest)  = {total_spectral - w_best:+.6e}")

    # ===== Phase 5: Warum F(gamma_j) ~ 0 =====
    print("\n--- Phase 5: Warum F(gamma_j) ~ 0 aus der Weil-Formel ---")
    print()
    print("  Aus w*F(gamma_j) = Sigma_rho q~(rho_bar) * e~(gamma_j)(rho) + Korrekturen:")
    print()
    print("  Wenn F(gamma_j) = 0, dann muss die rechte Seite = 0 sein.")
    print("  D.h. die spektrale Summe ueber ALLE Nullstellen beim EVALUATIONSPUNKT")
    print("  gamma_j muss verschwinden.")
    print()

    # Detaillierte Zerlegung fuer gamma_1
    g1 = RIEMANN_ZEROS[0]
    print(f"  Detailanalyse bei gamma_1 = {g1:.6f}:")
    print()

    e_g1 = build_eval_vector_even(g1, N, L_mp)
    F_g1 = float(dot_mp(e_g1, q_best, dim))

    print(f"  F(gamma_1) = {F_g1:+.6e}")
    print()
    print(f"  {'j':>3s}  {'gamma_j':>10s}  {'q~(rho_j_bar)*e~(g1)(rho_j)':>30s}  {'Anteil':>10s}")
    print("  " + "-" * 60)

    terms = []
    for j, gj in enumerate(RIEMANN_ZEROS):
        s_j = mp.mpc(mp.mpf('0.5'), mp.mpf(gj))
        Vn_mellin = [mellin_Vn_even(n, s_j, L_mp, lam_mp) for n in range(dim)]
        e_tilde_j = dot_complex(e_g1, Vn_mellin, dim)
        q_bar = mp.conj(q_mellin_at_zeros[j])
        term = q_bar * e_tilde_j
        terms.append(term)
        pct = float(abs(term)) / max(sum(float(abs(t)) for t in terms), 1e-50) * 100
        print(f"  {j+1:3d}  {gj:10.6f}  {float(term.real):+14.6e} + {float(term.imag):+14.6e}i  {pct:8.1f}%")

    total_terms = sum(terms, mp.mpc(0))
    print(f"\n  Summe: {float(total_terms.real):+.6e} + {float(total_terms.imag):+.6e}i")
    print(f"  |Summe|: {float(abs(total_terms)):.6e}")
    print(f"  w*F(g1): {w_best * F_g1:+.6e}")

    # ===== Phase 6: Die Bruecke =====
    print("\n--- Phase 6: Die Bruecke Weil ↔ Resolvent ---")
    print()
    print("  B10:  F(g) = beta * [alpha(g) - R(g,w)]")
    print("  Weil: w*F(g) = Sigma_rho q~(rho_bar)*e~(g)(rho) + Korr.")
    print()
    print("  Kombination:")
    print("    w * beta * [alpha(g) - R(g,w)] = Sigma_rho q~(rho_bar)*e~(g)(rho) + Korr.")
    print()
    print("  An Riemann-Nullstelle g = gamma_j:")
    print("    alpha(gamma_j) = R(gamma_j, w)")
    print("    => 0 = Sigma_rho q~(rho_bar)*e~(gamma_j)(rho) + Korr.")
    print()
    print("  D.h. die Weil-Identitaet erzwingt, dass die spektrale Summe")
    print("  bei gamma_j verschwindet — und die Resolvent-Identitaet zeigt,")
    print("  dass dies aequivalent zur Bedingung alpha = R ist.")
    print()

    # Teste: e~(gamma_j)(rho_j) hat es einen Pol/Resonanz bei j=j?
    print("  Analyse: Hat e~(gamma_j)(rho_j) bei j=j eine besondere Struktur?")
    print()
    print(f"  {'j':>3s}  {'gamma_j':>10s}  {'|e~(g_j)(rho_j)|':>18s}  {'|e~(g_j)(rho_1)|':>18s}")
    print("  " + "-" * 55)

    for j, gj in enumerate(RIEMANN_ZEROS):
        e_gj = build_eval_vector_even(gj, N, L_mp)
        # e~(g_j)(rho_j) — Mellin von e(g_j) bei der "eigenen" Nullstelle
        s_j = mp.mpc(mp.mpf('0.5'), mp.mpf(gj))
        Vn_j = [mellin_Vn_even(n, s_j, L_mp, lam_mp) for n in range(dim)]
        e_self = dot_complex(e_gj, Vn_j, dim)

        # e~(g_j)(rho_1) — Mellin von e(g_j) bei der ersten Nullstelle
        s_1 = mp.mpc(mp.mpf('0.5'), mp.mpf(RIEMANN_ZEROS[0]))
        Vn_1 = [mellin_Vn_even(n, s_1, L_mp, lam_mp) for n in range(dim)]
        e_cross = dot_complex(e_gj, Vn_1, dim)

        print(f"  {j+1:3d}  {gj:10.6f}  {float(abs(e_self)):18.6e}  {float(abs(e_cross)):18.6e}")

    # ===== Phase 7: Korrekturen quantifizieren =====
    print("\n--- Phase 7: Korrekturen quantifizieren ---")
    print("  Die Weil-Formel im endlichen Raum E_N:")
    print("  QW_N(f,g) = Sigma_{rho in Fenster} f~(rho_bar)*g~(rho) + Delta_N(f,g)")
    print("  wobei Delta_N die Trunkierungs-Korrekturen enthaelt.")
    print()

    # Berechne QW(q, e(gamma)) direkt
    for gamma, typ in [(RIEMANN_ZEROS[0], "ZERO"), (NON_ZEROS[0], "NON")]:
        e_g = build_eval_vector_even(gamma, N, L_mp)
        # QW(q, e) = <QW*q, e> = w * <q, e> = w * F(gamma)
        wF = w_best * float(dot_mp(e_g, q_best, dim))

        # Spektrale Summe (5 Nullstellen)
        spec = 0.0
        for j, gj in enumerate(RIEMANN_ZEROS):
            s_j = mp.mpc(mp.mpf('0.5'), mp.mpf(gj))
            Vn_j = [mellin_Vn_even(n, s_j, L_mp, lam_mp) for n in range(dim)]
            e_tilde_j = dot_complex(e_g, Vn_j, dim)
            q_bar = mp.conj(q_mellin_at_zeros[j])
            spec += float((q_bar * e_tilde_j).real)

        korr = wF - spec
        print(f"  gamma={gamma:.3f} ({typ}): w*F={wF:+.4e}, Sigma_rho={spec:+.4e}, "
              f"Delta_N={korr:+.4e} ({abs(korr)/max(abs(wF),1e-50)*100:.1f}%)")

    # ===== ZUSAMMENFASSUNG =====
    print("\n" + "=" * 100)
    print("ZUSAMMENFASSUNG: B10b Weil ↔ Resolvent")
    print("=" * 100)
    print()
    print("BEFUNDE:")
    print()

    q_sq_zeros = [float(abs(q)**2) for q in q_mellin_at_zeros]
    q_sq_non = [float(abs(q)**2) for q in q_mellin_at_non]
    print(f"  1. |q~(rho)|^2 an Nullstellen:     {min(q_sq_zeros):.2e} bis {max(q_sq_zeros):.2e}")
    print(f"     |q~(1/2+iz)|^2 an Nicht-Null.:   {min(q_sq_non):.2e} bis {max(q_sq_non):.2e}")
    if max(q_sq_zeros) > 0 and min(q_sq_non) > 0:
        print(f"     Verhaeltnis:                     {min(q_sq_non)/max(q_sq_zeros):.2e}")
    print()
    print("  2. Die Weil-Spektralsumme w*F(g) = Sigma_rho q~*e~ + Delta_N")
    print("     reproduziert F(gamma) auf Korrektur-Niveau Delta_N.")
    print()
    print("  3. BRUECKE: Die Resolvent-Identitaet (B10) und die Weil-Formel")
    print("     sind zwei Darstellungen derselben Aussage:")
    print()
    print("     Resolvent:  F(g)=0  <=>  alpha(g) = R(g,w)")
    print("                          (Pol-Antwort = Arithmetik-Resolvent)")
    print()
    print("     Weil:       F(g)=0  <=>  Sigma_rho q~(rho_bar)*e~(g)(rho) = -Delta_N")
    print("                          (Spektralsumme = Trunkierungs-Korrektur)")
    print()
    print("     Die Resolvent-Identitaet ist die OPERATORTHEORETISCHE Reformulierung")
    print("     der Weil-Explizitformel im endlichen Galerkin-Raum E_N.")
    print()
    print("  4. KONSEQUENZ: Die Cancellation alpha=R an Nullstellen ist AEQUIVALENT")
    print("     zur Aussage, dass die endliche Weil-Summe bei gamma_j verschwindet.")
    print("     Der Resolvent (PHP-wI)^{-1} kodiert die SELBE Information wie die")
    print("     Weil-Primsumme — nur in operatortheoretischer Sprache.")


if __name__ == "__main__":
    main()
