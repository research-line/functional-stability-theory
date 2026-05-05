"""
milestone_2_chi_v2_rho.py
=========================

Milestone 2_chi v2: vollstaendige Weil-Form mit rho-Term.

Vorgaenger `milestone_2_chi_eigenvalue_gap.py` hatte die Diagnose:
    QW = PW + Gamma_arch + Gamma_prime      (unvollstaendig)
gibt QW mit *negativen* Eigenwerten fuer nicht-triviale chi.

Korrektur: in der Weil-Quadratform fehlt der **rho-Term** (Summe ueber
L-Nullstellen). Richtige Form aus Weil 1952 / CCM 2024:
    QW_{lambda,chi}(f, g) = sum_{rho} f~(gamma_rho) g~(gamma_rho)
                            - Gamma_arch(chi)
                            - Gamma_prime(chi)

Die ersten Nullstellen werden aus LMFDB-Daten (fuer chi_5, chi_8, chi_33
aus dirichlet_atlas/_results/zeros_all_chars.json) bzw. aus Standard-
Literaturwerten (fuer chi_0 = Riemann, chi_4 = Dirichlet beta) genommen.

Ziel-Tests:
  (a) Wird QW_{lambda,chi} positiv-definit (oder positiv-semidefinit)
      nach Einbau des rho-Terms?
  (b) Zeigt die Defekt-Norm bzw. der Eigenwert-Vergleich eine schaerfere
      H1/H2-Signatur als in der rho-freien Version?

Autor: LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung)
Datum: 2026-04-18
"""

import json
import numpy as np
import time
from pathlib import Path
from scipy.special import digamma
from typing import Callable, Dict, Any, List

# -----------------------------------------------------------------
# 0. Konfiguration
# -----------------------------------------------------------------

LAMBDA = np.sqrt(14.0)
L = np.log(LAMBDA)
# Wichtig: Grid muss breit genug sein, damit die Nullstellen (gamma_k ~ 3-20)
# im Grid-Bereich liegen. Nur im Grid-Bereich wirkt der rho-Term.
T_WIDE = 25.0          # breit genug, 20 Riemann-Nullstellen bis 77 nicht alle, aber
                        # die unterste gamma_1=14.13 und viele weitere drin
N_GRID = 2400           # entsprechend feiner (50x mehr Punkte fuer dt ~0.02)
T_PW = 0.95 * L         # strikte Prolate-Bandbreite (Milestone 1_chi konform)
N_GALERKIN = 40
N_L_TERMS = 400
K_EIGENVALUES = 8

# Cutoff fuer rho-Summe: alle gamma <= T_WIDE (d.h. im Grid) werden beruecksichtigt
RHO_CUTOFF = 20.0       # 20.0 < T_WIDE, aber umfasst alle relevanten Nullstellen

RESULTS_DIR = Path(__file__).parent.parent / "_results"
RESULTS_DIR.mkdir(exist_ok=True)

# Riemann-Nullstellen (ersten 20, aus Odlyzko-Tabellen)
RIEMANN_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005150, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
]

# Dirichlet-Beta-Nullstellen (L(s, chi_4), ersten 10; LMFDB)
CHI4_ZEROS = [
    6.020948, 10.243766, 12.988096, 16.343297, 18.291996,
    21.428273, 23.265376, 26.068044, 28.106108, 30.296575,
]

# -----------------------------------------------------------------
# 1. Prolate-Basis (wiederverwendet)
# -----------------------------------------------------------------

def build_prolate_basis(n_grid, t_wide, t_pw, n_galerkin):
    t = np.linspace(-t_wide, t_wide, n_grid)
    dt = t[1] - t[0]
    diff = t[:, None] - t[None, :]
    K = np.where(np.abs(diff) < 1e-14,
                 t_pw / np.pi,
                 np.sin(t_pw * diff) / (np.pi * diff))
    K *= dt
    K = 0.5 * (K + K.T)
    eigvals, eigvecs = np.linalg.eigh(K)
    idx = np.argsort(eigvals)[::-1]
    H = eigvecs[:, idx][:, :n_galerkin] / np.sqrt(dt)
    lambdas = eigvals[idx][:n_galerkin]
    return {'t': t, 'dt': dt, 'H': H, 'lambdas': lambdas}

def project_diagonal(H, dt, diag):
    return (H.T * diag[None, :]) @ H * dt

def project_dense(H, dt, kernel):
    return H.T @ kernel @ H * dt**2

def interp_H_at(basis, gamma_vals):
    """
    Interpoliere Prolate-Basis-Werte h_n(t=gamma_k) fuer gegebene gamma_k
    via lineare Interpolation auf dem Grid.
    Rueckgabe: (len(gamma_vals), N_galerkin) Matrix.
    """
    t = basis['t']
    H = basis['H']  # (N_grid, N_gal)
    dt = basis['dt']
    # Lineare Interpolation jeder Basis-Funktion
    result = np.zeros((len(gamma_vals), H.shape[1]))
    for i, g in enumerate(gamma_vals):
        if g < t[0] or g > t[-1]:
            continue  # ausserhalb Grid -> 0
        # finde Index
        idx = np.searchsorted(t, g)
        if idx <= 0:
            result[i, :] = H[0, :]
        elif idx >= len(t):
            result[i, :] = H[-1, :]
        else:
            alpha = (g - t[idx-1]) / dt
            result[i, :] = (1 - alpha) * H[idx-1, :] + alpha * H[idx, :]
    return result

# -----------------------------------------------------------------
# 2. Gamma-Terme
# -----------------------------------------------------------------

def build_Gamma_arch(basis, parity):
    t = basis['t']
    if parity == +1:
        arch = 2.0 * np.real(digamma(0.25 + 1j * t / 2.0))
    else:
        arch = 2.0 * np.real(digamma(0.75 + 1j * t / 2.0))
    return project_diagonal(basis['H'], basis['dt'], arch)

def build_Gamma_prime(basis, chi, q):
    t = basis['t']; dt = basis['dt']; H = basis['H']
    n_grid = len(t)
    p_max = int(LAMBDA ** 2) + 1
    primes = [p for p in range(2, p_max + 1)
              if all(p % q_ != 0 for q_ in range(2, int(np.sqrt(p)) + 1))]
    diff = t[:, None] - t[None, :]
    kernel = np.zeros((n_grid, n_grid), dtype=complex)
    for p in primes:
        if q > 1 and p % q == 0:
            continue
        log_p = np.log(p); sqrt_p = np.sqrt(p)
        m_max = max(1, int(np.log(LAMBDA**2) / log_p))
        for m in range(1, m_max + 1):
            chi_pm = chi(p) ** m
            if chi_pm == 0: continue
            weight = -chi_pm * log_p / (sqrt_p ** m)
            kernel += weight * np.cos(diff * m * log_p)
    return project_dense(H, dt, kernel)

def build_Gamma_rho(basis, zeros_list, cutoff):
    """
    rho-Term: projektor-artige Matrix aus Delta-Distributionen an den Nullstellen.

    Gamma_rho[m, n] = sum_{|gamma| <= cutoff} [H(gamma)_m * H(gamma)_n
                                                + H(-gamma)_m * H(-gamma)_n]

    (Symmetrische Erweiterung: Weil-Form nimmt sowohl +gamma als auch -gamma
    da die L-Nullstellen auf der krit. Linie symmetrisch um t=0 liegen.)
    """
    valid = [g for g in zeros_list if abs(g) <= cutoff]
    if not valid:
        return np.zeros((basis['H'].shape[1], basis['H'].shape[1]), dtype=complex)
    gammas = np.array(valid + [-g for g in valid])
    H_at_gamma = interp_H_at(basis, gammas)  # (2K, N_gal)
    # sum over k of outer(H_at_gamma[k], H_at_gamma[k])
    Gamma_rho = H_at_gamma.T @ H_at_gamma
    return Gamma_rho.astype(complex)

# -----------------------------------------------------------------
# 3. Charaktere + Nullstellen
# -----------------------------------------------------------------

def make_chi(name, q, parity, chi_values, zeros):
    vals = np.zeros(q, dtype=complex)
    for n in range(q):
        if n in chi_values:
            vals[n] = chi_values[n]
    def chi(n): return vals[n % q]
    return {'name': name, 'q': q, 'parity': parity, 'chi': chi,
            'zeros': zeros, 'gamma1': zeros[0] if zeros else 0.0}

def load_chars():
    atlas_path = (Path(__file__).parent.parent.parent
                   / "dirichlet_atlas" / "_results" / "zeros_all_chars.json")
    zeros_atlas = json.load(open(atlas_path))

    chars = []
    # chi_0 (trivial, = Riemann zeta)
    c0 = make_chi('chi_0', 1, +1, {}, RIEMANN_ZEROS)
    def chi_triv(n): return 1.0 + 0j
    c0['chi'] = chi_triv
    chars.append(c0)

    # chi_4 (Dirichlet beta, odd)
    chars.append(make_chi('chi_4', 4, -1,
                          {1: 1.0+0j, 3: -1.0+0j},
                          CHI4_ZEROS))

    # chi_5 even
    chars.append(make_chi('chi_5', 5, +1,
                          {1: 1.0+0j, 2: -1.0+0j, 3: -1.0+0j, 4: 1.0+0j},
                          zeros_atlas['chi_5']))

    # chi_8 even
    chars.append(make_chi('chi_8', 8, +1,
                          {1: 1.0+0j, 3: -1.0+0j, 5: -1.0+0j, 7: 1.0+0j},
                          zeros_atlas['chi_8']))

    # chi_33 even
    chi_33_arr = [0, 1, 1, 0, 1, -1, 0, -1, 1, 0, -1, 0, 0, -1, -1, 0, 1,
                   1, 0, -1, -1, 0, 0, -1, 0, 1, -1, 0, -1, 1, 0, 1, 1]
    chi_33_vals = {n: float(chi_33_arr[n]) + 0j
                    for n in range(33) if chi_33_arr[n] != 0}
    chars.append(make_chi('chi_33', 33, +1, chi_33_vals,
                          zeros_atlas['chi_33']))

    return chars

# -----------------------------------------------------------------
# 4. QW-Konstruktion und Analyse
# -----------------------------------------------------------------

def build_QW_full(basis, chi_info, cutoff):
    G_arch = build_Gamma_arch(basis, chi_info['parity'])
    G_prim = build_Gamma_prime(basis, chi_info['chi'], chi_info['q'])
    G_rho = build_Gamma_rho(basis, chi_info['zeros'], cutoff)
    # Weil-Konvention: QW = Gamma_rho - Gamma_arch - Gamma_prime
    QW = G_rho - G_arch - G_prim
    return QW, G_arch, G_prim, G_rho

def analyze(basis, chi_info, PW_mat, cutoff):
    QW, G_arch, G_prim, G_rho = build_QW_full(basis, chi_info, cutoff)
    QW_herm = 0.5 * (QW + QW.conj().T)
    PW_herm = 0.5 * (PW_mat + PW_mat.conj().T)
    ew_QW = np.linalg.eigvalsh(QW_herm)
    ew_PW = np.linalg.eigvalsh(PW_herm)

    n_valid_zeros = sum(1 for g in chi_info['zeros'] if abs(g) <= cutoff)

    return {
        'chi_name': chi_info['name'],
        'q': chi_info['q'],
        'parity': chi_info['parity'],
        'gamma1': chi_info['gamma1'],
        'n_zeros_used': n_valid_zeros,
        'ew_QW_first_8': ew_QW[:K_EIGENVALUES].tolist(),
        'ew_QW_last_3': ew_QW[-3:].tolist(),
        'ew_QW_min': float(ew_QW.min()),
        'ew_QW_max': float(ew_QW.max()),
        'ew_PW_first_8': ew_PW[:K_EIGENVALUES].tolist(),
        'norm_G_arch': float(np.linalg.norm(G_arch, 2)),
        'norm_G_prim': float(np.linalg.norm(G_prim, 2)),
        'norm_G_rho':  float(np.linalg.norm(G_rho, 2)),
        'norm_QW':     float(np.linalg.norm(QW, 2)),
        'n_negative_eigenvalues': int(np.sum(ew_QW < -1e-8)),
        'positivity_gap': float(ew_QW.min()),  # <0 => not PSD
    }

# -----------------------------------------------------------------
# 5. Defekt-Norm mit volle QW
# -----------------------------------------------------------------

def L_values_on_line(t_grid, chi, n_terms=N_L_TERMS):
    result = np.zeros_like(t_grid, dtype=complex)
    s_grid = 0.5 + 1j * t_grid
    for n in range(1, n_terms + 1):
        cn = chi(n)
        if cn == 0: continue
        result += cn / (n ** s_grid)
    return result

def build_Psi(basis, chi):
    L_vals = L_values_on_line(basis['t'], chi)
    return project_diagonal(basis['H'], basis['dt'], L_vals)

def build_PW(basis):
    omega = basis['t']**2 + 1.0
    return project_diagonal(basis['H'], basis['dt'], omega)

def compute_defect(basis, chi_info, PW_mat, cutoff):
    QW, *_ = build_QW_full(basis, chi_info, cutoff)
    Psi = build_Psi(basis, chi_info['chi'])
    A = QW @ Psi
    B = Psi @ PW_mat
    tr_AB = np.trace(A.conj().T @ B)
    tr_BB = np.trace(B.conj().T @ B)
    mu_opt = tr_AB / tr_BB if abs(tr_BB) > 1e-14 else 1.0
    defect = A - mu_opt * B
    def_spec = np.linalg.norm(defect, 2)
    norm_B = np.linalg.norm(B, 2)
    norm_Psi_F = np.linalg.norm(Psi, 'fro')
    return {
        'chi_name': chi_info['name'],
        'gamma1': chi_info['gamma1'],
        'parity': chi_info['parity'],
        'defect_spec': float(def_spec),
        'rel_a_v1style': float(def_spec / max(norm_B, 1e-14)),
        'rel_b_polnorm': float(def_spec / max(norm_Psi_F, 1e-14)),
        'mu_opt_real': float(np.real(mu_opt)),
        'mu_opt_imag': float(np.imag(mu_opt)),
    }

# -----------------------------------------------------------------
# 6. Main
# -----------------------------------------------------------------

def main():
    cutoff = RHO_CUTOFF
    print("=" * 90)
    print("milestone_2_chi_v2_rho.py -- VOLLE Weil-Form: QW = Gamma_rho - Gamma_arch - Gamma_prime")
    print(f"lambda = {LAMBDA:.4f}, L = log(lambda) = {L:.4f}")
    print(f"N_GRID = {N_GRID}, N_GALERKIN = {N_GALERKIN}")
    print(f"rho-Cutoff: |gamma| <= {cutoff:.4f} (T_WIDE = {T_WIDE})")
    print("=" * 90)

    basis = build_prolate_basis(N_GRID, T_WIDE, T_PW, N_GALERKIN)
    PW_mat = build_PW(basis)
    chars = load_chars()

    # Analyse jedes Charakters
    print("\n[Eigenwerte von QW_{lambda,chi} (volle Weil-Form)]")
    print("-" * 90)
    results_ew = []
    results_def = []
    for c in chars:
        res = analyze(basis, c, PW_mat, cutoff)
        results_ew.append(res)
        dres = compute_defect(basis, c, PW_mat, cutoff)
        results_def.append(dres)
        print(f"\n{res['chi_name']} (q={res['q']}, parity={res['parity']:+d}, "
              f"gamma1={res['gamma1']:.3f}, {res['n_zeros_used']} Nullstellen im Cutoff):")
        print(f"  QW erste 5 EW: {['{:.4f}'.format(x) for x in res['ew_QW_first_8'][:5]]}")
        print(f"  QW letzte 3 EW: {['{:.4f}'.format(x) for x in res['ew_QW_last_3']]}")
        print(f"  min EW = {res['ew_QW_min']:+.4f}   (>=0 = positiv-definit)")
        print(f"  {res['n_negative_eigenvalues']} negative EW")
        print(f"  ||G_arch|| = {res['norm_G_arch']:.4f}  ||G_prim|| = {res['norm_G_prim']:.4f}  "
              f"||G_rho|| = {res['norm_G_rho']:.4f}  ||QW|| = {res['norm_QW']:.4f}")

    # Positivitaet-Zusammenfassung
    print("\n" + "=" * 90)
    print("POSITIVITAETS-CHECK")
    print("=" * 90)
    print(f"{'Charakter':10s} {'min EW':>10s} {'#negative':>10s} {'positiv?':>12s}")
    print("-" * 50)
    for res in results_ew:
        posit = "JA" if res['ew_QW_min'] >= -1e-6 else "nein"
        print(f"{res['chi_name']:10s} {res['ew_QW_min']:10.4f} "
              f"{res['n_negative_eigenvalues']:10d} {posit:>12s}")

    # Defekt-Norm mit voller QW
    print("\n" + "=" * 90)
    print("DEFEKT-NORM MIT VOLLER WEIL-FORM")
    print("=" * 90)
    print(f"{'Charakter':10s} {'gamma1':>8s} {'rel_a':>10s} {'rel_b':>10s} {'mu_opt_re':>10s}")
    print("-" * 55)
    for res in results_def:
        print(f"{res['chi_name']:10s} {res['gamma1']:8.3f} "
              f"{res['rel_a_v1style']:10.4f} {res['rel_b_polnorm']:10.4f} "
              f"{res['mu_opt_real']:10.4f}")

    # H1/H2-Log-Log auf Even-Familie
    even = [r for r in results_def if r['chi_name'] in ('chi_5', 'chi_8', 'chi_33')]
    if len(even) == 3:
        gs = np.array([r['gamma1'] for r in even])
        for metric in ('rel_a_v1style', 'rel_b_polnorm'):
            vals = np.array([r[metric] for r in even])
            if np.all(vals > 0):
                slope, _ = np.polyfit(np.log(gs), np.log(vals), 1)
                print(f"\nLog-Log-Slope ({metric} auf Even-Familie): {slope:+.3f}  (H1: -3, H2: -1)")

    # Speichern
    out = RESULTS_DIR / "MILESTONE_2_CHI_V2_RHO_2026-04-18.json"
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({
            'config': {
                'lambda': float(LAMBDA), 'L': float(L),
                'N_grid': N_GRID, 'N_galerkin': N_GALERKIN,
                'rho_cutoff': float(cutoff),
                'n_Riemann_zeros': len(RIEMANN_ZEROS),
                'n_chi4_zeros': len(CHI4_ZEROS),
            },
            'eigenvalue_analysis': results_ew,
            'defect_norms': results_def,
        }, f, indent=2)
    print(f"\nRohdaten: {out}")

if __name__ == "__main__":
    main()
