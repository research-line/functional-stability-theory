"""
milestone_2_chi_eigenvalue_gap.py
=================================

Milestone 2_chi: Min-Max-Vergleich der unteren Eigenwerte von QW_{lambda,chi}
und PW_lambda in der Prolate-Galerkin-Basis.

Aus DIRICHLET_CCM_TRANSFER.md Abschnitt 6:
  "Im Finit-dim E_N^chi (Galerkin-Einschraenkung): min-max-Vergleich der ersten O(1)
   Eigenwerte von QW_{lambda,chi} und PW_{lambda,chi}. Gibt obere Schranke fuer
   epsilon_{lambda,chi} am unteren Spektrumrand."

Strategie:
  1. Baue Prolate-Basis (via v2 sinc-Diagonalisierung).
  2. Baue PW_lambda (chi-unabhaengig) und QW_{lambda,chi} via Sonin-Zerlegung
     mit parity-korrektem Gamma_arch (Fix gegenueber v2).
  3. Eigendekomposition beider Operatoren in Prolate-Basis.
  4. Vergleiche die ersten 5 Eigenwerte fuer jeden Charakter.
  5. Teste Courant-Fischer-Schranke |lambda_k(QW) - lambda_k(PW)| <= ||Gamma_chi||.

Besonderheit: die unteren Eigenwerte von QW_{lambda,chi} sollten strukturell
die niedrigen L-Nullstellen gamma^(k)(chi) reflektieren, wenn die Sonin-
Zerlegung die Weil-Explizitformel korrekt traegt.

Autor: LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung -- Milestone 2_chi)
Datum: 2026-04-18
Ausfuehrung:
    PYTHONIOENCODING=utf-8 python milestone_2_chi_eigenvalue_gap.py
"""

import json
import numpy as np
import time
from pathlib import Path
from scipy.special import digamma
from typing import Callable, Dict, Any, List

# Wiederverwendbare Setup-Konstanten (analog zu v2)
LAMBDA = np.sqrt(14.0)
L = np.log(LAMBDA)
N_GRID = 600
T_WIDE = 3.0 * L
T_PW = 0.95 * L
N_GALERKIN = 40       # groesser als v2 fuer bessere Eigenwert-Aufloesung
N_L_TERMS = 400
K_EIGENVALUES = 8     # Anzahl unterer Eigenwerte zum Vergleich

RESULTS_DIR = Path(__file__).parent.parent / "_results"
RESULTS_DIR.mkdir(exist_ok=True)

# =================================================================
# 1. Prolate-Basis (analog zu v2)
# =================================================================

def build_prolate_basis(n_grid: int, t_wide: float, t_pw: float,
                         n_galerkin: int) -> Dict[str, np.ndarray]:
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

# =================================================================
# 2. Charaktere (mit korrekter Parity)
# =================================================================

def make_chi(name: str, q: int, parity: int, gamma1: float,
              chi_values: Dict[int, complex]):
    vals = np.zeros(q, dtype=complex)
    for n in range(q):
        if n in chi_values:
            vals[n] = chi_values[n]
    def chi(n):
        return vals[n % q]
    return {'name': name, 'q': q, 'parity': parity, 'chi': chi, 'gamma1': gamma1}

def load_chars() -> List[Dict]:
    atlas_path = (Path(__file__).parent.parent.parent
                   / "dirichlet_atlas" / "_results" / "zeros_all_chars.json")
    zeros = json.load(open(atlas_path))

    chars = []
    # chi_0 (trivial)
    chars.append(make_chi('chi_0', 1, +1, 14.1347, {0: 1.0+0j}))  # dummy, chi(n)=1 immer
    # Spezialbehandlung im Abruf unten (siehe chi-Lambda unten)

    # chi_4: odd, Kronecker(-4)
    chars.append(make_chi('chi_4', 4, -1, 6.0209,
                          {1: 1.0+0j, 3: -1.0+0j}))

    # chi_5: even, Kronecker(+5): chi(1)=1, chi(2)=-1, chi(3)=-1, chi(4)=1
    chars.append(make_chi('chi_5', 5, +1, zeros['chi_5'][0],
                          {1: 1.0+0j, 2: -1.0+0j, 3: -1.0+0j, 4: 1.0+0j}))

    # chi_8: even, Kronecker(+8): chi(1)=1, chi(3)=-1, chi(5)=-1, chi(7)=+1
    chars.append(make_chi('chi_8', 8, +1, zeros['chi_8'][0],
                          {1: 1.0+0j, 3: -1.0+0j, 5: -1.0+0j, 7: 1.0+0j}))

    # chi_33: even, Kronecker(+33), via sympy verifiziert
    chi_33_vals_arr = [0, 1, 1, 0, 1, -1, 0, -1, 1, 0, -1, 0, 0, -1, -1, 0, 1,
                        1, 0, -1, -1, 0, 0, -1, 0, 1, -1, 0, -1, 1, 0, 1, 1]
    chi_33_vals = {n: float(chi_33_vals_arr[n]) + 0j
                    for n in range(33) if chi_33_vals_arr[n] != 0}
    chars.append(make_chi('chi_33', 33, +1, zeros['chi_33'][0], chi_33_vals))

    # chi_0 Fix: chi(n)=1 fuer alle n
    chars[0] = make_chi('chi_0', 1, +1, 14.1347, {})
    def chi_trivial(n): return 1.0 + 0j
    chars[0]['chi'] = chi_trivial

    return chars

# =================================================================
# 3. Operatoren
# =================================================================

def build_PW(basis):
    t = basis['t']
    omega = t**2 + 1.0  # Mellin-Bild PW: t^2 + Konstante
    return project_diagonal(basis['H'], basis['dt'], omega)

def build_Gamma_arch(basis, parity: int):
    """
    Archimedischer Weil-Term mit parity-Korrektur:
        parity = +1 (even): Gamma_arch(t) = 2*Re(psi(1/4 + it/2))
        parity = -1 (odd):  Gamma_arch(t) = 2*Re(psi(3/4 + it/2))
    (CCM 2024 semilokale Trace-Formel; parity-abhaengige Gamma-Faktoren
    aus der Funktionalgleichung.)
    """
    t = basis['t']
    if parity == +1:
        arch = 2.0 * np.real(digamma(0.25 + 1j * t / 2.0))
    else:  # parity == -1
        arch = 2.0 * np.real(digamma(0.75 + 1j * t / 2.0))
    return project_diagonal(basis['H'], basis['dt'], arch)

def build_Gamma_prime(basis, chi: Callable[[int], complex], q: int):
    t = basis['t']
    dt = basis['dt']
    H = basis['H']
    n_grid = len(t)
    p_max = int(LAMBDA ** 2) + 1
    primes = [p for p in range(2, p_max + 1)
              if all(p % q_ != 0 for q_ in range(2, int(np.sqrt(p)) + 1))]
    diff = t[:, None] - t[None, :]
    kernel = np.zeros((n_grid, n_grid), dtype=complex)
    for p in primes:
        if q > 1 and p % q == 0:
            continue
        log_p = np.log(p)
        sqrt_p = np.sqrt(p)
        m_max = max(1, int(np.log(LAMBDA**2) / log_p))
        for m in range(1, m_max + 1):
            chi_pm = chi(p) ** m
            if chi_pm == 0:
                continue
            weight = -chi_pm * log_p / (sqrt_p ** m)
            kernel += weight * np.cos(diff * m * log_p)
    return project_dense(H, dt, kernel)

def build_QW(basis, chi_info, PW_mat):
    G_arch = build_Gamma_arch(basis, chi_info['parity'])
    G_prim = build_Gamma_prime(basis, chi_info['chi'], chi_info['q'])
    return PW_mat + G_arch + G_prim, G_arch, G_prim

# =================================================================
# 4. Eigenwert-Vergleich
# =================================================================

def eigenvalue_gap_analysis(basis, chi_info, PW_mat):
    QW, G_arch, G_prim = build_QW(basis, chi_info, PW_mat)

    # Hermitesch sichern (numerische Symmetrisierung)
    QW_herm = 0.5 * (QW + QW.conj().T)
    PW_herm = 0.5 * (PW_mat + PW_mat.conj().T)

    # Eigenwerte (aufsteigend)
    ew_QW = np.linalg.eigvalsh(QW_herm)
    ew_PW = np.linalg.eigvalsh(PW_herm)

    # Gamma-Groessen
    norm_G_arch = np.linalg.norm(G_arch, ord=2)
    norm_G_prim = np.linalg.norm(G_prim, ord=2)
    norm_G_tot = np.linalg.norm(G_arch + G_prim, ord=2)

    # Min-Max-Check: |lambda_k(QW) - lambda_k(PW)| <= ||G||_2  (Weyl-Ungleichung)
    diffs = ew_QW - ew_PW
    max_abs_diff = np.max(np.abs(diffs))
    weyl_ok = max_abs_diff <= norm_G_tot + 1e-10

    return {
        'chi_name': chi_info['name'],
        'q': chi_info['q'],
        'parity': chi_info['parity'],
        'gamma1': chi_info['gamma1'],
        'ew_QW': ew_QW.tolist(),
        'ew_PW': ew_PW.tolist(),
        'diffs_first_k': diffs[:K_EIGENVALUES].tolist(),
        'norm_G_arch': float(norm_G_arch),
        'norm_G_prim': float(norm_G_prim),
        'norm_G_tot': float(norm_G_tot),
        'max_abs_diff': float(max_abs_diff),
        'weyl_bound_satisfied': bool(weyl_ok),
    }

# =================================================================
# 5. Hauptprogramm
# =================================================================

def main():
    print("=" * 85)
    print("milestone_2_chi_eigenvalue_gap.py")
    print(f"lambda = {LAMBDA:.4f}, L = log(lambda) = {L:.4f}")
    print(f"N_GRID = {N_GRID}, N_GALERKIN = {N_GALERKIN}, K_EIGENVALUES = {K_EIGENVALUES}")
    print(f"Mit parity-korrektem Gamma_arch (FIX gegenueber v2).")
    print("=" * 85)

    t0 = time.time()
    basis = build_prolate_basis(N_GRID, T_WIDE, T_PW, N_GALERKIN)
    print(f"\n[1] Prolate-Basis gebaut ({time.time()-t0:.1f}s).")
    print(f"    Erste 5 Eigenwerte: {basis['lambdas'][:5]}")
    print(f"    Letzter (n={N_GALERKIN-1}): {basis['lambdas'][-1]:.4e}")

    PW_mat = build_PW(basis)
    PW_ew = np.linalg.eigvalsh(0.5*(PW_mat + PW_mat.conj().T))
    print(f"\n[2] PW-Operator gebaut.")
    print(f"    Erste 8 PW-Eigenwerte: {PW_ew[:K_EIGENVALUES]}")

    chars = load_chars()
    print(f"\n[3] {len(chars)} Charaktere geladen.")

    results = []
    for c in chars:
        res = eigenvalue_gap_analysis(basis, c, PW_mat)
        results.append(res)

    # Report
    print("\n" + "=" * 85)
    print("EIGENWERT-VERGLEICH (QW_{lambda,chi} vs PW_lambda)")
    print("=" * 85)
    for res in results:
        print(f"\n{res['chi_name']} (q={res['q']}, parity={res['parity']:+d}, gamma1={res['gamma1']:.3f}):")
        print(f"  QW erste 5 EW: {['{:.4f}'.format(x) for x in res['ew_QW'][:5]]}")
        print(f"  PW erste 5 EW: {['{:.4f}'.format(x) for x in res['ew_PW'][:5]]}")
        print(f"  Diffs       : {['{:+.4f}'.format(x) for x in res['diffs_first_k'][:5]]}")
        print(f"  ||G_arch||  = {res['norm_G_arch']:.4f}")
        print(f"  ||G_prim||  = {res['norm_G_prim']:.4f}")
        print(f"  ||G_tot||   = {res['norm_G_tot']:.4f}")
        print(f"  max|diff|   = {res['max_abs_diff']:.4f}")
        print(f"  Weyl OK?    = {res['weyl_bound_satisfied']}")

    # Zusammenfassung: H1/H2-Check auf unterstem Eigenwert
    print("\n" + "=" * 85)
    print("H1/H2-SIGNATUR AM UNTEREN EIGENWERT")
    print("=" * 85)

    # Baseline chi_0
    base = [r for r in results if r['chi_name'] == 'chi_0'][0]
    base_ew = base['ew_QW'][0]
    base_diff = base['diffs_first_k'][0]
    print(f"\nchi_0 QW-EW[0] = {base_ew:.4f}, diff(QW-PW)[0] = {base_diff:+.4f}")

    print(f"\n{'Charakter':10s} {'gamma1':>8s} {'QW_EW[0]':>10s} {'PW_EW[0]':>10s} "
          f"{'diff[0]':>10s} {'|G_prim|':>10s} {'ratio(G_prim)':>15s}")
    print("-" * 85)
    base_norm = base['norm_G_prim']
    for res in results:
        gp = res['norm_G_prim']
        ratio_gp = gp / max(base_norm, 1e-14)
        print(f"{res['chi_name']:10s} {res['gamma1']:8.3f} {res['ew_QW'][0]:10.4f} "
              f"{res['ew_PW'][0]:10.4f} {res['diffs_first_k'][0]:+10.4f} "
              f"{gp:10.4f} {ratio_gp:15.4f}")

    # H1/H2 auf Even-Familie (chi_5, chi_8, chi_33)
    even_chars = [r for r in results if r['chi_name'] in ('chi_5', 'chi_8', 'chi_33')]
    if len(even_chars) >= 3:
        gs = np.array([r['gamma1'] for r in even_chars])
        # Mehrere Target-Variablen testen
        print(f"\nLog-Log-Regression auf Even-Familie (chi_5, chi_8, chi_33):")
        for metric, values in [
            ('|G_prim|', [r['norm_G_prim'] for r in even_chars]),
            ('|diff[0]|', [abs(r['diffs_first_k'][0]) for r in even_chars]),
            ('|diff[1]|', [abs(r['diffs_first_k'][1]) for r in even_chars]),
            ('QW_EW[0]', [r['ew_QW'][0] for r in even_chars]),
            ('|QW_EW[0]-PW_EW[0]|', [abs(r['ew_QW'][0]-r['ew_PW'][0]) for r in even_chars]),
        ]:
            vals = np.array(values)
            if np.all(vals > 0):
                slope, intercept = np.polyfit(np.log(gs), np.log(vals), 1)
                print(f"  {metric:22s} ~ gamma^{slope:+.3f}   (H1: -3, H2: -1)")
            else:
                print(f"  {metric:22s}: nicht alle positiv, log-reg nicht moeglich")

    # Speichern
    out = RESULTS_DIR / "MILESTONE_2_CHI_EIGENVALUE_2026-04-18.json"
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({
            'config': {
                'lambda': float(LAMBDA), 'L': float(L),
                'N_grid': N_GRID, 'N_galerkin': N_GALERKIN,
                'T_PW': float(T_PW), 'T_wide': float(T_WIDE),
                'K_eigenvalues': K_EIGENVALUES,
            },
            'PW_eigenvalues_first_k': PW_ew[:K_EIGENVALUES].tolist(),
            'results': results,
        }, f, indent=2)
    print(f"\nRohdaten: {out}")

if __name__ == "__main__":
    main()
