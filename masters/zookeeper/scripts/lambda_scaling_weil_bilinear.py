"""
lambda_scaling_weil_bilinear.py
===============================

Erster echter lambda-Asymptotik-Test der chi-twisteten Weil-Form.

Zentrale Erkenntnis aus MILESTONE_2_CHI_V2_RHO_2026-04-18.md:
  Die Matrix-Form "QW = G_rho - G_arch - G_prim" ist NICHT die
  Weil-Positivitaets-Aussage. Diese ist bilinear-quadratisch:
    Q_Weil(f) = sum_rho |f~(gamma_rho)|^2 - arch(f,f) - prim(f,f)
  Ziel: Q_Weil(f) >= 0 fuer alle Testfunktionen f im PW_{log(lambda)}-Raum.

Dieses Skript testet die Skalen-Abhaengigkeit:

  Fuer jedes lambda in {sqrt(14), 2*sqrt(14), 4*sqrt(14), 8*sqrt(14)}:
    - Baue PW-Prolate-Basis der Dimension N_galerkin
    - Fuer jede Basisfunktion f = h_n:
        Q_n = sum_rho |f~(gamma_rho)|^2 - arch(f,f) - prim(f,f)
    - Min Q_n ueber alle n = "Positivitaets-Indikator"
    - Summe Q_n / N = "Mittlerer Defekt"
  Teste:
    (a) Wird min(Q_n) >= 0 (oder wenigstens monoton ansteigend) mit lambda?
    (b) Skaliert der mittlere Defekt wie 1/lambda, oder anders?

Besonders getestet: chi_21 (gamma1=2.32, wird bei lambda~e^2.3=10 sichtbar)
und chi_33 (gamma1=3.0, wird bei lambda~e^3=20 sichtbar).

Autor: LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung 3)
Datum: 2026-04-18
Ausfuehrung:
    PYTHONIOENCODING=utf-8 python lambda_scaling_weil_bilinear.py
"""

import json
import numpy as np
import time
from pathlib import Path
from scipy.special import digamma

# -----------------------------------------------------------------
# 0. Skalen-Schema
# -----------------------------------------------------------------

LAMBDA_VALUES = [
    np.sqrt(14.0),          # ~3.74
    2.0 * np.sqrt(14.0),    # ~7.48
    4.0 * np.sqrt(14.0),    # ~14.97
    8.0 * np.sqrt(14.0),    # ~29.93
]

N_GALERKIN = 30
N_L_TERMS = 400
RHO_CUTOFF = 40.0          # nimm alle Nullstellen bis 40

RESULTS_DIR = Path(__file__).parent.parent / "_results"

# Mehr Riemann-Nullstellen (30 Stueck)
RIEMANN_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005150, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
    79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
    92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
]
CHI4_ZEROS = [
    6.020948, 10.243766, 12.988096, 16.343297, 18.291996,
    21.428273, 23.265376, 26.068044, 28.106108, 30.296575,
    31.765775, 34.479648, 35.872660, 37.967690, 40.340230,
]

# -----------------------------------------------------------------
# 1. Prolate-Basis bei Skala lambda
# -----------------------------------------------------------------

def build_basis(lam, n_galerkin, rho_cutoff):
    """
    Baue Prolate-Basis mit Paley-Wiener-Bandbreite T_PW = 0.95*log(lambda),
    Grid breit genug, um Nullstellen bis rho_cutoff zu enthalten.
    """
    L_val = np.log(lam)
    T_PW = 0.95 * L_val
    # Grid: mindestens so breit wie rho_cutoff + Sicherheitsrand
    T_wide = max(rho_cutoff + 5.0, 3.0 * L_val)
    # Grid-Feinheit: dt klein genug, dass bei T_PW*dt << 1 sinc aufloest
    # Aber auch dt << 1/max(gamma), damit Interpolation praezise
    n_grid = min(3000, max(1200, int(40 * T_wide / L_val)))
    t = np.linspace(-T_wide, T_wide, n_grid)
    dt = t[1] - t[0]
    diff = t[:, None] - t[None, :]
    K = np.where(np.abs(diff) < 1e-14,
                 T_PW / np.pi,
                 np.sin(T_PW * diff) / (np.pi * diff))
    K *= dt
    K = 0.5 * (K + K.T)
    ev, evv = np.linalg.eigh(K)
    idx = np.argsort(ev)[::-1]
    H = evv[:, idx][:, :n_galerkin] / np.sqrt(dt)
    lambdas = ev[idx][:n_galerkin]
    return {'t': t, 'dt': dt, 'H': H, 'lambdas': lambdas,
            'T_PW': T_PW, 'T_wide': T_wide, 'L': L_val, 'n_grid': n_grid,
            'lambda_val': lam}

# -----------------------------------------------------------------
# 2. Bilineare Weil-Form
# -----------------------------------------------------------------

def mellin_transform_at(basis, gamma_vals):
    """
    Fuer jede Basisfunktion h_n, wert h_n(t=gamma) via Interpolation.
    Das ist im Grunde f~(gamma) bis auf Normalisierung.
    """
    t = basis['t']; H = basis['H']; dt = basis['dt']
    n_grid, n_gal = H.shape
    result = np.zeros((len(gamma_vals), n_gal))
    for i, g in enumerate(gamma_vals):
        if g < t[0] or g > t[-1]:
            continue
        idx = np.searchsorted(t, g)
        if idx <= 0:
            result[i, :] = H[0, :]
        elif idx >= n_grid:
            result[i, :] = H[-1, :]
        else:
            alpha = (g - t[idx-1]) / dt
            result[i, :] = (1-alpha)*H[idx-1,:] + alpha*H[idx,:]
    return result

def arch_form(basis, parity, q):
    """
    Weil-Archimedischer Term mit Conductor-Korrektur fuer L(s, chi):

      arch(h, h) = - log(q/pi) * |h~(0)|^2
                   + int h(t) [2*Re digamma(arg + it/2)] h(t) dt

    Der log(q/pi)-Term ist der Conductor-Beitrag (fehlte in frueherer Version).
    Ohne diesen Term ist die Weil-Form fuer nicht-triviale chi nicht korrekt
    ausbalanciert.

    Vorzeichen-Konvention: der arch-Term wird in Q = rho - arch - prim
    SUBTRAHIERT, also ein positiver log(q/pi) macht die Form positiver
    (= guenstiger fuer Weil-Positivitaet).
    """
    t = basis['t']
    arg = 0.25 if parity == +1 else 0.75
    v = 2.0 * np.real(digamma(arg + 1j*t/2.0))
    # Diagonal-Integral
    arch_matrix = (basis['H'].T * v[None, :]) @ basis['H'] * basis['dt']
    # Conductor-Term: log(q/pi) skaliert mit Rang-1-Projektion auf h~(0)
    # Naeherung: h~(0) = int h(t) dt, also h_n-Summe * dt
    h_0 = np.sum(basis['H'], axis=0) * basis['dt']  # (N_gal,) = h_n~(t=0)
    # Das gehoert eigentlich zu einer Testfunktion Mellin-Transform an 1/2
    # In unserer t-Koordinate (Im(s)), t=0 entspricht Re(s)=1/2, Im(s)=0, also s=1/2.
    # h_n~(s=1/2) = int h_n(t) e^{-(1/2) * log ?} dt
    # Fuer Einfachheit: nehme h_n(t=0)-Wert (diskret)
    # Aber h_n im Mellin-Bild, t=0 entspricht Re(s)=1/2, Im(s)=0.
    # Der log(q/pi)-Term traegt log(q/pi) * |int h_n|^2
    conductor_offset = np.log(q / np.pi)
    conductor_matrix = conductor_offset * np.outer(h_0, h_0)
    # Vorzeichen-Konvention Weil: arch-Term enthaelt -log(q/pi) * |h~(1/2)|^2
    # d.h. fuer grosses q: arch wird negativer (in der Subtraktions-Konvention),
    # was die Gesamtform positiver macht.
    # Hier: arch(h,h) explizit inkl. -conductor_offset
    return arch_matrix - conductor_matrix

def prime_form(basis, chi, q, lam):
    """
    prim(h_m, h_n) = int int h_m(t) K_prim(t,t') h_n(t') dt dt'
    Diese ist die Projektion des Weil-Prim-Kerns.
    """
    t = basis['t']; dt = basis['dt']; H = basis['H']
    n_grid = len(t)
    p_max = int(lam**2) + 1
    primes = [p for p in range(2, p_max+1)
              if all(p % q_ != 0 for q_ in range(2, int(np.sqrt(p))+1))]
    diff = t[:, None] - t[None, :]
    K = np.zeros((n_grid, n_grid), dtype=complex)
    for p in primes:
        if q > 1 and p % q == 0: continue
        lp = np.log(p); sp = np.sqrt(p)
        m_max = max(1, int(np.log(lam**2) / lp))
        for m in range(1, m_max+1):
            cpm = chi(p) ** m
            if cpm == 0: continue
            w = -cpm * lp / (sp**m)
            K += w * np.cos(diff * m * lp)
    return H.T @ K @ H * dt**2

def rho_form_diagonal(basis, zeros, cutoff):
    """
    Q_rho(h_n, h_n) = sum_{gamma} |h_n~(gamma)|^2 + |h_n~(-gamma)|^2
    als diagonale Werte (n-te Basisfunktion).

    Rueckgabe: array of shape (N_galerkin,) mit Q_rho-Diagonalwerten.

    Im Matrix-Sinn: Q_rho[m,n] = sum_gamma h_m(gamma)*h_n(gamma) + h_m(-gamma)*h_n(-gamma)
    Das ist ein Rank-K-Operator, die Diagonalwerte sind direkt die Beitraege.
    """
    valid = [g for g in zeros if abs(g) <= cutoff]
    if not valid:
        return np.zeros(basis['H'].shape[1])
    gammas = np.array(valid + [-g for g in valid])
    Hg = mellin_transform_at(basis, gammas)  # (2K, N_gal)
    # Diagonal: sum_k |Hg[k,n]|^2
    diag = np.sum(Hg * Hg, axis=0)
    return diag

# -----------------------------------------------------------------
# 3. Q_Weil-Testroutine
# -----------------------------------------------------------------

def test_weil_bilinear(basis, chi_info):
    """
    Fuer jede Prolate-Basisfunktion h_n:
      Q_n = Q_rho(h_n, h_n) - arch(h_n, h_n) - prim(h_n, h_n)

    Wenn Q_n >= 0 fuer alle n, ist Weil-Positivitaet auf diesem Raum erfuellt.
    Andernfalls: Defekt-Mass = min_n Q_n.
    """
    arch_mat = arch_form(basis, chi_info['parity'], chi_info['q'])
    prim_mat = prime_form(basis, chi_info['chi'], chi_info['q'], basis['lambda_val'])
    rho_diag = rho_form_diagonal(basis, chi_info['zeros'], RHO_CUTOFF)

    # Diagonal-Werte (Q_n fuer jede Basisfunktion)
    arch_diag = np.real(np.diag(arch_mat))
    prim_diag = np.real(np.diag(prim_mat))
    # Wichtig: Vorzeichen von arch/prim ist Konvention.
    # Weil-Formel: Q_Weil(f) = rho-Term - arch - prim
    Q = rho_diag - arch_diag - prim_diag

    return {
        'chi_name': chi_info['name'],
        'gamma1': chi_info['gamma1'],
        'Q_min': float(Q.min()),
        'Q_max': float(Q.max()),
        'Q_mean': float(Q.mean()),
        'Q_positive_fraction': float(np.mean(Q >= 0)),
        'n_positive': int(np.sum(Q >= 0)),
        'n_total': len(Q),
        'rho_diag_first_5': rho_diag[:5].tolist(),
        'arch_diag_first_5': arch_diag[:5].tolist(),
        'prim_diag_first_5': prim_diag[:5].tolist(),
        'Q_first_5': Q[:5].tolist(),
    }

# -----------------------------------------------------------------
# 4. Charaktere
# -----------------------------------------------------------------

def load_chars():
    zf = Path(__file__).parent.parent.parent / "dirichlet_atlas" / "_results" / "zeros_all_chars.json"
    za = json.load(open(zf))

    def mkchi(name, q, parity, vals, zeros):
        arr = np.zeros(q, dtype=complex)
        for n, v in vals.items(): arr[n] = v
        def chi(n): return arr[n % q]
        return {'name': name, 'q': q, 'parity': parity, 'chi': chi, 'zeros': zeros,
                'gamma1': zeros[0] if zeros else 0}

    chars = []
    # chi_0
    c0 = mkchi('chi_0', 1, +1, {}, RIEMANN_ZEROS)
    def chi_t(n): return 1.0 + 0j
    c0['chi'] = chi_t
    chars.append(c0)
    # chi_4 odd
    chars.append(mkchi('chi_4', 4, -1, {1: 1.0+0j, 3: -1.0+0j}, CHI4_ZEROS))
    # chi_5 even
    chars.append(mkchi('chi_5', 5, +1,
                       {1: 1.0+0j, 2: -1.0+0j, 3: -1.0+0j, 4: 1.0+0j},
                       za['chi_5']))
    # chi_21 even (gamma1=2.32 - niedrigste!)
    # Kronecker(21, n) via sympy verifiziert: 21 = 3*7
    chi_21_arr = [0, 1, -1, 0, 1, 1, 0, 0, -1, 0, -1, -1, 0, 1, 0, 0, 1, 1, 0, -1, 1]
    chars.append(mkchi('chi_21', 21, +1,
                       {n: float(chi_21_arr[n])+0j for n in range(21) if chi_21_arr[n] != 0},
                       za['chi_21']))
    # chi_33 even
    arr33 = [0,1,1,0,1,-1,0,-1,1,0,-1,0,0,-1,-1,0,1,1,0,-1,-1,0,0,-1,0,1,-1,0,-1,1,0,1,1]
    chars.append(mkchi('chi_33', 33, +1,
                       {n: float(arr33[n])+0j for n in range(33) if arr33[n] != 0},
                       za['chi_33']))
    return chars

# -----------------------------------------------------------------
# 5. Main
# -----------------------------------------------------------------

def main():
    print("=" * 90)
    print("lambda_scaling_weil_bilinear.py")
    print("Bilineare Weil-Form Q(f) = sum |f~(gamma)|^2 - arch(f,f) - prim(f,f)")
    print(f"N_GALERKIN = {N_GALERKIN}, rho_cutoff = {RHO_CUTOFF}")
    print("=" * 90)

    chars = load_chars()
    print(f"\nCharaktere: {[c['name'] for c in chars]}")
    g1_strs = ['{:.2f}'.format(c['gamma1']) for c in chars]
    print(f"gamma^(1): {g1_strs}")

    all_results = []
    for lam in LAMBDA_VALUES:
        L_val = np.log(lam)
        print(f"\n{'-'*90}")
        print(f"lambda = {lam:.4f}  (log lambda = {L_val:.4f})")
        print(f"{'-'*90}")
        t0 = time.time()
        basis = build_basis(lam, N_GALERKIN, RHO_CUTOFF)
        print(f"  Basis: N_grid = {basis['n_grid']}, T_wide = {basis['T_wide']:.2f}, "
              f"T_PW = {basis['T_PW']:.4f}")
        print(f"  Prolate EW erste 3: {basis['lambdas'][:3]}, letzte: {basis['lambdas'][-1]:.2e}")
        print(f"  Build-Zeit: {time.time()-t0:.1f}s")

        print(f"\n  {'chi':8s} {'gamma1':>7s} {'Q_min':>10s} {'Q_max':>10s} "
              f"{'Q_mean':>10s} {'pos%':>6s} {'rho_1':>9s} {'arch_1':>9s} {'prim_1':>9s}")
        results_this_lambda = []
        for c in chars:
            t1 = time.time()
            res = test_weil_bilinear(basis, c)
            res['lambda'] = float(lam)
            res['log_lambda'] = float(L_val)
            res['time'] = time.time() - t1
            # Check: gamma^(1) im PW-Bereich?
            visible = c['gamma1'] < basis['T_PW']
            res['gamma1_in_PW'] = bool(visible)
            results_this_lambda.append(res)
            print(f"  {c['name']:8s} {c['gamma1']:7.3f} {res['Q_min']:10.4f} "
                  f"{res['Q_max']:10.4f} {res['Q_mean']:10.4f} "
                  f"{res['Q_positive_fraction']*100:5.1f}% "
                  f"{res['rho_diag_first_5'][0]:9.4f} "
                  f"{res['arch_diag_first_5'][0]:9.4f} "
                  f"{res['prim_diag_first_5'][0]:9.4f}"
                  + ("  <-- gamma1 sichtbar" if visible else ""))
        all_results.append({'lambda': float(lam), 'log_lambda': float(L_val),
                            'basis_info': {'N_grid': basis['n_grid'],
                                           'T_PW': float(basis['T_PW']),
                                           'T_wide': float(basis['T_wide'])},
                            'chars': results_this_lambda})

    # Skalen-Zusammenfassung pro Charakter
    print(f"\n{'=' * 90}")
    print("SKALENANALYSE PRO CHARAKTER")
    print(f"{'=' * 90}")
    for c_name in [c['name'] for c in chars]:
        print(f"\n{c_name}:")
        print(f"  {'lambda':>10s} {'log(L)':>8s} {'gamma1_in_PW':>14s} "
              f"{'Q_min':>10s} {'Q_mean':>10s} {'pos%':>6s}")
        for res_l in all_results:
            cr = next((r for r in res_l['chars'] if r['chi_name'] == c_name), None)
            if cr is None: continue
            vis = "ja" if cr['gamma1_in_PW'] else "nein"
            print(f"  {res_l['lambda']:10.3f} {res_l['log_lambda']:8.3f} {vis:>14s} "
                  f"{cr['Q_min']:10.4f} {cr['Q_mean']:10.4f} "
                  f"{cr['Q_positive_fraction']*100:5.1f}%")

    # Speichern
    out = RESULTS_DIR / "LAMBDA_SCALING_WEIL_BILINEAR_2026-04-18.json"
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({'config': {'N_galerkin': N_GALERKIN, 'rho_cutoff': RHO_CUTOFF,
                              'lambda_values': [float(x) for x in LAMBDA_VALUES]},
                   'results': all_results}, f, indent=2)
    print(f"\nRohdaten: {out}")

if __name__ == "__main__":
    main()
