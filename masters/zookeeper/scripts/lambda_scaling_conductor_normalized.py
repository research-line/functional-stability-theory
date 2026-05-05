"""
lambda_scaling_conductor_normalized.py
======================================

Conductor-bereinigte Version des lambda_scaling_weil_bilinear.py.

Kernerkenntnis aus LAMBDA_SCALING_WEIL_BILINEAR_2026-04-18.md:
  Der Conductor-Term -log(q/pi) * |h_n(0)|^2 im archimedischen Teil dominiert
  das Q_Weil-Signal bei nicht-trivialen chi (q > 1) und kleinem lambda. Das
  macht chi_21, chi_33 scheinbar "positiv" bei kleinem lambda, waehrend die
  echte Weil-Struktur dahinter versteckt bleibt.

Drei Versionen der Q_Weil-Form werden parallel berechnet:

  (A) Q_raw           = rho_diag - arch_without_conductor_diag - prim_diag
      (Matrix-Form ohne Conductor-Term)
  (B) Q_with_cond     = rho_diag - arch_full_diag - prim_diag
      (Matrix-Form mit Conductor-Term, wie frueher)
  (C) Q_residual      = Q_with_cond + log(q/pi) * |h_n(0)|^2
      (Conductor-bereinigt; sollte chi-weise vergleichbar sein)

Die Version (C) entfernt den Conductor-Artefakt und isoliert das "echte"
Weil-Positivitaets-Signal.

Zusaetzlich: Test mit erweiterter Lambda-Reihe (bis 16*sqrt(14)).

Autor: LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung 4)
Datum: 2026-04-18
"""

import json
import numpy as np
import time
from pathlib import Path
from scipy.special import digamma

LAMBDA_VALUES = [np.sqrt(14.0) * f for f in (1.0, 2.0, 4.0, 8.0, 16.0)]
N_GALERKIN = 30
N_L_TERMS = 400
RHO_CUTOFF = 50.0

RESULTS_DIR = Path(__file__).parent.parent / "_results"

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
    42.496097, 44.478293, 46.751013, 48.716374, 50.712570,
]

# =================================================================
# 1. Prolate-Basis
# =================================================================

def build_basis(lam, n_galerkin, rho_cutoff):
    L_val = np.log(lam)
    T_PW = 0.95 * L_val
    T_wide = max(rho_cutoff + 5.0, 3.0 * L_val)
    n_grid = min(3500, max(1200, int(40 * T_wide / max(L_val, 0.5))))
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
    return {'t': t, 'dt': dt, 'H': H, 'lambdas': ev[idx][:n_galerkin],
            'T_PW': T_PW, 'T_wide': T_wide, 'L': L_val, 'n_grid': n_grid, 'lambda_val': lam}

# =================================================================
# 2. Terme (getrennt nach Conductor)
# =================================================================

def arch_form_no_conductor(basis, parity):
    """Nur der digamma-Teil, ohne log(q/pi)-Term."""
    t = basis['t']
    arg = 0.25 if parity == +1 else 0.75
    v = 2.0 * np.real(digamma(arg + 1j*t/2.0))
    return (basis['H'].T * v[None, :]) @ basis['H'] * basis['dt']

def conductor_coupling(basis, q):
    """
    -log(q/pi) * |h_n(0)|^2 als Diagonalvektor (nur Diagonalbeitrag).
    h_n(0) = Wert der n-ten Basisfunktion am t=0 (Mellin-Variable).
    """
    t = basis['t']
    # Finde Index t=0
    idx_0 = np.argmin(np.abs(t))
    h_at_0 = basis['H'][idx_0, :]  # (N_gal,)
    conductor_offset = np.log(q / np.pi)
    return -conductor_offset * (h_at_0 ** 2)  # (N_gal,) Diagonalbeitrag

def prim_form(basis, chi, q, lam):
    t = basis['t']; dt = basis['dt']; H = basis['H']
    n = len(t)
    p_max = int(lam**2) + 1
    primes = [p for p in range(2, p_max+1)
              if all(p % q_ != 0 for q_ in range(2, int(np.sqrt(p))+1))]
    diff = t[:, None] - t[None, :]
    K = np.zeros((n, n), dtype=complex)
    for p in primes:
        if q > 1 and p % q == 0: continue
        lp = np.log(p); sp = np.sqrt(p)
        m_max = max(1, int(np.log(lam**2) / lp))
        for m in range(1, m_max+1):
            cpm = chi(p)**m
            if cpm == 0: continue
            w = -cpm * lp / (sp**m)
            K += w * np.cos(diff * m * lp)
    return H.T @ K @ H * dt**2

def interp_H(basis, gammas):
    t = basis['t']; H = basis['H']; dt = basis['dt']
    n_grid, n_gal = H.shape
    result = np.zeros((len(gammas), n_gal))
    for i, g in enumerate(gammas):
        if g < t[0] or g > t[-1]: continue
        idx = np.searchsorted(t, g)
        if idx <= 0: result[i, :] = H[0, :]
        elif idx >= n_grid: result[i, :] = H[-1, :]
        else:
            alpha = (g - t[idx-1]) / dt
            result[i, :] = (1-alpha)*H[idx-1,:] + alpha*H[idx,:]
    return result

def rho_diag(basis, zeros, cutoff):
    valid = [g for g in zeros if abs(g) <= cutoff]
    if not valid:
        return np.zeros(basis['H'].shape[1])
    gammas = np.array(valid + [-g for g in valid])
    Hg = interp_H(basis, gammas)
    return np.sum(Hg * Hg, axis=0)

# =================================================================
# 3. Q_Weil-Drei-Version-Test
# =================================================================

def test_Q_three(basis, chi_info):
    arch_mat = arch_form_no_conductor(basis, chi_info['parity'])
    prim_mat = prim_form(basis, chi_info['chi'], chi_info['q'], basis['lambda_val'])
    r_diag = rho_diag(basis, chi_info['zeros'], RHO_CUTOFF)
    cond_diag = conductor_coupling(basis, chi_info['q'])

    arch_d = np.real(np.diag(arch_mat))
    prim_d = np.real(np.diag(prim_mat))

    # (A) ohne Conductor
    Q_raw = r_diag - arch_d - prim_d

    # (B) mit Conductor (= Version aus LAMBDA_SCALING_WEIL_BILINEAR)
    Q_with_cond = Q_raw + cond_diag  # cond_diag = -log(q/pi)*|h_n(0)|^2, d.h. Q_with = Q_raw - log(q/pi)*|h(0)|^2

    # (C) Conductor-bereinigt (Q_with_cond nach Rueckkorrektur)
    # Diese sollte NICHT mehr den Conductor-Bias haben
    # Aber das ist algebraisch = Q_raw (cond herausgerechnet)
    # Genaugenommen: falls Q_with_cond = Q_raw - log(q/pi)*|h(0)|^2,
    # dann Q_residual = Q_with_cond + log(q/pi)*|h(0)|^2 = Q_raw.
    # Also ist (C) identisch mit (A). Das ist die "correct" Conductor-Bereinigung.
    # Der Unterschied: (A) ignoriert Conductor ganz, (C) entfernt ihn nach Einbau.
    # Fuer unseren Test: nimm (A) = (C).
    Q_residual = Q_raw.copy()

    return {
        'chi_name': chi_info['name'],
        'q': chi_info['q'],
        'gamma1': chi_info['gamma1'],
        'Q_raw_min': float(Q_raw.min()),
        'Q_raw_mean': float(Q_raw.mean()),
        'Q_raw_pos_frac': float(np.mean(Q_raw >= 0)),
        'Q_with_cond_min': float(Q_with_cond.min()),
        'Q_with_cond_mean': float(Q_with_cond.mean()),
        'Q_with_cond_pos_frac': float(np.mean(Q_with_cond >= 0)),
        'conductor_offset': float(np.log(chi_info['q'] / np.pi)),
        'mean_conductor_contrib': float(np.mean(cond_diag)),
    }

# =================================================================
# 4. Charaktere
# =================================================================

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
    c0 = mkchi('chi_0', 1, +1, {}, RIEMANN_ZEROS)
    def chi_t(n): return 1.0 + 0j
    c0['chi'] = chi_t
    chars.append(c0)
    chars.append(mkchi('chi_4', 4, -1, {1: 1.0+0j, 3: -1.0+0j}, CHI4_ZEROS))
    chars.append(mkchi('chi_5', 5, +1,
                       {1: 1.0+0j, 2: -1.0+0j, 3: -1.0+0j, 4: 1.0+0j},
                       za['chi_5']))
    chi_21_arr = [0, 1, -1, 0, 1, 1, 0, 0, -1, 0, -1, -1, 0, 1, 0, 0, 1, 1, 0, -1, 1]
    chars.append(mkchi('chi_21', 21, +1,
                       {n: float(chi_21_arr[n])+0j for n in range(21) if chi_21_arr[n] != 0},
                       za['chi_21']))
    arr33 = [0,1,1,0,1,-1,0,-1,1,0,-1,0,0,-1,-1,0,1,1,0,-1,-1,0,0,-1,0,1,-1,0,-1,1,0,1,1]
    chars.append(mkchi('chi_33', 33, +1,
                       {n: float(arr33[n])+0j for n in range(33) if arr33[n] != 0},
                       za['chi_33']))
    return chars

# =================================================================
# 5. Main
# =================================================================

def main():
    print("=" * 100)
    print("lambda_scaling_conductor_normalized.py")
    print("Q_Weil OHNE Conductor-Term vs MIT Conductor-Term vs BEREINIGT")
    print(f"N_GALERKIN={N_GALERKIN}, rho_cutoff={RHO_CUTOFF}")
    print(f"lambda_values: {[f'{l:.2f}' for l in LAMBDA_VALUES]}")
    print("=" * 100)

    chars = load_chars()
    all_res = []
    for lam in LAMBDA_VALUES:
        t0 = time.time()
        basis = build_basis(lam, N_GALERKIN, RHO_CUTOFF)
        print(f"\nlambda={lam:.3f}, log lambda={basis['L']:.3f}, N_grid={basis['n_grid']}, "
              f"T_wide={basis['T_wide']:.1f}, build={time.time()-t0:.1f}s")
        res_l = {'lambda': float(lam), 'log_lambda': float(basis['L']), 'chars': []}
        print(f"  {'chi':8s} {'gamma1':>7s} {'Q_raw_mean':>11s} {'Q_cond_mean':>12s} "
              f"{'Q_raw_min':>10s} {'Q_raw_pos%':>10s} {'cond_mean':>10s}")
        for c in chars:
            r = test_Q_three(basis, c)
            r['lambda'] = float(lam)
            r['log_lambda'] = float(basis['L'])
            r['gamma1_in_PW'] = bool(c['gamma1'] < basis['T_PW'])
            res_l['chars'].append(r)
            print(f"  {r['chi_name']:8s} {r['gamma1']:7.3f} {r['Q_raw_mean']:11.4f} "
                  f"{r['Q_with_cond_mean']:12.4f} {r['Q_raw_min']:10.4f} "
                  f"{r['Q_raw_pos_frac']*100:9.1f}% {r['mean_conductor_contrib']:10.4f}")
        all_res.append(res_l)

    # Skalenanalyse Q_raw (OHNE Conductor)
    print(f"\n{'='*100}")
    print("Q_RAW (OHNE Conductor-Term) pro Charakter:")
    print(f"{'='*100}")
    for c_name in [c['name'] for c in chars]:
        print(f"\n{c_name}:")
        print(f"  {'lambda':>10s} {'log L':>7s} {'Q_raw_min':>10s} {'Q_raw_mean':>11s} {'pos%':>6s}")
        for rl in all_res:
            cr = next((r for r in rl['chars'] if r['chi_name'] == c_name), None)
            if cr is None: continue
            print(f"  {rl['lambda']:10.3f} {rl['log_lambda']:7.3f} {cr['Q_raw_min']:10.4f} "
                  f"{cr['Q_raw_mean']:11.4f} {cr['Q_raw_pos_frac']*100:5.1f}%")

    # H2-Check bei groesstem lambda
    if len(all_res) >= 2:
        print(f"\n{'='*100}")
        print(f"H2-CHECK bei lambda={all_res[-1]['lambda']:.2f}:")
        print(f"{'='*100}")
        last = all_res[-1]['chars']
        evens = [r for r in last if r['chi_name'] in ('chi_5', 'chi_21', 'chi_33')]
        if len(evens) == 3:
            print(f"  {'chi':8s} {'gamma1':>7s} {'Q_raw_mean':>11s} {'Q_raw_min':>10s}")
            for r in sorted(evens, key=lambda x: -x['gamma1']):
                print(f"  {r['chi_name']:8s} {r['gamma1']:7.3f} {r['Q_raw_mean']:11.4f} "
                      f"{r['Q_raw_min']:10.4f}")
            gs = np.array([r['gamma1'] for r in evens])
            # Nicht-log, da Q_raw kann negativ sein. Nimm -Q_raw_mean (ist durchweg negativ)
            vals = np.array([-r['Q_raw_mean'] for r in evens])
            if np.all(vals > 0):
                slope, _ = np.polyfit(np.log(gs), np.log(vals), 1)
                print(f"\n  Slope log(-Q_raw_mean) vs log(gamma1): {slope:+.3f}   (H2 erwartet +1, H1 erwartet +3)")
                print(f"  Interpretation: kleineres gamma1 -> weniger negativ -> positiver slope erwartet.")

    out = RESULTS_DIR / "LAMBDA_SCALING_CONDUCTOR_NORMALIZED_2026-04-18.json"
    json.dump({'config': {'N_galerkin': N_GALERKIN, 'rho_cutoff': RHO_CUTOFF,
                          'lambda_values': [float(x) for x in LAMBDA_VALUES]},
               'results': all_res}, open(out, 'w'), indent=2)
    print(f"\nRohdaten: {out}")

if __name__ == "__main__":
    main()
