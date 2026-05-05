"""
lambda_scaling_full_atlas.py
============================

VOLLVERSION des lambda-Skalierungs-Tests mit:

  1. Conductor-Erweiterung im arch-Term (log(q/pi) * |h~(0)|^2 als Rank-1).
  2. Parity-korrektes digamma-Argument (1/4 fuer even, 3/4 fuer odd).
  3. Alle 10 Atlas-Charaktere (chi_5, chi_8, chi_12, chi_13, chi_17,
     chi_21, chi_24, chi_29, chi_33, chi_60) mit Kronecker-Werten
     via sympy verifiziert.
  4. Erweiterte Lambda-Skala (bis 128*sqrt(14)).
  5. Zwei Positivitaets-Tests parallel:
     - Matrix-Test:  min_Eigenvalue(M) mit M = G_rho - G_arch - G_prim
     - Diagonal-Test: min_n Q_n (wie in vorherigen Varianten)

Weil-Form (Iwaniec-Kowalski Thm 5.12 umformuliert):

  Fuer primitiven Charakter chi mod q, parity epsilon = chi(-1):
    Q_Weil(f, f) = sum_rho |f~(gamma_rho)|^2
                   - log(q/pi) * |f~(0)|^2
                   - int |f(t)|^2 * 2 Re[psi(arg_par + it/2)] dt
                   - prim-Summe mit chi(p)^m
  wobei arg_par = 1/4 (even) oder 3/4 (odd).

Rechenaufwand-Schaetzung:
  - N_grid = 2000-2500, eigh ~ 10-15s pro lambda
  - 5 lambdas * 10 chars * 3s pro char = 150s = 2.5min
  - Gesamt: ca. 3-5 min lokal

Autor: LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung 5)
Datum: 2026-04-18
"""

import json
import numpy as np
import time
from pathlib import Path
from scipy.special import digamma

LAMBDA_VALUES = [np.sqrt(14.0) * f for f in (1.0, 4.0, 16.0, 64.0, 128.0)]
N_GALERKIN = 35
N_L_TERMS = 400
RHO_CUTOFF = 50.0

RESULTS_DIR = Path(__file__).parent.parent / "_results"

# 50 Riemann-Nullstellen
RIEMANN_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005150, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
    79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
    92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
    103.725538, 105.446623, 107.168612, 111.029536, 111.874659,
    114.320221, 116.226680, 118.790783, 121.370125, 122.946829,
]
CHI4_ZEROS = [
    6.020948, 10.243766, 12.988096, 16.343297, 18.291996,
    21.428273, 23.265376, 26.068044, 28.106108, 30.296575,
    31.765775, 34.479648, 35.872660, 37.967690, 40.340230,
    42.496097, 44.478293, 46.751013, 48.716374, 50.712570,
]

# Kronecker-Werte via sympy verifiziert (D > 0 ist Discriminante, even parity)
KRONECKER_TABLES = {
    5:  [0, 1, -1, -1, 1],
    8:  [0, 1, 0, -1, 0, -1, 0, 1],
    12: [0, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 1],
    13: [0, 1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, 1],
    17: [0, 1, 1, -1, 1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1],
    21: [0, 1, -1, 0, 1, 1, 0, 0, -1, 0, -1, -1, 0, -1, 0, 0, 1, 1, 0, -1, 1],
    24: [0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0, -1, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1],
    29: [0, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, -1, -1, 1],
    33: [0, 1, 1, 0, 1, -1, 0, -1, 1, 0, -1, 0, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 0, -1, 0, 1, -1, 0, -1, 1, 0, 1, 1],
    60: [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 1, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
}

# =================================================================

def build_basis(lam, n_galerkin, rho_cutoff):
    L_val = np.log(lam)
    T_PW = 0.95 * L_val
    T_wide = max(rho_cutoff + 5.0, 3.0 * L_val)
    n_grid = min(2500, max(1200, int(25 * T_wide / max(L_val, 0.5))))
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

def arch_matrix_full(basis, parity, q):
    """
    Weil-Archimedisch mit Conductor-Term:
        arch_matrix = log(q/pi) * H^T H + int h_m(t) * 2*Re[psi(arg+it/2)] h_n(t) dt
    wobei arg = 1/4 (even) oder 3/4 (odd), und H^T H hat Rank = 1 (Projektor auf konstante).
    """
    t = basis['t']; H = basis['H']; dt = basis['dt']
    arg = 0.25 if parity == +1 else 0.75
    v = 2.0 * np.real(digamma(arg + 1j*t/2.0))
    # Integral-Teil
    arch_integral = (H.T * v[None, :]) @ H * dt
    # Conductor-Teil: log(q/pi) * h_m~(0) * h_n~(0) mit h~(0) = sum h(t_k) * dt
    h_at_0_tilde = np.sum(H, axis=0) * dt  # (N_gal,)
    conductor_matrix = np.log(q / np.pi) * np.outer(h_at_0_tilde, h_at_0_tilde)
    return arch_integral + conductor_matrix

def prim_matrix(basis, chi, q, lam, p_cap=200):
    """
    Prim-Matrix via cos+sin Fourier-Komponenten (KORREKT fuer beliebige Basen).

    Aus int int h_m(t) cos(u*(t-t')) h_n(t') dt dt'
       = [int h_m(t) cos(ut) dt][int h_n(t') cos(ut') dt']
       + [int h_m(t) sin(ut) dt][int h_n(t') sin(ut') dt']

    Fuer Prolate-Basis mit geraden UND ungeraden Moden sind beide Integrale
    noetig (vorherige Version mit nur cos hat sin-Beitrag fuer ungerade Moden
    ignoriert).

    Cutoff: p <= min(lambda^2, p_cap) um Laufzeit zu begrenzen.
    """
    t = basis['t']; dt = basis['dt']; H = basis['H']
    N_gal = H.shape[1]
    p_max = min(int(lam**2) + 1, p_cap)
    primes = [p for p in range(2, p_max+1)
              if all(p % q_ != 0 for q_ in range(2, int(np.sqrt(p))+1))]
    K_proj = np.zeros((N_gal, N_gal), dtype=complex)
    for p in primes:
        if q > 1 and p % q == 0: continue
        lp = np.log(p); sp = np.sqrt(p)
        m_max = max(1, int(np.log(lam**2) / lp))
        for m in range(1, m_max+1):
            cpm = chi(p)**m
            if cpm == 0: continue
            w = cpm * lp / (sp**m)
            u = m * lp
            cos_vec = np.cos(u * t)
            sin_vec = np.sin(u * t)
            F_cos = (H * cos_vec[:, None]).sum(axis=0) * dt
            F_sin = (H * sin_vec[:, None]).sum(axis=0) * dt
            K_proj += w * (np.outer(F_cos, F_cos) + np.outer(F_sin, F_sin))
    return K_proj

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

def rho_matrix(basis, zeros, cutoff):
    valid = [g for g in zeros if abs(g) <= cutoff]
    if not valid:
        N = basis['H'].shape[1]
        return np.zeros((N, N))
    gammas = np.array(valid + [-g for g in valid])
    Hg = interp_H(basis, gammas)  # (2K, N_gal)
    return Hg.T @ Hg  # (N_gal, N_gal) rank-2K projektor

def test_Q(basis, chi_info):
    """
    Q_Weil = G_rho - G_arch - G_prim  (als Matrix)
    Zwei Tests:
      (A) Diagonal-Positivitaet: min_n Q[n,n] >= 0
      (B) Matrix-Positivitaet:  min EW(Q) >= 0

    (B) ist strenger als (A); (A) ist interpretiert als Weil-Positivitaet
    fuer jede einzelne Prolate-Mode als Testfunktion.
    """
    arch = arch_matrix_full(basis, chi_info['parity'], chi_info['q'])
    prim = prim_matrix(basis, chi_info['chi'], chi_info['q'], basis['lambda_val'])
    rho  = rho_matrix(basis, chi_info['zeros'], RHO_CUTOFF)
    # Q = rho - arch - prim
    Q_mat = rho - np.real(arch) - np.real(prim)
    Q_mat_herm = 0.5 * (Q_mat + Q_mat.T)  # Symmetrisierung

    Q_diag = np.real(np.diag(Q_mat_herm))
    ew = np.linalg.eigvalsh(Q_mat_herm)

    n_rho = sum(1 for g in chi_info['zeros'] if abs(g) <= RHO_CUTOFF)

    return {
        'chi_name': chi_info['name'],
        'q': chi_info['q'],
        'parity': chi_info['parity'],
        'gamma1': chi_info['gamma1'],
        'n_rho_used': n_rho,
        # Diagonal-Test
        'Q_diag_min': float(Q_diag.min()),
        'Q_diag_mean': float(Q_diag.mean()),
        'Q_diag_pos_frac': float(np.mean(Q_diag >= 0)),
        # Matrix-EW-Test (strenger)
        'Q_matrix_ew_min': float(ew.min()),
        'Q_matrix_ew_max': float(ew.max()),
        'n_negative_ew': int(np.sum(ew < -1e-8)),
    }

def mkchi_from_kronecker(name, D, zeros):
    """D > 0 = Discriminante, even parity. Kronecker-Tabelle aus KRONECKER_TABLES."""
    vals = KRONECKER_TABLES[D]
    q = D
    arr = np.array([complex(v) for v in vals])
    def chi(n): return arr[n % q]
    return {'name': name, 'q': q, 'parity': +1, 'chi': chi, 'zeros': zeros,
            'gamma1': zeros[0] if zeros else 0}

def load_chars():
    zf = Path(__file__).parent.parent.parent / "dirichlet_atlas" / "_results" / "zeros_all_chars.json"
    za = json.load(open(zf))

    chars = []
    # chi_0 (Riemann)
    def chi_triv(n): return 1.0 + 0j
    chars.append({'name': 'chi_0', 'q': 1, 'parity': +1, 'chi': chi_triv,
                  'zeros': RIEMANN_ZEROS, 'gamma1': RIEMANN_ZEROS[0]})
    # chi_4 (odd)
    arr4 = np.array([0.0+0j, 1.0+0j, 0.0+0j, -1.0+0j])
    def chi4(n): return arr4[n % 4]
    chars.append({'name': 'chi_4', 'q': 4, 'parity': -1, 'chi': chi4,
                  'zeros': CHI4_ZEROS, 'gamma1': CHI4_ZEROS[0]})
    # Atlas-Charaktere
    atlas_names = ['chi_5', 'chi_8', 'chi_12', 'chi_13', 'chi_17',
                    'chi_21', 'chi_24', 'chi_29', 'chi_33', 'chi_60']
    for name in atlas_names:
        D = int(name.split('_')[1])
        chars.append(mkchi_from_kronecker(name, D, za[name]))
    return chars

# =================================================================

def main():
    print("="*110)
    print("lambda_scaling_full_atlas.py")
    print(f"Voller Q_Weil mit Conductor, {len(KRONECKER_TABLES)+2} Charakteren, {len(LAMBDA_VALUES)} lambda-Werten.")
    print(f"N_GALERKIN={N_GALERKIN}, rho_cutoff={RHO_CUTOFF}")
    print("="*110)

    chars = load_chars()
    all_res = []
    for lam in LAMBDA_VALUES:
        t0 = time.time()
        basis = build_basis(lam, N_GALERKIN, RHO_CUTOFF)
        print(f"\nlambda={lam:.3f} (log={basis['L']:.3f}), N_grid={basis['n_grid']}, "
              f"T_wide={basis['T_wide']:.1f}, T_PW={basis['T_PW']:.3f}, "
              f"build={time.time()-t0:.1f}s")
        print(f"  {'chi':8s} {'q':>3s} {'par':>4s} {'gamma1':>7s} {'n_rho':>6s} "
              f"{'Q_diag_min':>11s} {'Q_ew_min':>10s} {'pos_diag%':>10s} {'#neg_EW':>8s}")
        res_l = {'lambda': float(lam), 'log_lambda': float(basis['L']),
                 'T_PW': float(basis['T_PW']), 'chars': []}
        for c in chars:
            r = test_Q(basis, c)
            r['lambda'] = float(lam)
            r['log_lambda'] = float(basis['L'])
            r['gamma1_in_PW'] = bool(c['gamma1'] < basis['T_PW'])
            res_l['chars'].append(r)
            mark = " PW" if r['gamma1_in_PW'] else "   "
            print(f"  {r['chi_name']:8s} {r['q']:3d} {r['parity']:+4d} {r['gamma1']:7.3f} "
                  f"{r['n_rho_used']:6d} {r['Q_diag_min']:11.3f} {r['Q_matrix_ew_min']:10.3f} "
                  f"{r['Q_diag_pos_frac']*100:9.1f}% {r['n_negative_ew']:8d}{mark}")
        all_res.append(res_l)

    # Zusammenfassung pro Charakter
    print(f"\n{'='*110}")
    print("Q_diag_min und Q_ew_min pro Charakter ueber lambda:")
    print(f"{'='*110}")
    for c_name in [c['name'] for c in chars]:
        print(f"\n{c_name}:")
        print(f"  {'lambda':>10s} {'log_L':>7s} {'Q_diag_min':>12s} {'Q_ew_min':>10s} "
              f"{'Q_diag_mean':>12s} {'pos_diag%':>10s} {'#neg_EW':>8s}")
        for rl in all_res:
            cr = next((r for r in rl['chars'] if r['chi_name'] == c_name), None)
            if cr is None: continue
            print(f"  {rl['lambda']:10.2f} {rl['log_lambda']:7.2f} "
                  f"{cr['Q_diag_min']:12.4f} {cr['Q_matrix_ew_min']:10.4f} "
                  f"{cr['Q_diag_mean']:12.4f} {cr['Q_diag_pos_frac']*100:9.1f}% "
                  f"{cr['n_negative_ew']:8d}")

    # H2-Analyse: bei groesstem lambda, sortiert nach gamma1
    print(f"\n{'='*110}")
    print(f"H2-STRUKTUR bei lambda={all_res[-1]['lambda']:.2f}:")
    print(f"{'='*110}")
    last = sorted(all_res[-1]['chars'], key=lambda r: r['gamma1'])
    print(f"  {'chi':8s} {'gamma1':>7s} {'1/gamma1':>9s} {'Q_diag_min':>11s} {'Q_ew_min':>10s}")
    for r in last:
        print(f"  {r['chi_name']:8s} {r['gamma1']:7.3f} {1/r['gamma1']:9.4f} "
              f"{r['Q_diag_min']:11.4f} {r['Q_matrix_ew_min']:10.4f}")

    # Speichern
    out = RESULTS_DIR / "LAMBDA_SCALING_FULL_ATLAS_2026-04-18.json"
    json.dump({'config': {'N_galerkin': N_GALERKIN, 'rho_cutoff': RHO_CUTOFF,
                          'lambda_values': [float(x) for x in LAMBDA_VALUES]},
               'results': all_res}, open(out, 'w'), indent=2)
    print(f"\nRohdaten: {out}")

if __name__ == "__main__":
    main()
