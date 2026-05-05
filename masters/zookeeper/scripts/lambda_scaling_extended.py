"""
lambda_scaling_extended.py
==========================

Erweitert den Test auf groessere lambda-Werte: 32*sqrt(14), 64*sqrt(14),
um zu pruefen ob:
  (a) Q_raw_min fuer chi_0 monoton gegen 0 ... und dann positiv bleibt.
  (b) H2-Slope bei Even-Familie sich bei grossem lambda stabilisiert.

Basis: lambda_scaling_conductor_normalized.py (ohne Conductor-Term).

Rechenaufwand: N_grid waechst mit T_wide = max(rho_cutoff+5, 3*L). Bei
lambda=64*sqrt(14), log lambda ~ 5.48, T_wide ~ 55 (wie vorher), dt ~ 0.046
bei N_grid=2400, eigh(2400)^3 ~ 15s. Akzeptabel lokal.

Fuer lambda=256*sqrt(14) (log=6.87): T_wide braucht hoeher, N_grid ~ 3500,
eigh ~ 45s. Auch lokal ok.

Autor: LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung 5)
Datum: 2026-04-18
"""

import json
import numpy as np
import time
from pathlib import Path
from scipy.special import digamma

LAMBDA_VALUES = [np.sqrt(14.0) * f for f in (1.0, 4.0, 16.0, 32.0)]
N_GALERKIN = 30
N_L_TERMS = 400
RHO_CUTOFF = 40.0

RESULTS_DIR = Path(__file__).parent.parent / "_results"

# 50 Riemann-Nullstellen (Odlyzko-Tabelle)
RIEMANN_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005150, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
    79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
    92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
    103.725538, 105.446623, 107.168612, 111.029536, 111.874659,
    114.320221, 116.226680, 118.790783, 121.370125, 122.946829,
    124.256819, 127.516684, 129.578704, 131.087689, 133.497737,
    134.756510, 138.116042, 139.736209, 141.123707, 143.111846,
]
CHI4_ZEROS = [
    6.020948, 10.243766, 12.988096, 16.343297, 18.291996,
    21.428273, 23.265376, 26.068044, 28.106108, 30.296575,
    31.765775, 34.479648, 35.872660, 37.967690, 40.340230,
    42.496097, 44.478293, 46.751013, 48.716374, 50.712570,
    52.430985, 54.518293, 56.336024, 58.055497, 59.853076,
]

def build_basis(lam, n_galerkin, rho_cutoff):
    L_val = np.log(lam)
    T_PW = 0.95 * L_val
    T_wide = max(rho_cutoff + 5.0, 3.0 * L_val)
    # n_grid adaptive: kleiner halten fuer Performance
    n_grid = min(2500, max(1200, int(30 * T_wide / max(L_val, 0.5))))
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

def arch_form(basis, parity):
    t = basis['t']
    arg = 0.25 if parity == +1 else 0.75
    v = 2.0 * np.real(digamma(arg + 1j*t/2.0))
    return (basis['H'].T * v[None, :]) @ basis['H'] * basis['dt']

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

def test_Q(basis, chi_info):
    arch_mat = arch_form(basis, chi_info['parity'])
    prim_mat = prim_form(basis, chi_info['chi'], chi_info['q'], basis['lambda_val'])
    r_diag = rho_diag(basis, chi_info['zeros'], RHO_CUTOFF)
    arch_d = np.real(np.diag(arch_mat))
    prim_d = np.real(np.diag(prim_mat))
    Q = r_diag - arch_d - prim_d
    n_rho_used = sum(1 for g in chi_info['zeros'] if abs(g) <= RHO_CUTOFF)
    return {
        'chi_name': chi_info['name'],
        'q': chi_info['q'],
        'gamma1': chi_info['gamma1'],
        'n_rho_used': n_rho_used,
        'Q_min': float(Q.min()),
        'Q_max': float(Q.max()),
        'Q_mean': float(Q.mean()),
        'Q_pos_frac': float(np.mean(Q >= 0)),
        'Q_first_5': Q[:5].tolist(),
    }

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

def main():
    print("="*100)
    print("lambda_scaling_extended.py (bis 64*sqrt(14))")
    print(f"N_GALERKIN={N_GALERKIN}, rho_cutoff={RHO_CUTOFF}")
    print(f"lambdas: {[f'{l:.2f}' for l in LAMBDA_VALUES]}")
    print("="*100)
    chars = load_chars()
    all_res = []
    for lam in LAMBDA_VALUES:
        t0 = time.time()
        basis = build_basis(lam, N_GALERKIN, RHO_CUTOFF)
        print(f"\nlambda={lam:.3f} (log={basis['L']:.3f}), N_grid={basis['n_grid']}, "
              f"T_wide={basis['T_wide']:.1f}, build={time.time()-t0:.1f}s")
        print(f"  {'chi':8s} {'gamma1':>7s} {'n_rho':>6s} {'Q_min':>10s} {'Q_mean':>10s} {'pos%':>6s}")
        res_l = {'lambda': float(lam), 'log_lambda': float(basis['L']), 'T_PW': float(basis['T_PW']),
                 'chars': []}
        for c in chars:
            r = test_Q(basis, c)
            r['lambda'] = float(lam); r['log_lambda'] = float(basis['L'])
            r['gamma1_in_PW'] = bool(c['gamma1'] < basis['T_PW'])
            res_l['chars'].append(r)
            mark = " <--gamma1 in PW" if r['gamma1_in_PW'] else ""
            print(f"  {r['chi_name']:8s} {r['gamma1']:7.3f} {r['n_rho_used']:6d} "
                  f"{r['Q_min']:10.4f} {r['Q_mean']:10.4f} {r['Q_pos_frac']*100:5.1f}%{mark}")
        all_res.append(res_l)

    # Q_min-Monotonie pro Charakter
    print(f"\n{'='*100}")
    print("Q_min vs lambda:")
    print(f"{'='*100}")
    for c_name in [c['name'] for c in chars]:
        print(f"\n{c_name}:")
        Q_mins = []
        Q_means = []
        logLs = []
        for rl in all_res:
            cr = next((r for r in rl['chars'] if r['chi_name'] == c_name), None)
            if cr is None: continue
            print(f"  lambda={rl['lambda']:8.2f}  log(L)={rl['log_lambda']:5.2f}  "
                  f"Q_min={cr['Q_min']:+10.4f}  Q_mean={cr['Q_mean']:+10.4f}  "
                  f"pos={cr['Q_pos_frac']*100:5.1f}%  n_rho={cr['n_rho_used']}")
            Q_mins.append(cr['Q_min'])
            Q_means.append(cr['Q_mean'])
            logLs.append(rl['log_lambda'])
        # Slope von Q_min gegen log(lambda)
        Q_mins = np.array(Q_mins); logLs = np.array(logLs)
        if len(Q_mins) >= 3:
            slope, intercept = np.polyfit(logLs, Q_mins, 1)
            # Extrapolation: bei welchem log(L) wird Q_min = 0?
            if slope > 1e-6:
                log_lambda_zero = -intercept / slope
                lambda_zero = np.exp(log_lambda_zero)
                print(f"  -> Q_min-Slope vs log(lambda) = {slope:+.4f}; extrapoliert Q_min=0 bei lambda={lambda_zero:.1f}")
            else:
                print(f"  -> Q_min-Slope vs log(lambda) = {slope:+.4f} (kein monotones Wachstum)")

    out = RESULTS_DIR / "LAMBDA_SCALING_EXTENDED_2026-04-18.json"
    json.dump({'config': {'N_galerkin': N_GALERKIN, 'rho_cutoff': RHO_CUTOFF,
                          'lambda_values': [float(x) for x in LAMBDA_VALUES]},
               'results': all_res}, open(out, 'w'), indent=2)
    print(f"\nRohdaten: {out}")

if __name__ == "__main__":
    main()
