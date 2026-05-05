"""
c2bd_full_n25_run.py

Vollstaendiger Derivative-Freezing-Lauf bei N=25, NGRID=30 fuer alle 5 Lambda-Werte.
Ziel: konsistente K_lambda-Werte auf gleichem Trunkierungsniveau.
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import c2bd_derivative_matching as mod
mod.NGRID = 30

DPS = int(os.environ.get("DPS", 25))

from c2bd_derivative_matching import compute_full_data


def run_collect(lam_base, N, dlam=0.01):
    L_base = 2 * np.log(lam_base)
    dL = 2 * dlam / lam_base

    print(f"\n  lam={lam_base:.1f}, N={N}, NGRID={mod.NGRID}, dlam={dlam}, dL={dL:.6e}")

    t0 = time.time()
    d_m = compute_full_data(lam_base - dlam, N)
    d_0 = compute_full_data(lam_base, N)
    d_p = compute_full_data(lam_base + dlam, N)
    elapsed = time.time() - t0
    print(f"    3 Punkte in {elapsed:.0f}s", flush=True)

    dL_T2bd = (d_p['T2_boundary'] - d_m['T2_boundary']) / (2 * dL)
    dL_Rbulk = (d_p['R_bulk'] - d_m['R_bulk']) / (2 * dL)
    dL_mua = (d_p['mu_alpha'] - d_m['mu_alpha']) / (2 * dL)
    dL_u0 = (d_p['u0'] - d_m['u0']) / (2 * dL)
    dL_Ct = (d_p['Ctilde'] - d_m['Ctilde']) / (2 * dL)

    u0_base = d_0['u0']
    mH_m = sum(b**2 / (h - u0_base) for b, h in zip(d_m['betas'], d_m['h_vals'])
               if abs(h - u0_base) > 1e-30)
    mH_p = sum(b**2 / (h - u0_base) for b, h in zip(d_p['betas'], d_p['h_vals'])
               if abs(h - u0_base) > 1e-30)
    dL_mH = (mH_p - mH_m) / (2 * dL)
    mHp = d_0['mHp']
    Ct = d_0['Ctilde']
    rhs = dL_Ct / Ct**2

    cancel = abs(dL_mua) / (abs(dL_T2bd) + abs(dL_Rbulk)) if (abs(dL_T2bd) + abs(dL_Rbulk)) > 0 else float('nan')
    ratio = dL_T2bd / dL_Rbulk if abs(dL_Rbulk) > 1e-40 else float('nan')
    K_lam = -dL_T2bd / (lam_base**3 * dL_mH) if abs(lam_base**3 * dL_mH) > 1e-40 else float('nan')

    return {
        'lam': lam_base, 'L': L_base,
        'dL_T2bd': dL_T2bd, 'dL_Rbulk': dL_Rbulk, 'dL_mua': dL_mua,
        'dL_mH': dL_mH, 'mHp_u0p': mHp * dL_u0, 'Ct_p_Ct2': rhs,
        'K_lam': K_lam, 'cancel_pct': (1 - cancel) * 100, 'ratio': ratio,
        'T2_bd': d_0['T2_boundary'], 'R_bulk': d_0['R_bulk'],
        'mu_alpha': d_0['mu_alpha'], 'alpha': d_0['alpha'],
    }


def main():
    print(f"C2bd: Full N=25, NGRID=30 Run")
    print(f"DPS={DPS}")

    results = []
    for lam in [3.0, 4.0, 5.0, 7.0, 9.0]:
        r = run_collect(lam, N=25, dlam=0.01)
        results.append(r)

    print(f"\n{'='*90}")
    print(f"  ZUSAMMENFASSUNG (N=25, NGRID=30)")
    print(f"{'='*90}")

    print(f"\n  (1) DERIVATIVE FREEZING")
    print(f"  {'lam':>5} {'Cancel%':>10} {'d/dL(T2_bd)':>14} {'d/dL(R_bulk)':>14} {'d/dL(mu/a)':>14} {'Ratio':>12}")
    for r in results:
        print(f"  {r['lam']:5.1f} {r['cancel_pct']:10.1f}% {r['dL_T2bd']:+14.6e} {r['dL_Rbulk']:+14.6e} "
              f"{r['dL_mua']:+14.6e} {r['ratio']:+12.6f}")

    print(f"\n  (2) K_LAMBDA = -d(T2_bd)/dL / (lam^3 * dL_mH)")
    print(f"  {'lam':>5} {'K_lam':>14} {'lam^3':>8} {'dL_mH':>14} {'sign(dT2bd)':>12}")
    for r in results:
        sgn = "+" if r['dL_T2bd'] > 0 else "-"
        print(f"  {r['lam']:5.1f} {r['K_lam']:+14.6e} {r['lam']**3:8.0f} {r['dL_mH']:+14.6e} {sgn:>12}")

    k_arr = np.array([r['K_lam'] for r in results])
    good = [(r['K_lam'], r['lam']) for r in results if r['cancel_pct'] >= 99.5]
    if good:
        k_good = np.array([g[0] for g in good])
        lams_good = [g[1] for g in good]
        print(f"\n  Punkte mit Cancel >= 99.5%: lam = {lams_good}")
        print(f"  K_mean = {np.mean(k_good):+.6e} +/- {np.std(k_good):.2e}")
        print(f"  CV = {np.std(k_good)/abs(np.mean(k_good)):.4f} ({np.std(k_good)/abs(np.mean(k_good))*100:.1f}%)")

    print(f"\n  (3) SECULAR vs BOUNDARY-ABLEITUNGEN")
    print(f"  {'lam':>5} {'dL_mH':>14} {'mHp*u0p':>14} {'Ct_p/Ct2':>14} {'dL_mua/Ct_p':>14}")
    for r in results:
        cr = r['dL_mua'] / r['Ct_p_Ct2'] if abs(r['Ct_p_Ct2']) > 1e-40 else float('nan')
        print(f"  {r['lam']:5.1f} {r['dL_mH']:+14.6e} {r['mHp_u0p']:+14.6e} {r['Ct_p_Ct2']:+14.6e} {cr:+14.6e}")

    # Power-law fit of |K| vs lambda
    print(f"\n  (4) SKALIERUNGS-ANALYSE")
    K_abs = np.array([abs(r['K_lam']) for r in results])
    K_sign = np.array([np.sign(r['K_lam']) for r in results])
    lams = np.array([r['lam'] for r in results])
    print(f"  {'lam':>5} {'|K|':>14} {'sign(K)':>8} {'|dL_T2bd|':>14} {'|dL_mH|':>14}")
    for r in results:
        print(f"  {r['lam']:5.1f} {abs(r['K_lam']):14.6e} {'+' if r['K_lam']>0 else '-':>8} "
              f"{abs(r['dL_T2bd']):14.6e} {abs(r['dL_mH']):14.6e}")

    # Power-law exponent of |dL_T2bd|
    print(f"\n  Power-law von |d/dL(T2_bd)| vs lam:")
    for i in range(len(results) - 1):
        r1, r2 = results[i], results[i+1]
        if abs(r1['dL_T2bd']) > 1e-40 and abs(r2['dL_T2bd']) > 1e-40:
            p = np.log(abs(r2['dL_T2bd']) / abs(r1['dL_T2bd'])) / np.log(r2['lam'] / r1['lam'])
            print(f"    lam={r1['lam']:.0f}->{r2['lam']:.0f}: p = {p:.3f}")

    # Versuch: Alternative Normalisierungen
    print(f"\n  (5) ALTERNATIVE NORMALISIERUNGEN")
    print(f"  d/dL(T2_bd) / (lam^p * dL_mH) fuer verschiedene p:")
    for p in [1, 2, 3, 4]:
        vals = [r['dL_T2bd'] / (r['lam']**p * r['dL_mH']) for r in results]
        cv = np.std(vals) / abs(np.mean(vals)) if abs(np.mean(vals)) > 1e-40 else float('nan')
        print(f"    p={p}: {[f'{v:+.4e}' for v in vals]}  CV={cv:.4f}")

    # Normalisierung mit T2_bd selbst
    print(f"\n  d/dL(T2_bd) / (T2_bd * lam^p) fuer verschiedene p:")
    for p in [0, 1, 2, 3]:
        vals = [r['dL_T2bd'] / (r['T2_bd'] * r['lam']**p) for r in results]
        cv = np.std(vals) / abs(np.mean(vals)) if abs(np.mean(vals)) > 1e-40 else float('nan')
        print(f"    p={p}: {[f'{v:+.4e}' for v in vals]}  CV={cv:.4f}")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
