"""
c2bd_highres_validation.py

Wiederholung der Derivative-Freezing-Messung bei lam=4 und lam=9
mit hoeherer Aufloesung (N=25, NGRID=30), um Trunkierungseffekte
auszuschliessen.
"""

import os, sys, time
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import mpmath as mp
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from c2bd_derivative_matching import compute_full_data, run

DPS = int(os.environ.get("DPS", 25))

# Override NGRID via monkey-patching
import c2bd_derivative_matching as mod
mod.NGRID = 30

def run_highres(lam_base, N, dlam=0.01):
    """Run mit erhoehtem N und NGRID=30."""
    L_base = 2 * np.log(lam_base)
    dL = 2 * dlam / lam_base

    print(f"\n{'='*90}")
    print(f"  HIGH-RES DERIVATIVE MATCHING: lam={lam_base}, N={N}, dlam={dlam}")
    print(f"  L={L_base:.6f}, dL={dL:.6e}, NGRID={mod.NGRID}")
    print(f"{'='*90}")

    t0 = time.time()
    print(f"  lam = {lam_base - dlam:.4f} ...", flush=True)
    d_m = compute_full_data(lam_base - dlam, N)
    print(f"    done ({time.time()-t0:.0f}s)", flush=True)

    t1 = time.time()
    print(f"  lam = {lam_base:.4f} ...", flush=True)
    d_0 = compute_full_data(lam_base, N)
    print(f"    done ({time.time()-t1:.0f}s)", flush=True)

    t2 = time.time()
    print(f"  lam = {lam_base + dlam:.4f} ...", flush=True)
    d_p = compute_full_data(lam_base + dlam, N)
    print(f"    done ({time.time()-t2:.0f}s)", flush=True)

    # Finite differences
    dL_T2bd = (d_p['T2_boundary'] - d_m['T2_boundary']) / (2 * dL)
    dL_Rbulk = (d_p['R_bulk'] - d_m['R_bulk']) / (2 * dL)
    dL_mua = (d_p['mu_alpha'] - d_m['mu_alpha']) / (2 * dL)
    dL_u0 = (d_p['u0'] - d_m['u0']) / (2 * dL)
    dL_Ct = (d_p['Ctilde'] - d_m['Ctilde']) / (2 * dL)

    # Secular
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

    print(f"\n  --- ERGEBNIS ---")
    print(f"  T2_bd = {d_0['T2_boundary']:+.6e}, R_bulk = {d_0['R_bulk']:+.6e}")
    print(f"  mu/a = {d_0['mu_alpha']:+.6e}")
    print(f"  d/dL(T2_bd)  = {dL_T2bd:+.6e}")
    print(f"  d/dL(R_bulk) = {dL_Rbulk:+.6e}")
    print(f"  d/dL(mu/a)   = {dL_mua:+.6e}")
    print(f"  Cancel-Ratio: {cancel:.4f} ({cancel*100:.1f}%)")
    print(f"  => {(1-cancel)*100:.1f}% Cancellation")
    print(f"\n  dL_mH        = {dL_mH:+.6e}")
    print(f"  mHp*u0p      = {mHp * dL_u0:+.6e}")
    print(f"  Ct'/Ct^2     = {rhs:+.6e}")
    print(f"\n  --- SCHLUESSELGROESSEN ---")
    print(f"  K_lam        = {K_lam:+.6e}")
    print(f"  T2bd/Rbulk   = {ratio:+.8f}  (Abw. von -1: {1+ratio:+.4e})")
    print(f"  dL_mua/Ct'   = {dL_mua/rhs:+.6e}")

    return {
        'lam': lam_base, 'L': L_base, 'N': N, 'NGRID': mod.NGRID,
        'K_lam': K_lam, 'cancel_pct': (1-cancel)*100,
        'ratio': ratio, 'dL_mua': dL_mua,
        'dL_T2bd': dL_T2bd, 'dL_Rbulk': dL_Rbulk,
        'dL_mH': dL_mH, 'mHp_u0p': mHp * dL_u0,
    }


def main():
    print(f"C2bd: High-Res Validation (N=25, NGRID=30)")
    print(f"DPS={DPS}")

    results = []
    for lam in [4.0, 9.0]:
        r = run_highres(lam, N=25, dlam=0.01)
        results.append(r)

    # Vergleich mit N=15 Daten
    print(f"\n{'='*90}")
    print(f"  VERGLEICH: N=15 vs N=25")
    print(f"{'='*90}")
    n15_data = {
        4.0: {'K_lam': 5.369e-6, 'cancel': 99.0, 'ratio': -0.979},
        9.0: {'K_lam': 2.917e-6, 'cancel': 94.0, 'ratio': -0.887},
    }
    print(f"  {'lam':>5} {'K(N=15)':>14} {'K(N=25)':>14} {'Cancel(15)':>12} {'Cancel(25)':>12} {'Ratio(15)':>12} {'Ratio(25)':>12}")
    for r in results:
        old = n15_data.get(r['lam'], {})
        print(f"  {r['lam']:5.1f} {old.get('K_lam', float('nan')):+14.6e} {r['K_lam']:+14.6e} "
              f"{old.get('cancel', float('nan')):12.1f}% {r['cancel_pct']:12.1f}% "
              f"{old.get('ratio', float('nan')):+12.6f} {r['ratio']:+12.6f}")

    # Erwartung: K bei guten Punkten
    K_ref = 1.902e-5
    print(f"\n  Referenz K (aus lam=3,5,7 @N=15): {K_ref:.3e}")
    for r in results:
        dev = abs(r['K_lam'] - K_ref) / K_ref * 100
        print(f"  lam={r['lam']:.0f}: K = {r['K_lam']:+.6e}, Abw. von Ref: {dev:.1f}%")


if __name__ == "__main__":
    main()
    print("\nDone.", flush=True)
