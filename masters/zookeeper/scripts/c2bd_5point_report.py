"""
c2bd_5point_report.py

Berechnet die 5 Schluesselgroessen aus dem Derivative-Freezing-Lauf
gemaess User-Spezifikation.
"""

import numpy as np

data = [
    {'lam': 3.0, 'L': 2.197225,
     'dL_T2bd': -6.231962e-05, 'dL_Rbulk': +6.223796e-05,
     'dL_mH': +1.328829e-01, 'mHp_u0p': -3.334231e-02,
     'dL_mua': -8.166291e-08, 'Ct_p_Ct2': +1.127202e-01,
     'cancel_pct': 99.9},
    {'lam': 4.0, 'L': 2.772589,
     'dL_T2bd': -2.853987e-05, 'dL_Rbulk': +2.914323e-05,
     'dL_mH': +8.305946e-02, 'mHp_u0p': -1.342888e-02,
     'dL_mua': +6.033515e-07, 'Ct_p_Ct2': +7.345309e-02,
     'cancel_pct': 99.0},
    {'lam': 5.0, 'L': 3.218876,
     'dL_T2bd': -1.480290e-04, 'dL_Rbulk': +1.485447e-04,
     'dL_mH': +6.035663e-02, 'mHp_u0p': -7.314953e-03,
     'dL_mua': +5.157073e-07, 'Ct_p_Ct2': +5.598545e-02,
     'cancel_pct': 99.8},
    {'lam': 7.0, 'L': 3.891820,
     'dL_T2bd': -2.892139e-04, 'dL_Rbulk': +2.881267e-04,
     'dL_mH': +4.199135e-02, 'mHp_u0p': -3.196922e-03,
     'dL_mua': -1.087284e-06, 'Ct_p_Ct2': +3.956561e-02,
     'cancel_pct': 99.8},
    {'lam': 9.0, 'L': 4.394449,
     'dL_T2bd': -6.958650e-05, 'dL_Rbulk': +7.848264e-05,
     'dL_mH': +3.272533e-02, 'mHp_u0p': -1.796919e-03,
     'dL_mua': +8.896142e-06, 'Ct_p_Ct2': +3.148559e-02,
     'cancel_pct': 94.0},
]

print("=" * 90)
print("  5-PUNKT-REPORT: Derivative Freezing Transferfaktor")
print("  N=15, NGRID=15, DPS=25, dlam=0.01")
print("=" * 90)

print(f"\n  (1) K_lambda := -d(T2_bd)/dL / (lam^3 * dL_mH)")
print(f"  {'lam':>5} {'lam^3':>8} {'e^(3L/2)':>12} {'K_lam':>14} {'Cancel%':>10} {'Flag':>8}")
k_vals = []
for d in data:
    lam3 = d['lam']**3
    e3L2 = np.exp(1.5 * d['L'])
    K = -d['dL_T2bd'] / (lam3 * d['dL_mH'])
    k_vals.append(K)
    flag = "OK" if d['cancel_pct'] >= 99.5 else "SCHWACH" if d['cancel_pct'] >= 97 else "!!TRUNC"
    print(f"  {d['lam']:5.1f} {lam3:8.0f} {e3L2:12.2f} {K:+14.6e} {d['cancel_pct']:10.1f}% {flag:>8}")
k_arr = np.array(k_vals)
k_good = [k for k, d in zip(k_vals, data) if d['cancel_pct'] >= 99.5]
if len(k_good) >= 2:
    k_good = np.array(k_good)
    print(f"\n  Nur lam mit Cancel >= 99.5%: K = {np.mean(k_good):+.6e} +/- {np.std(k_good):.2e}")
    print(f"  CV = {np.std(k_good)/abs(np.mean(k_good)):.4f} ({np.std(k_good)/abs(np.mean(k_good))*100:.1f}%)")
print(f"  Alle 5 Punkte: K = {np.mean(k_arr):+.6e} +/- {np.std(k_arr):.2e}")
print(f"  CV = {np.std(k_arr)/abs(np.mean(k_arr)):.4f}")

print(f"\n  (2) d(T2_bd)/dL / d(R_bulk)/dL  — Locking auf -1")
print(f"  {'lam':>5} {'Ratio':>14} {'Abw. von -1':>14}")
for d in data:
    r = d['dL_T2bd'] / d['dL_Rbulk']
    print(f"  {d['lam']:5.1f} {r:+14.8f} {1+r:+14.6e}")

print(f"\n  (3) d/dL(mu/a) / (Ct'/Ct^2)  — Secular controls residual")
print(f"  {'lam':>5} {'Ratio':>14}")
for d in data:
    r = d['dL_mua'] / d['Ct_p_Ct2']
    print(f"  {d['lam']:5.1f} {r:+14.6e}")

print(f"\n  (4) Darstellung in L:  lam^3 = e^(3L/2)")
print(f"  {'lam':>5} {'L':>10} {'lam^3':>10} {'e^(3L/2)':>12} {'Delta':>12}")
for d in data:
    lam3 = d['lam']**3
    e3L2 = np.exp(1.5 * d['L'])
    print(f"  {d['lam']:5.1f} {d['L']:10.6f} {lam3:10.1f} {e3L2:12.1f} {abs(lam3-e3L2):12.2e}")

print(f"\n  (5) E_lambda := d/dL(R_bulk) + K_mean * lam^3 * dL_mH")
# Verwende K_mean nur aus guten Punkten
K_mean = np.mean(k_good) if len(k_good) >= 2 else np.mean(k_arr)
print(f"      (K_mean = {K_mean:+.6e}, berechnet aus lam=3,5,7)")
print(f"  {'lam':>5} {'d/dL(R_bulk)':>14} {'K*lam^3*dL_mH':>14} {'E_lam':>14} {'E/d_Rbulk':>12}")
for d in data:
    pred = K_mean * d['lam']**3 * d['dL_mH']
    E = d['dL_Rbulk'] + pred
    rel = E / d['dL_Rbulk'] if abs(d['dL_Rbulk']) > 1e-40 else float('nan')
    print(f"  {d['lam']:5.1f} {d['dL_Rbulk']:+14.6e} {pred:+14.6e} {E:+14.6e} {rel:+12.4e}")

# Zusatz: Diagnosik ob N=15 genuegt
print(f"\n  --- TRUNCATION-DIAGNOSTIK ---")
print(f"  {'lam':>5} {'L':>8} {'L/N':>8} {'NGRID/L':>8} {'dL_T1':>14} {'|dL_mua|':>14}")
for d in data:
    dL_T1 = None
    print(f"  {d['lam']:5.1f} {d['L']:8.3f} {d['L']/15:8.4f} {15/d['L']:8.2f} {'—':>14} {abs(d['dL_mua']):14.4e}")

print(f"\n  FAZIT:")
good_lams = [d['lam'] for d in data if d['cancel_pct'] >= 99.5]
weak_lams = [d['lam'] for d in data if d['cancel_pct'] < 99.5]
print(f"  Robuste Punkte (Cancel >= 99.5%): lam = {good_lams}")
print(f"  Schwache Punkte: lam = {weak_lams}")
if weak_lams:
    print(f"  => lam={weak_lams} muessen mit groesserem N und NGRID wiederholt werden")
print(f"  Lambda^3-Hypothese (aus robusten Punkten): K = {np.mean(k_good):+.6e}, CV = {np.std(k_good)/abs(np.mean(k_good))*100:.1f}%")

print("\nDone.")
