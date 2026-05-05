"""
c2bd_scaling_analysis.py

Detaillierte Skalierungsanalyse der Derivative-Freezing-Daten.
Frage: Wie haengen d/dL(T2_bd) und d/dL(R_bulk) mit den Secular-Termen zusammen,
und welche Lambda-abhaengige Funktion vermittelt die Kopplung?
"""

import numpy as np

# Daten aus c2bd_derivative_matching.py (Multi-Lambda-Lauf)
data = [
    {
        'lam': 3.0, 'L': 2.197225,
        'dL_mua': -8.166291e-08,
        'dL_T2bd': -6.231962e-05, 'dL_Rbulk': +6.223796e-05,
        'dL_mH': +1.328829e-01, 'mHp_u0p': -3.334231e-02, 'Ct_p_Ct2': +1.127202e-01,
        'alpha': 0.654612, 'mu_alpha': 4.601596e-07,
        'T2_bd': 1.606141e-06, 'R_bulk': -1.145981e-06,
        'T1': 4.901832e-02, 'T2': -4.901786e-02,
        'Ct': None, 'mHp': None, 'dL_u0': None,
    },
    {
        'lam': 5.0, 'L': 3.218876,
        'dL_mua': +5.157073e-07,
        'dL_T2bd': -1.480290e-04, 'dL_Rbulk': +1.485447e-04,
        'dL_mH': +6.035663e-02, 'mHp_u0p': -7.314953e-03, 'Ct_p_Ct2': +5.598545e-02,
        'alpha': 0.662095, 'mu_alpha': 5.071737e-07,
        'T2_bd': 3.200568e-06, 'R_bulk': -2.693394e-06,
        'T1': 4.148265e-02, 'T2': -4.148214e-02,
    },
    {
        'lam': 7.0, 'L': 3.891820,
        'dL_mua': -1.087284e-06,
        'dL_T2bd': -2.892139e-04, 'dL_Rbulk': +2.881267e-04,
        'dL_mH': +4.199135e-02, 'mHp_u0p': -3.196922e-03, 'Ct_p_Ct2': +3.956561e-02,
        'alpha': 0.612087, 'mu_alpha': 5.696030e-07,
        'T2_bd': 5.234979e-06, 'R_bulk': -4.665376e-06,
        'T1': 3.815950e-02, 'T2': -3.815893e-02,
    },
]

lams = np.array([d['lam'] for d in data])
Ls = np.array([d['L'] for d in data])

print("=" * 80)
print("  C2bd SKALIERUNGSANALYSE")
print("=" * 80)

# 1. Power-law fits for ratios
print("\n  --- 1. Power-Law Exponent: d/dL(T2_bd) / dL_mH = C * lam^p ---")
ratio_bd = np.array([d['dL_T2bd'] / d['dL_mH'] for d in data])
print(f"  Ratios: {[f'{r:.4e}' for r in ratio_bd]}")
for i in range(len(data) - 1):
    p = np.log(abs(ratio_bd[i+1] / ratio_bd[i])) / np.log(lams[i+1] / lams[i])
    print(f"    lam={lams[i]:.0f}->{lams[i+1]:.0f}: p = {p:.3f}")

print("\n  --- 2. Power-Law Exponent: d/dL(R_bulk) / (mHp*u0p) = C * lam^p ---")
ratio_rb = np.array([d['dL_Rbulk'] / d['mHp_u0p'] for d in data])
print(f"  Ratios: {[f'{r:.4e}' for r in ratio_rb]}")
for i in range(len(data) - 1):
    p = np.log(abs(ratio_rb[i+1] / ratio_rb[i])) / np.log(lams[i+1] / lams[i])
    print(f"    lam={lams[i]:.0f}->{lams[i+1]:.0f}: p = {p:.3f}")

# 2. Normalisierung: Welche Vorfaktoren machen die Ratios lam-unabhaengig?
print("\n  --- 3. Vorfaktor-Kandidaten ---")
candidates = {
    'T2_bd': np.array([d['T2_bd'] for d in data]),
    'R_bulk': np.array([abs(d['R_bulk']) for d in data]),
    'mu/alpha': np.array([d['mu_alpha'] for d in data]),
    'T1': np.array([d['T1'] for d in data]),
    'alpha': np.array([d['alpha'] for d in data]),
    'alpha^2': np.array([d['alpha']**2 for d in data]),
    'L': Ls,
    'L^2': Ls**2,
    'lam': lams,
    'lam^2': lams**2,
    'lam^3': lams**3,
    '1/lam': 1.0 / lams,
    'lam*L': lams * Ls,
    'T2_bd/alpha': np.array([d['T2_bd'] / d['alpha'] for d in data]),
    'T1*alpha': np.array([d['T1'] * d['alpha'] for d in data]),
}

# Test: d/dL(T2_bd) / (dL_mH * candidate) = const?
print(f"\n  d/dL(T2_bd) / (dL_mH * f(lam)) = const?")
print(f"  {'Kandidat':<16} {'lam=3':>12} {'lam=5':>12} {'lam=7':>12} {'CV':>8}")
for name, vals in candidates.items():
    normed = ratio_bd / vals
    cv = np.std(normed) / abs(np.mean(normed)) if abs(np.mean(normed)) > 1e-40 else float('nan')
    print(f"  {name:<16} {normed[0]:+12.4e} {normed[1]:+12.4e} {normed[2]:+12.4e} {cv:8.4f}")

print(f"\n  d/dL(R_bulk) / (mHp*u0p * f(lam)) = const?")
print(f"  {'Kandidat':<16} {'lam=3':>12} {'lam=5':>12} {'lam=7':>12} {'CV':>8}")
for name, vals in candidates.items():
    normed = ratio_rb / vals
    cv = np.std(normed) / abs(np.mean(normed)) if abs(np.mean(normed)) > 1e-40 else float('nan')
    print(f"  {name:<16} {normed[0]:+12.4e} {normed[1]:+12.4e} {normed[2]:+12.4e} {cv:8.4f}")

# 3. Alternative: d/dL(T2_bd) = f(lam) * dL_mH + g(lam) * mHp*u0p
# Fuer jeden lambda-Punkt gibt es nur EINE Gleichung und ZWEI Unbekannte
# Aber wir koennen die Differenz und Summe testen
print("\n  --- 4. Direkte Skalierung der Ableitungen ---")
print(f"  {'lam':>5} {'|dL_T2bd|':>12} {'|dL_Rbulk|':>12} {'|dL_mH|':>12} {'|mHp_u0p|':>12} {'T2bd/Rbulk':>12}")
for d in data:
    print(f"  {d['lam']:5.1f} {abs(d['dL_T2bd']):12.4e} {abs(d['dL_Rbulk']):12.4e} "
          f"{abs(d['dL_mH']):12.4e} {abs(d['mHp_u0p']):12.4e} "
          f"{d['dL_T2bd']/d['dL_Rbulk']:+12.6f}")

# 4. Kreuz-Skalierung: d/dL(T2_bd) vs d/dL(R_bulk)
print(f"\n  --- 5. Kreuz-Verhaeltnis ---")
print(f"  d/dL(T2_bd) / d/dL(R_bulk) bei allen lam:")
for d in data:
    r = d['dL_T2bd'] / d['dL_Rbulk']
    print(f"    lam={d['lam']:.0f}: {r:+.8f}  (Abweichung von -1: {1+r:+.4e})")

# 5. Test: Ist d/dL(T2_bd) proportional zu T2_bd selbst?
# (logarithmische Ableitung: d/dL log(T2_bd) = const?)
print(f"\n  --- 6. Logarithmische Ableitung ---")
print(f"  d/dL(T2_bd) / T2_bd = ?")
for d in data:
    if abs(d['T2_bd']) > 1e-40:
        r = d['dL_T2bd'] / d['T2_bd']
        print(f"    lam={d['lam']:.0f}: {r:+.6e}  (d/dL log|T2_bd|)")
print(f"\n  d/dL(R_bulk) / R_bulk = ?")
for d in data:
    if abs(d['R_bulk']) > 1e-40:
        r = d['dL_Rbulk'] / d['R_bulk']
        print(f"    lam={d['lam']:.0f}: {r:+.6e}  (d/dL log|R_bulk|)")
print(f"\n  d/dL(mu/a) / mu/a = ?")
for d in data:
    if abs(d['mu_alpha']) > 1e-40:
        r = d['dL_mua'] / d['mu_alpha']
        print(f"    lam={d['lam']:.0f}: {r:+.6e}  (d/dL log|mu/a|)")

# 6. Resolvent-Interpretation
# Wenn m_H(u0) = -1/Ct, dann d/dL(T2_bd) koennte mit d/dL(ub_b/alpha) zusammenhaengen
# ub_b = boundary_u (Projektion von B_L-Fourierkoeffizienten auf ut)
# T2_bd = ub_b/alpha
# => d/dL(T2_bd) = d/dL(ub_b/alpha) = [d/dL(ub_b) * alpha - ub_b * d/dL(alpha)] / alpha^2
print(f"\n  --- 7. Dimensionsanalyse ---")
print(f"  d/dL(T2_bd) hat Dimension [Eigenvalue/L]")
print(f"  dL_mH hat Dimension [1/(Eigenvalue*L)]")
print(f"  => Ratio hat Dimension [Eigenvalue^2]")
print(f"  T1 = w_bar - u0 hat Dimension [Eigenvalue]")
print(f"  T1^2:")
for d in data:
    t1sq = d['T1']**2
    pred = ratio_bd[list(lams).index(d['lam'])] / t1sq
    print(f"    lam={d['lam']:.0f}: T1^2 = {t1sq:.6e}, ratio/T1^2 = {pred:.6e}")

print(f"\n  --- 8. Produkt-Ansatz: d/dL(T2_bd) = -kappa * T2_bd * dL_mH ---")
for d in data:
    if abs(d['T2_bd'] * d['dL_mH']) > 1e-40:
        kappa = -d['dL_T2bd'] / (d['T2_bd'] * d['dL_mH'])
        print(f"    lam={d['lam']:.0f}: kappa = {kappa:+.6e}")

print(f"\n  --- 9. Produkt-Ansatz: d/dL(R_bulk) = -kappa * R_bulk * mHp*u0p ---")
for d in data:
    if abs(d['R_bulk'] * d['mHp_u0p']) > 1e-40:
        kappa = -d['dL_Rbulk'] / (d['R_bulk'] * d['mHp_u0p'])
        print(f"    lam={d['lam']:.0f}: kappa = {kappa:+.6e}")

print(f"\n  --- 10. Gemischter Produkt-Ansatz: d/dL(T2_bd) = -kappa * R_bulk * dL_mH ---")
for d in data:
    if abs(d['R_bulk'] * d['dL_mH']) > 1e-40:
        kappa = -d['dL_T2bd'] / (d['R_bulk'] * d['dL_mH'])
        print(f"    lam={d['lam']:.0f}: kappa = {kappa:+.6e}")

print("\nDone.")
