"""
FST Dark Energy w(z) vs DESI Constraints
=========================================
Implementiert die tanh-Saettigungsdynamik aus dem FST-DE Paper und
vergleicht sie mit den DESI DR2 (2025) Constraints in der CPL-Parametrisierung.

Paper-Formel (Eq. weff):
    Omega_Phi(a) = Phi_0 * tanh[kappa * (a - a_trans)]

    w_eff(a) = -1 - (1/3) * d(ln Omega_Phi)/d(ln a)

    d(ln Omega_Phi)/d(ln a) = a*kappa / (Phi_0 * tanh[...]) * sech^2[...]

Paper-Parameter (illustrative, Sec. 4.2):
    kappa = 1.5 bis 2.0 (Mittelwert 1.75)
    a_trans = 0.50
    Phi_0 = 0.685

DESI DR2 (2025) - Werte aus dem Paper (Zeile 818-819):
    w0 = -0.42 +/- 0.21
    wa = -1.79 +/- 0.67

DESI Year 1 (2024) BAO+CMB+PantheonPlus (als Referenz):
    w0 = -0.727 +/- 0.067
    wa = -1.05 (+0.31/-0.27)
"""

import sys
import os
import numpy as np

# Sicherstellen dass matplotlib verfuegbar ist
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MATPLOTLIB = True
except ImportError:
    print("matplotlib nicht verfuegbar. Installiere...")
    os.system(f"{sys.executable} -m pip install matplotlib")
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MATPLOTLIB = True

# ============================================================
# Parameter
# ============================================================

# FST-Modell Paper-Parameter (illustrative values)
KAPPA_LOW  = 1.5
KAPPA_MID  = 1.75   # Mittelwert
KAPPA_HIGH = 2.0
A_TRANS    = 0.50
PHI_0      = 0.685  # = Omega_Lambda

# DESI DR2 (2025) -- Werte direkt aus dem Paper
# "DESI DR2 analysis: w0 = -0.42 +/- 0.21, wa = -1.79 +/- 0.67"
DESI_W0_DR2    = -0.42
DESI_W0_ERR_DR2 = 0.21
DESI_WA_DR2    = -1.79
DESI_WA_ERR_DR2 = 0.67

# DESI Year 1 (2024) BAO+CMB+PantheonPlus (Referenzwerte aus Aufgabe)
DESI_W0_Y1    = -0.727
DESI_W0_ERR_Y1 = 0.067
DESI_WA_Y1    = -1.05
DESI_WA_ERR_Y1 = 0.29  # Symmetrisiert: (0.31 + 0.27) / 2

# ============================================================
# FST-Modell: w_eff(z)
# ============================================================

def omega_phi(a, kappa=KAPPA_MID, phi_0=PHI_0, a_trans=A_TRANS):
    """Saettigungs-Loesung fuer Omega_Phi(a)."""
    return phi_0 * np.tanh(kappa * (a - a_trans))

def w_eff_fst(a, kappa=KAPPA_MID, phi_0=PHI_0, a_trans=A_TRANS, eps=1e-6):
    """
    Effektive Zustandsgleichung aus dem FST-Paper (Eq. weff):
        w_eff(a) = -1 - (1/3) * d(ln Omega_Phi)/d(ln a)

    d(ln Omega_Phi)/d(ln a) = a * kappa * sech^2[kappa*(a-a_trans)]
                               / (phi_0 * tanh[kappa*(a-a_trans)])

    HINWEIS: Bei a = a_trans ist tanh = 0 (Singularitaet).
    In diesem Bereich gilt Omega_Phi = 0 (Materie-DE-Gleichgewicht).
    Das Modell ist hier nicht definiert -- wir maskieren.
    Fuer a < a_trans ist Omega_Phi < 0 (unphysikalisch im Sinne der DE).
    Physikalisch relevanter Bereich: a > a_trans (d.h. z < z_trans = 1/a_trans - 1).
    """
    arg = kappa * (a - a_trans)
    tanh_val = np.tanh(arg)
    sech2_val = 1.0 - tanh_val**2  # = sech^2(arg)

    # Maske: Bereiche wo |tanh| < eps als NaN markieren
    valid = np.abs(tanh_val) > eps
    safe_tanh = np.where(valid, tanh_val, np.nan)

    dlnO_dlna = a * kappa * sech2_val / (phi_0 * safe_tanh)
    result = -1.0 - dlnO_dlna / 3.0

    # Unphysikalische Bereiche (a << a_trans, Omega_Phi << 0) ebenfalls maskieren
    # Physikalisch relevant: a > 0.3 (z < 2.33) und Omega_Phi > 0.01
    omega = omega_phi(a, kappa=kappa, phi_0=phi_0, a_trans=a_trans)
    result = np.where(omega > 0.01, result, np.nan)

    return result

def w_eff_fst_at_z(z, kappa=KAPPA_MID, phi_0=PHI_0, a_trans=A_TRANS):
    """w_eff als Funktion von z (z = 1/a - 1, also a = 1/(1+z))."""
    a = 1.0 / (1.0 + z)
    return w_eff_fst(a, kappa=kappa, phi_0=phi_0, a_trans=a_trans)

# ============================================================
# CPL-Parametrisierung: w(z) = w0 + wa * z/(1+z)
# ============================================================

def w_cpl(z, w0, wa):
    """CPL-Parametrisierung (Chevallier-Polarski-Linder)."""
    return w0 + wa * z / (1.0 + z)

# ============================================================
# CPL-Fit an FST: Linearer Fit bei z=0
# ============================================================

def fst_to_cpl(kappa=KAPPA_MID, phi_0=PHI_0, a_trans=A_TRANS):
    """
    Linearisierung des FST-Modells bei z=0 in CPL-Parameter:
        w0 = w_eff(z=0)
        wa = -dw/dz |_{z=0}  (CPL: w(z) = w0 + wa*z/(1+z), dw/dz|0 = wa)

    Numerische Ableitung.
    """
    z0 = 1e-4
    w0_fst = float(w_eff_fst_at_z(0.0, kappa=kappa, phi_0=phi_0, a_trans=a_trans))
    w_dz   = float(w_eff_fst_at_z(z0,  kappa=kappa, phi_0=phi_0, a_trans=a_trans))
    # CPL: dw/dz|_{z=0} = wa (da z/(1+z) -> z fuer kleine z)
    wa_fst = (w_dz - w0_fst) / z0
    return w0_fst, wa_fst

# ============================================================
# Kompatibilitaet: Chi2-artige Metrik
# ============================================================

def compatibility_sigma(w0_fst, wa_fst,
                         w0_ref, w0_err, wa_ref, wa_err):
    """
    Einfache quadratische Abweichung (ohne Kovarianz):
        delta_sigma = sqrt( ((w0_fst - w0_ref)/w0_err)^2
                          + ((wa_fst - wa_ref)/wa_err)^2 )
    """
    delta_w0 = (w0_fst - w0_ref) / w0_err
    delta_wa = (wa_fst - wa_ref) / wa_err
    return np.sqrt(delta_w0**2 + delta_wa**2), delta_w0, delta_wa

# ============================================================
# Hauptberechnung
# ============================================================

def main():
    print("=" * 65)
    print("FST Dark Energy w(z) vs DESI Constraints")
    print("=" * 65)

    # z-Array
    z_arr = np.linspace(0.0, 3.0, 600)

    # --- FST-Modell (3 kappa-Werte) ---
    w_fst_mid  = w_eff_fst_at_z(z_arr, kappa=KAPPA_MID)
    w_fst_low  = w_eff_fst_at_z(z_arr, kappa=KAPPA_LOW)
    w_fst_high = w_eff_fst_at_z(z_arr, kappa=KAPPA_HIGH)

    # --- CPL-Modelle ---
    # DESI DR2
    w_desi_dr2   = w_cpl(z_arr, DESI_W0_DR2, DESI_WA_DR2)
    # 1-sigma Bänder: Fehlerfortpflanzung
    # sigma_w = sqrt( (dw/dw0)^2 * sigma_w0^2 + (dw/dwa)^2 * sigma_wa^2 )
    #         = sqrt( 1 * sigma_w0^2 + (z/(1+z))^2 * sigma_wa^2 )
    f_z = z_arr / (1.0 + z_arr)
    w_desi_dr2_1s_lo = w_desi_dr2 - np.sqrt(DESI_W0_ERR_DR2**2 + (f_z * DESI_WA_ERR_DR2)**2)
    w_desi_dr2_1s_hi = w_desi_dr2 + np.sqrt(DESI_W0_ERR_DR2**2 + (f_z * DESI_WA_ERR_DR2)**2)
    w_desi_dr2_2s_lo = w_desi_dr2 - 2 * np.sqrt(DESI_W0_ERR_DR2**2 + (f_z * DESI_WA_ERR_DR2)**2)
    w_desi_dr2_2s_hi = w_desi_dr2 + 2 * np.sqrt(DESI_W0_ERR_DR2**2 + (f_z * DESI_WA_ERR_DR2)**2)

    # DESI Year 1
    w_desi_y1   = w_cpl(z_arr, DESI_W0_Y1, DESI_WA_Y1)
    w_desi_y1_1s_lo = w_desi_y1 - np.sqrt(DESI_W0_ERR_Y1**2 + (f_z * DESI_WA_ERR_Y1)**2)
    w_desi_y1_1s_hi = w_desi_y1 + np.sqrt(DESI_W0_ERR_Y1**2 + (f_z * DESI_WA_ERR_Y1)**2)
    w_desi_y1_2s_lo = w_desi_y1 - 2 * np.sqrt(DESI_W0_ERR_Y1**2 + (f_z * DESI_WA_ERR_Y1)**2)
    w_desi_y1_2s_hi = w_desi_y1 + 2 * np.sqrt(DESI_W0_ERR_Y1**2 + (f_z * DESI_WA_ERR_Y1)**2)

    # --- CPL-Linearisierung des FST-Modells ---
    w0_fst_mid,  wa_fst_mid  = fst_to_cpl(kappa=KAPPA_MID)
    w0_fst_low,  wa_fst_low  = fst_to_cpl(kappa=KAPPA_LOW)
    w0_fst_high, wa_fst_high = fst_to_cpl(kappa=KAPPA_HIGH)

    # --- Kompatibilitaet mit DESI DR2 ---
    sigma_dr2, dw0_dr2, dwa_dr2 = compatibility_sigma(
        w0_fst_mid, wa_fst_mid,
        DESI_W0_DR2, DESI_W0_ERR_DR2,
        DESI_WA_DR2, DESI_WA_ERR_DR2
    )
    # --- Kompatibilitaet mit DESI Year 1 ---
    sigma_y1, dw0_y1, dwa_y1 = compatibility_sigma(
        w0_fst_mid, wa_fst_mid,
        DESI_W0_Y1, DESI_W0_ERR_Y1,
        DESI_WA_Y1, DESI_WA_ERR_Y1
    )

    # ============================================================
    # Plot
    # ============================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("FST Dark Energy $w(z)$ vs DESI Constraints", fontsize=14, fontweight='bold')

    # ---- Panel 1: vs DESI DR2 (2025) ----
    ax1 = axes[0]

    # 2-sigma Band
    ax1.fill_between(z_arr, w_desi_dr2_2s_lo, w_desi_dr2_2s_hi,
                     alpha=0.15, color='royalblue', label='DESI DR2 $2\\sigma$')
    # 1-sigma Band
    ax1.fill_between(z_arr, w_desi_dr2_1s_lo, w_desi_dr2_1s_hi,
                     alpha=0.30, color='royalblue', label='DESI DR2 $1\\sigma$')
    # DESI best-fit
    ax1.plot(z_arr, w_desi_dr2, color='royalblue', lw=2.0, ls='--',
             label=f'DESI DR2 best-fit ($w_0={DESI_W0_DR2}$, $w_a={DESI_WA_DR2}$)')

    # FST-Modell (kappa-Bereich als Band)
    ax1.fill_between(z_arr, w_fst_low, w_fst_high,
                     alpha=0.25, color='crimson', label='FST ($\\kappa=1.5\\text{--}2.0$)')
    ax1.plot(z_arr, w_fst_mid, color='crimson', lw=2.5,
             label=f'FST best ($\\kappa={KAPPA_MID}$, $a_{{\\rm trans}}={A_TRANS}$, $\\Phi_0={PHI_0}$)')

    # Lambda-CDM
    ax1.axhline(-1.0, color='gray', ls=':', lw=1.5, label='$\\Lambda$CDM ($w=-1$)')

    ax1.set_xlabel('Redshift $z$', fontsize=12)
    ax1.set_ylabel('$w(z)$', fontsize=12)
    ax1.set_title('vs DESI DR2 (2025)', fontsize=11)
    ax1.set_xlim(0, 3)
    ax1.set_ylim(-4.0, 0.5)
    ax1.legend(fontsize=7.5, loc='lower left')
    ax1.grid(True, alpha=0.3)

    # ---- Panel 2: vs DESI Year 1 (2024) ----
    ax2 = axes[1]

    ax2.fill_between(z_arr, w_desi_y1_2s_lo, w_desi_y1_2s_hi,
                     alpha=0.15, color='steelblue', label='DESI Y1 $2\\sigma$')
    ax2.fill_between(z_arr, w_desi_y1_1s_lo, w_desi_y1_1s_hi,
                     alpha=0.30, color='steelblue', label='DESI Y1 $1\\sigma$')
    ax2.plot(z_arr, w_desi_y1, color='steelblue', lw=2.0, ls='--',
             label=f'DESI Y1 best-fit ($w_0={DESI_W0_Y1}$, $w_a={DESI_WA_Y1}$)')

    ax2.fill_between(z_arr, w_fst_low, w_fst_high,
                     alpha=0.25, color='crimson', label='FST ($\\kappa=1.5\\text{--}2.0$)')
    ax2.plot(z_arr, w_fst_mid, color='crimson', lw=2.5,
             label=f'FST best ($\\kappa={KAPPA_MID}$)')

    ax2.axhline(-1.0, color='gray', ls=':', lw=1.5, label='$\\Lambda$CDM ($w=-1$)')

    ax2.set_xlabel('Redshift $z$', fontsize=12)
    ax2.set_ylabel('$w(z)$', fontsize=12)
    ax2.set_title('vs DESI Year 1 (2024)', fontsize=11)
    ax2.set_xlim(0, 3)
    ax2.set_ylim(-4.0, 0.5)
    ax2.legend(fontsize=7.5, loc='lower left')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Speicherpfad
    out_path = r"C:\Users\User\OneDrive\.RESEARCH\Natur&Technik\3 Folgebeweise\Dark Energy\compute_w_vs_desi.png"
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"\nPlot gespeichert: {out_path}")
    plt.close()

    # ============================================================
    # Zusammenfassung
    # ============================================================
    print("\n" + "=" * 65)
    print("ZUSAMMENFASSUNG")
    print("=" * 65)

    print("\n--- FST-Modell Parameter (Paper, illustrativ) ---")
    print(f"  kappa        = {KAPPA_MID} (Bereich: {KAPPA_LOW}--{KAPPA_HIGH})")
    print(f"  a_trans      = {A_TRANS}")
    print(f"  Phi_0        = {PHI_0}")

    def fmt_w(val):
        v = float(val)
        if np.isnan(v):
            return "NaN (Singularitaet nahe a_trans oder Omega_Phi<0)"
        return f"{v:.4f}"

    print("\n--- FST w_eff Werte (kappa = 1.75) ---")
    print(f"  w_eff(z=0)   = {fmt_w(w_eff_fst_at_z(0.0))}  [Paper: -1.4 bis -1.5]")
    print(f"  w_eff(z=0.5) = {fmt_w(w_eff_fst_at_z(0.5))}")
    print(f"  w_eff(z=1.0) = {fmt_w(w_eff_fst_at_z(1.0))}  [= a_trans: Singularitaet]")
    print(f"  w_eff(z=2.0) = {fmt_w(w_eff_fst_at_z(2.0))}")
    print(f"  w_eff(z=3.0) = {fmt_w(w_eff_fst_at_z(3.0))}")

    print("\n--- CPL-Linearisierung des FST-Modells (bei z=0) ---")
    print(f"  FST w0 (kappa=1.5)  = {w0_fst_low:.4f}")
    print(f"  FST w0 (kappa=1.75) = {w0_fst_mid:.4f}")
    print(f"  FST w0 (kappa=2.0)  = {w0_fst_high:.4f}")
    print(f"  FST wa (kappa=1.75) = {wa_fst_mid:.4f}  [dw/dz|_{{z=0}}]")

    print("\n--- DESI DR2 (2025) --- [Werte aus dem Paper]")
    print(f"  w0 = {DESI_W0_DR2} +/- {DESI_W0_ERR_DR2}")
    print(f"  wa = {DESI_WA_DR2} +/- {DESI_WA_ERR_DR2}")

    print("\n--- DESI Year 1 (2024) ---")
    print(f"  w0 = {DESI_W0_Y1} +/- {DESI_W0_ERR_Y1}")
    print(f"  wa = {DESI_WA_Y1} +/- {DESI_WA_ERR_Y1}")

    print("\n--- Kompatibilitaet FST vs DESI DR2 ---")
    print(f"  Abweichung w0: {dw0_dr2:+.2f} sigma")
    print(f"  Abweichung wa: {dwa_dr2:+.2f} sigma")
    print(f"  Gesamt:        {sigma_dr2:.2f} sigma")
    compat_dr2 = "JA (< 2 sigma)" if sigma_dr2 < 2.0 else ("GRENZWERTIG (2-3 sigma)" if sigma_dr2 < 3.0 else "NEIN (> 3 sigma)")
    print(f"  Kompatibel:    {compat_dr2}")

    print("\n--- Kompatibilitaet FST vs DESI Year 1 ---")
    print(f"  Abweichung w0: {dw0_y1:+.2f} sigma")
    print(f"  Abweichung wa: {dwa_y1:+.2f} sigma")
    print(f"  Gesamt:        {sigma_y1:.2f} sigma")
    compat_y1 = "JA (< 2 sigma)" if sigma_y1 < 2.0 else ("GRENZWERTIG (2-3 sigma)" if sigma_y1 < 3.0 else "NEIN (> 3 sigma)")
    print(f"  Kompatibel:    {compat_y1}")

    print("\n--- Hinweis ---")
    print("  Die illustrativen Paper-Parameter sind noch nicht")
    print("  an SN Ia + BAO Daten gefittet. Ein dedizierter")
    print("  numerischer Fit (MCMC) koennte die Kompatibilitaet")
    print("  deutlich verbessern. Das Paper verweist selbst auf")
    print("  diesen offenen Punkt (Sec. 4.2, nach Gl. w0-prediction).")
    print("  WICHTIG: w_eff(z) != w_DE(z). Das Paper unterscheidet")
    print("  klar zwischen dem effektiven und dem intrinsischen w.")
    print("  Der Vergleich mit CPL ist daher nur orientierend.")

    print("\n" + "=" * 65)
    print("Berechnung abgeschlossen.")
    print("=" * 65)


if __name__ == "__main__":
    main()
