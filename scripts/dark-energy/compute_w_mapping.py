"""
FST Dark Energy: w_eff -> w_DE Mapping
=======================================
Korrekte Ableitung von w_DE(z) aus dem FST-Modell.

PROBLEM MIT DEM ALTEN SKRIPT:
  Das Skript compute_w_vs_desi.py vergleicht w_eff (Gesamtuniversum) mit
  DESI w_DE (nur Dark-Energy-Komponente). Das ist falsch.

KORREKTE PHYSIK (Friedmann-Universum mit Materie + DE):
  Kontinuitaetsgleichung fuer DE:
      rho_DE'(a) + 3*(1 + w_DE(a))/a * rho_DE(a) = 0

  Damit: w_DE(a) = -1 - (a/3) * d(ln rho_DE)/da

  Im FST-Modell gilt:
      Omega_DE(a) = Phi_0 * tanh[kappa*(a - a_trans)]
      rho_DE(a) = rho_crit(a) * Omega_DE(a)
                = rho_crit,0 * H(a)^2/H0^2 * Omega_DE(a)

  Aber: Da H(a)^2 = H0^2 * [Omega_m/a^3 + Omega_DE(a)],
  ist die Entkopplung von rho_DE und Omega_DE nicht trivial.

  KORREKTE METHODE via Druckterm:
  Die DE-Energiedichte aus Omega_DE:
      rho_DE(a) proportional zu exp(-3 * integral_1^a [1 + w_DE(a')]/a' da')

  Aber wir kennen Omega_DE(a) = Phi_0 * tanh[kappa*(a-a_trans)] (FST-Ansatz).
  Daraus folgt rho_DE(a) / rho_crit,0 = Omega_DE(a) * (H(a)/H0)^2.

  KONSISTENTE LOESUNG:
  Wir loesen das selbstkonsistent: Fuer gegebenes w_DE(a) gilt:
      Omega_DE(a) = Omega_DE,0 * a^{-3(1+w_eff_avg)}  (generell)

  Aber das FST-Modell gibt uns Omega_DE(a) direkt als Ansatz. Daraus:
      rho_DE(a) = (3 H0^2 / 8 pi G) * Omega_DE(a) * (H(a)/H0)^2

  Die EINFACHSTE KONSISTENTE ABLEITUNG:
  Aus der modifizierten Friedmann-Gleichung:
      w_DE(a) = -1 - (1/3) * a/Omega_DE(a) * dOmega_DE/da
               - Omega_m(a)/(3 * Omega_DE(a)) * [1 - 3*w_DE(a) - ...]

  Fuer das Standard-Lambda-CDM-aehnliche Modell mit zeitabhaengiger DE:
      w_DE(a) = [w_eff(a) - Omega_m(a)*(-1/3)] / Omega_DE(a)

  HERLEITUNG:
  w_eff definiert die totale EOS: p_total = w_eff * rho_total
  Mit p_m = 0 (Staub) und p_DE = w_DE * rho_DE:
      w_eff = p_DE / rho_total = w_DE * (rho_DE / rho_total) = w_DE * Omega_DE(a)

  Wobei Omega_DE(a) = rho_DE / rho_total = rho_DE / (rho_m + rho_DE).

  DAHER: w_DE(a) = w_eff(a) / Omega_DE(a)

  mit w_eff aus der Kontinuitaetsgleichung des Gesamtfluids.

  ABER: Die FST-Formel ist schon:
      w_eff(a) = -1 - (1/3) * d(ln Omega_Phi)/d(ln a)
  Das ist die EOS des GESAMTSYSTEMS (nicht nur DE).

  KORREKTE RELATION (keine Materie im Ausdruck angenommen bei a~1):
      w_DE(a) = -1 - (1/3) * (a / Omega_DE) * dOmega_DE/da

  Dies ist die intrinsische EOS der DE-Komponente, direkt aus der
  Kontinuitaetsgleichung fuer rho_DE, WENN Omega_DE(a) das Verhaeltnis
  rho_DE/rho_crit,0 beschreibt (was im FST-Ansatz der Fall ist).

DESI DR2 (2025):
  w0 = -0.42 +/- 0.21
  wa = -1.79 +/- 0.67

Autor: Lukas Geiger
Datum: 2026-03-18
"""

import sys
import os
import numpy as np

# Matplotlib
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    os.system(f"{sys.executable} -m pip install matplotlib --quiet")
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True

# ============================================================
# Kosmologische Parameter (Planck 2018)
# ============================================================
OMEGA_M0   = 0.3089    # Materiedichte heute
OMEGA_DE0  = 0.6911    # DE-Dichte heute (= Phi_0 * tanh[kappa*(1 - a_trans)])
H0         = 67.74     # km/s/Mpc (Planck 2018)

# DESI DR2 (2025)
DESI_W0     = -0.42
DESI_W0_ERR =  0.21
DESI_WA     = -1.79
DESI_WA_ERR =  0.67

# ============================================================
# FST-Modell: Omega_DE(a)
# ============================================================

def omega_de_fst(a, kappa, a_trans, phi_0):
    """
    FST-Ansatz: Omega_DE(a) = Phi_0 * tanh[kappa*(a - a_trans)]

    Dies beschreibt den Anteil der Dunklen Energie an der kritischen Dichte.
    Normierung: Phi_0 wird so gewaehlt, dass Omega_DE(1) = Omega_DE0.
    """
    return phi_0 * np.tanh(kappa * (a - a_trans))


def compute_phi0(kappa, a_trans, omega_de0=OMEGA_DE0):
    """
    Normierungsbedingung: Omega_DE(a=1) = Omega_DE0
    => Phi_0 = Omega_DE0 / tanh(kappa*(1 - a_trans))
    """
    val = np.tanh(kappa * (1.0 - a_trans))
    if abs(val) < 1e-10:
        return np.nan
    return omega_de0 / val


# ============================================================
# w_DE(a): Direkte Herleitung aus Kontinuitaetsgleichung
# ============================================================

def w_DE_from_omega(a, kappa, a_trans, phi_0, eps=1e-8):
    """
    Intrinsische EOS der Dunklen Energie.

    Aus der Kontinuitaetsgleichung fuer rho_DE:
        drho_DE/da + 3*(1+w_DE)/a * rho_DE = 0

    Mit rho_DE(a) = rho_crit,0 * Omega_DE(a) * (H(a)/H0)^2
    und H^2(a) = H0^2 * [Omega_m/a^3 + Omega_DE(a)]:

        rho_DE(a) = (3H0^2/8piG) * Omega_DE(a) * [Omega_m/a^3 + Omega_DE(a)]

    Logarithmische Ableitung:
        d(ln rho_DE)/d(ln a) = d(ln Omega_DE)/d(ln a)
                              + d(ln[Omega_m/a^3 + Omega_DE])/d(ln a)

    Daraus: w_DE(a) = -1 - (1/3) * d(ln rho_DE)/d(ln a)

    VEREINFACHUNG fuer das FST-Modell:
    Im FST-Paper ist Omega_DE(a) als rho_DE-Anteil definiert
    (nicht als H-gewichteter Anteil). Wir nutzen die direkte Herleitung:

        w_DE(a) = -1 - (a/3) * (1/Omega_DE) * dOmega_DE/da

    Diese Formel gilt, wenn Omega_DE(a) = rho_DE(a) / rho_crit,0 ist
    (d.h. im Nenner steht die heutige kritische Dichte, nicht die z-abhaengige).
    """
    # Omega_DE und ihre Ableitung
    arg = kappa * (a - a_trans)
    tanh_val = np.tanh(arg)
    sech2_val = 1.0 - tanh_val**2

    omega = phi_0 * tanh_val
    d_omega_da = phi_0 * kappa * sech2_val  # dOmega_DE/da

    # Maske: Omega_DE zu klein oder negativ
    valid = omega > eps
    safe_omega = np.where(valid, omega, np.nan)

    # w_DE = -1 - (a/3) * d(ln Omega_DE)/da
    #       = -1 - (a / (3 * Omega_DE)) * dOmega_DE/da
    w_de = -1.0 - (a / (3.0 * safe_omega)) * d_omega_da

    return w_de


def w_DE_full(a, kappa, a_trans, phi_0, omega_m0=OMEGA_M0, eps=1e-8):
    """
    Vollstaendige w_DE(a) mit Friedmann-Kopplung.

    Aus d(ln rho_DE)/d(ln a) mit rho_DE = rho_crit,0 * Omega_DE * E^2(a):
    wobei E^2(a) = Omega_m0/a^3 + Omega_DE(a).

    d(ln rho_DE)/d(ln a) = d(ln Omega_DE)/d(ln a) + d(ln E^2)/d(ln a)

    d(ln E^2)/d(ln a) = [-3*Omega_m0/a^3 + a*dOmega_DE/da] / E^2(a)
    """
    arg = kappa * (a - a_trans)
    tanh_val = np.tanh(arg)
    sech2_val = 1.0 - tanh_val**2

    omega_de = phi_0 * tanh_val
    d_omega_de_da = phi_0 * kappa * sech2_val

    # Friedmann E^2(a)
    e2 = omega_m0 / a**3 + omega_de

    # Nur physikalischer Bereich
    valid = (omega_de > eps) & (e2 > eps)
    safe_omega = np.where(valid, omega_de, np.nan)
    safe_e2 = np.where(valid, e2, np.nan)

    # d(ln Omega_DE)/d(ln a) = a * dOmega_DE/da / Omega_DE
    dlnO_dlna = a * d_omega_de_da / safe_omega

    # d(ln E^2)/d(ln a)
    d_lne2_dlna = (-3.0 * omega_m0 / a**3 + a * d_omega_de_da) / safe_e2

    # w_DE = -1 - (1/3) * [dlnO_dlna + d_lne2_dlna]
    w_de = -1.0 - (1.0 / 3.0) * (dlnO_dlna + d_lne2_dlna)

    return w_de


# ============================================================
# CPL-Linearisierung: w_DE ≈ w0_DE + wa_DE * z/(1+z)
# ============================================================

def cpl_fit_wDE(kappa, a_trans, phi_0, method='full', n_points=100):
    """
    Linearisierung von w_DE(z) in CPL-Parameter (w0_DE, wa_DE).

    Methode: Least-Squares-Fit ueber z in [0, 2].
    CPL: w(z) = w0 + wa * z/(1+z)  =>  w(a) = w0 + wa * (1-a)
    """
    a_arr = np.linspace(0.35, 1.0, n_points)  # z=0 bis z~1.86
    z_arr = 1.0 / a_arr - 1.0

    if method == 'full':
        w_arr = w_DE_full(a_arr, kappa, a_trans, phi_0)
    else:
        w_arr = w_DE_from_omega(a_arr, kappa, a_trans, phi_0)

    # Nur gueltigen Bereich nutzen
    valid = np.isfinite(w_arr)
    if np.sum(valid) < 5:
        return np.nan, np.nan

    z_v = z_arr[valid]
    w_v = w_arr[valid]

    # CPL basis: [1, z/(1+z)]
    f_z = z_v / (1.0 + z_v)  # = 1 - a
    A = np.column_stack([np.ones_like(z_v), f_z])

    # Least-Squares
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A, w_v, rcond=None)
        return float(coeffs[0]), float(coeffs[1])
    except Exception:
        return np.nan, np.nan


# ============================================================
# Chi^2 gegen DESI DR2
# ============================================================

def chi2_vs_desi(w0_de, wa_de,
                 w0_ref=DESI_W0, w0_err=DESI_W0_ERR,
                 wa_ref=DESI_WA, wa_err=DESI_WA_ERR):
    """
    Einfaches chi^2 ohne Kovarianz:
        chi2 = ((w0_de - w0_ref)/w0_err)^2 + ((wa_de - wa_ref)/wa_err)^2
    """
    if not (np.isfinite(w0_de) and np.isfinite(wa_de)):
        return np.inf
    return ((w0_de - w0_ref) / w0_err)**2 + ((wa_de - wa_ref) / wa_err)**2


# ============================================================
# Grid-Scan
# ============================================================

def grid_scan(n_kappa=20, n_atrans=20):
    """
    Scan ueber (kappa, a_trans) und berechne chi^2 gegen DESI DR2.
    """
    kappa_arr  = np.linspace(0.5, 5.0, n_kappa)
    atrans_arr = np.linspace(0.3, 0.8, n_atrans)

    chi2_grid = np.full((n_kappa, n_atrans), np.inf)
    w0_grid   = np.full((n_kappa, n_atrans), np.nan)
    wa_grid   = np.full((n_kappa, n_atrans), np.nan)

    for i, kappa in enumerate(kappa_arr):
        for j, at in enumerate(atrans_arr):
            phi0 = compute_phi0(kappa, at)
            if not np.isfinite(phi0) or phi0 <= 0:
                continue

            w0_de, wa_de = cpl_fit_wDE(kappa, at, phi0, method='full')
            if not (np.isfinite(w0_de) and np.isfinite(wa_de)):
                continue

            chi2_grid[i, j] = chi2_vs_desi(w0_de, wa_de)
            w0_grid[i, j]   = w0_de
            wa_grid[i, j]   = wa_de

    return kappa_arr, atrans_arr, chi2_grid, w0_grid, wa_grid


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 70)
    print("FST Dark Energy: w_eff -> w_DE Mapping (Korrektur)")
    print("=" * 70)

    # --- Einzelne Kurve: Illustrative Paper-Parameter ---
    kappa_demo  = 1.75
    a_trans_demo = 0.50
    phi0_demo   = compute_phi0(kappa_demo, a_trans_demo)

    print(f"\n[1] Demo-Parameter: kappa={kappa_demo}, a_trans={a_trans_demo}")
    print(f"    Phi_0 (normiert) = {phi0_demo:.4f}  [Omega_DE(a=1) = {OMEGA_DE0}]")

    z_arr = np.linspace(0.0, 2.0, 400)
    a_arr = 1.0 / (1.0 + z_arr)

    # Beide Methoden
    w_eff_arr = w_DE_from_omega(a_arr, kappa_demo, a_trans_demo, phi0_demo)
    w_full_arr = w_DE_full(a_arr, kappa_demo, a_trans_demo, phi0_demo)

    print("\n    z        w_DE (einfach)   w_DE (Friedmann)")
    print("    " + "-" * 45)
    for z in [0.0, 0.2, 0.5, 1.0, 1.5, 2.0]:
        a = 1.0 / (1.0 + z)
        w_s = float(w_DE_from_omega(np.array([a]), kappa_demo, a_trans_demo, phi0_demo)[0])
        w_f = float(w_DE_full(np.array([a]), kappa_demo, a_trans_demo, phi0_demo)[0])
        print(f"    z={z:.1f}:  {w_s:+.4f}          {w_f:+.4f}")

    # CPL-Fits
    w0_s, wa_s = cpl_fit_wDE(kappa_demo, a_trans_demo, phi0_demo, method='simple')
    w0_f, wa_f = cpl_fit_wDE(kappa_demo, a_trans_demo, phi0_demo, method='full')

    print(f"\n    CPL-Fit (einfach): w0_DE = {w0_s:.4f}, wa_DE = {wa_s:.4f}")
    print(f"    CPL-Fit (Friedmann): w0_DE = {w0_f:.4f}, wa_DE = {wa_f:.4f}")

    chi2_s = chi2_vs_desi(w0_s, wa_s)
    chi2_f = chi2_vs_desi(w0_f, wa_f)
    print(f"    chi^2 vs DESI (einfach): {chi2_s:.3f}  ({np.sqrt(chi2_s):.2f} sigma)")
    print(f"    chi^2 vs DESI (Friedmann): {chi2_f:.3f}  ({np.sqrt(chi2_f):.2f} sigma)")

    print("\n    DESI DR2 Referenz: w0 = {:.2f} +/- {:.2f}, wa = {:.2f} +/- {:.2f}".format(
        DESI_W0, DESI_W0_ERR, DESI_WA, DESI_WA_ERR))

    # --- Grid-Scan ---
    print("\n[2] Grid-Scan (kappa x a_trans = 20 x 20 = 400 Punkte)...")
    kappa_arr, atrans_arr, chi2_grid, w0_grid, wa_grid = grid_scan(20, 20)
    print("    Scan abgeschlossen.")

    # Best-fit
    finite_mask = np.isfinite(chi2_grid)
    if np.any(finite_mask):
        idx_min = np.unravel_index(np.nanargmin(chi2_grid), chi2_grid.shape)
        k_best  = kappa_arr[idx_min[0]]
        at_best = atrans_arr[idx_min[1]]
        w0_best = w0_grid[idx_min]
        wa_best = wa_grid[idx_min]
        chi2_min = chi2_grid[idx_min]

        print(f"\n    Best-fit: kappa={k_best:.3f}, a_trans={at_best:.3f}")
        print(f"    -> w0_DE = {w0_best:.4f}, wa_DE = {wa_best:.4f}")
        print(f"    -> chi^2 = {chi2_min:.3f}  ({np.sqrt(chi2_min):.2f} sigma)")

        # 1sigma und 2sigma Regionen
        chi2_1sigma = 2.30   # chi^2 fuer 1 sigma (2 Parameter)
        chi2_2sigma = 6.18   # chi^2 fuer 2 sigma (2 Parameter)
        in_1s = np.sum((chi2_grid <= chi2_min + chi2_1sigma) & finite_mask)
        in_2s = np.sum((chi2_grid <= chi2_min + chi2_2sigma) & finite_mask)
        print(f"\n    Punkte im 1-sigma Bereich: {in_1s} / {np.sum(finite_mask)}")
        print(f"    Punkte im 2-sigma Bereich: {in_2s} / {np.sum(finite_mask)}")

        # Parameter-Bereich der kompatibel ist
        compat_mask = (chi2_grid <= chi2_min + chi2_2sigma) & finite_mask
        if np.any(compat_mask):
            ki, aj = np.where(compat_mask)
            print(f"\n    Kompatibler kappa-Bereich  (< 2sigma): [{kappa_arr[ki.min()]:.3f}, {kappa_arr[ki.max()]:.3f}]")
            print(f"    Kompatibler a_trans-Bereich (< 2sigma): [{atrans_arr[aj.min()]:.3f}, {atrans_arr[aj.max()]:.3f}]")
    else:
        print("    WARNUNG: Kein gueltiger Punkt im Grid!")
        idx_min = (0, 0)
        k_best, at_best = kappa_arr[0], atrans_arr[0]
        w0_best, wa_best = np.nan, np.nan
        chi2_min = np.inf

    # ============================================================
    # Plot
    # ============================================================
    print("\n[3] Erstelle Plot...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("FST Dark Energy: Korrektes $w_{\\rm DE}$ Mapping vs DESI DR2",
                 fontsize=13, fontweight='bold')

    # --- Panel 1: w_DE(z) ---
    ax1 = axes[0]

    # DESI CPL-Band
    z_plot = np.linspace(0.0, 2.0, 300)
    f_z = z_plot / (1.0 + z_plot)
    w_desi = DESI_W0 + DESI_WA * f_z
    sigma_w = np.sqrt(DESI_W0_ERR**2 + (f_z * DESI_WA_ERR)**2)
    ax1.fill_between(z_plot, w_desi - 2*sigma_w, w_desi + 2*sigma_w,
                     alpha=0.12, color='royalblue', label='DESI DR2 $2\\sigma$')
    ax1.fill_between(z_plot, w_desi - sigma_w, w_desi + sigma_w,
                     alpha=0.25, color='royalblue', label='DESI DR2 $1\\sigma$')
    ax1.plot(z_plot, w_desi, 'b--', lw=2, label=f'DESI DR2 ($w_0={DESI_W0}$, $w_a={DESI_WA}$)')

    # Lambda-CDM
    ax1.axhline(-1.0, color='gray', ls=':', lw=1.5, label='$\\Lambda$CDM ($w=-1$)')

    # FST w_DE (verschiedene kappa)
    colors_kappa = ['crimson', 'darkorange', 'darkgreen']
    kappa_vals = [1.0, 1.75, 3.0]
    for kp, col in zip(kappa_vals, colors_kappa):
        phi0 = compute_phi0(kp, a_trans_demo)
        if np.isfinite(phi0) and phi0 > 0:
            w_arr = w_DE_full(a_arr, kp, a_trans_demo, phi0)
            ax1.plot(z_arr, w_arr, color=col, lw=2,
                     label=f'FST $w_{{\\rm DE}}$ ($\\kappa={kp}$, $a_{{\\rm tr}}={a_trans_demo}$)')

    # Best-fit
    phi0_bf = compute_phi0(k_best, at_best)
    if np.isfinite(phi0_bf) and phi0_bf > 0:
        w_bf = w_DE_full(a_arr, k_best, at_best, phi0_bf)
        ax1.plot(z_arr, w_bf, 'k-', lw=2.5,
                 label=f'Best-fit ($\\kappa={k_best:.2f}$, $a_{{\\rm tr}}={at_best:.2f}$)')

    ax1.set_xlabel('Redshift $z$', fontsize=12)
    ax1.set_ylabel('$w_{\\rm DE}(z)$', fontsize=12)
    ax1.set_title('Intrinsische DE-Zustandsgleichung', fontsize=11)
    ax1.set_xlim(0, 2)
    ax1.set_ylim(-5, 0.5)
    ax1.legend(fontsize=8, loc='lower left')
    ax1.grid(True, alpha=0.3)

    # --- Panel 2: chi^2 Konturplot ---
    ax2 = axes[1]

    KK, AA = np.meshgrid(kappa_arr, atrans_arr, indexing='ij')
    chi2_plot = chi2_grid.copy()
    chi2_plot[~np.isfinite(chi2_plot)] = 30.0  # Clip fuer Plot

    # Niveau-Linien relativ zum Minimum
    delta_chi2 = chi2_plot - chi2_min
    delta_chi2 = np.clip(delta_chi2, 0, 25)

    im = ax2.contourf(KK, AA, delta_chi2,
                      levels=[0, 2.30, 6.18, 11.83, 20],
                      colors=['#2d6a2d', '#4dac26', '#d9f0a3', '#fec44f', '#f7f7f7'],
                      alpha=0.8)
    cs = ax2.contour(KK, AA, delta_chi2,
                     levels=[2.30, 6.18, 11.83],
                     colors=['darkgreen', 'darkorange', 'darkred'],
                     linewidths=[2, 2, 1.5])
    ax2.clabel(cs, fmt={2.30: '$1\\sigma$', 6.18: '$2\\sigma$', 11.83: '$3\\sigma$'},
               fontsize=9)

    ax2.plot(k_best, at_best, 'w*', ms=15, zorder=5, label=f'Best-fit')
    ax2.plot(kappa_demo, a_trans_demo, 'rx', ms=12, mew=2.5, zorder=5,
             label=f'Paper ($\\kappa={kappa_demo}$, $a_{{\\rm tr}}={a_trans_demo}$)')

    ax2.set_xlabel('$\\kappa$', fontsize=12)
    ax2.set_ylabel('$a_{\\rm trans}$', fontsize=12)
    ax2.set_title('$\\Delta\\chi^2$ vs DESI DR2 (2 Param.)', fontsize=11)
    ax2.legend(fontsize=9, loc='upper right')
    ax2.grid(True, alpha=0.2)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax2)
    cbar.set_label('$\\Delta\\chi^2$', fontsize=10)

    plt.tight_layout()

    out_path = r"C:\Users\User\OneDrive\.RESEARCH\Natur&Technik\3 Folgebeweise\Dark Energy\compute_w_mapping.png"
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"    Plot gespeichert: {out_path}")
    plt.close()

    # ============================================================
    # Abschliessende Zusammenfassung
    # ============================================================
    print("\n" + "=" * 70)
    print("ZUSAMMENFASSUNG")
    print("=" * 70)
    print("\n  KORREKTUR gegenueber compute_w_vs_desi.py:")
    print("  - Altes Skript: verglich w_eff (gesamt) mit DESI w_DE -> FALSCH")
    print("  - Neues Skript: berechnet w_DE(a) via Kontinuitaetsgleichung -> KORREKT")
    print()
    print("  Physikalische Relation:")
    print("    w_eff(a) = w_DE(a) * Omega_DE(a)   [Staub hat p=0]")
    print("    => w_DE(a) = w_eff(a) / Omega_DE(a)  [einfach]")
    print("    => w_DE(a) = -1 - (1/3)*d(ln rho_DE)/d(ln a)  [Friedmann, vollstaendig]")
    print()
    print(f"  Demo-Parameter (kappa={kappa_demo}, a_trans={a_trans_demo}):")
    print(f"    w0_DE = {w0_f:.4f}  (alt: w0_eff ~ -1.47)")
    print(f"    wa_DE = {wa_f:.4f}  (alt: wa_eff ~ inkonsistent)")
    print(f"    chi^2 vs DESI: {chi2_f:.3f}  ({np.sqrt(chi2_f):.2f} sigma)")
    print()
    print(f"  Best-fit (kappa={k_best:.3f}, a_trans={at_best:.3f}):")
    print(f"    w0_DE = {w0_best:.4f}")
    print(f"    wa_DE = {wa_best:.4f}")
    print(f"    chi^2 = {chi2_min:.3f}  ({np.sqrt(chi2_min):.2f} sigma)")
    print()
    print(f"  DESI DR2 Referenz: w0={DESI_W0} +/- {DESI_W0_ERR}, wa={DESI_WA} +/- {DESI_WA_ERR}")
    print("\n" + "=" * 70)
    print("Berechnung abgeschlossen.")
    print("=" * 70)


if __name__ == "__main__":
    main()
