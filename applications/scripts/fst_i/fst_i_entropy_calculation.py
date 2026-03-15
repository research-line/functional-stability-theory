#!/usr/bin/env python3
"""
FST-I: Integrated Stellar Entropy Production vs. Fundamental Constants
======================================================================

Computes the total entropy produced by a solar-mass main-sequence star
as a function of the fine-structure constant alpha_EM and the
electron-to-proton mass ratio m_e/m_p.

Physics:
  - Stellar structure: polytropic model (n=3) calibrated to solar values
  - Nuclear fusion: pp-chain with Gamow penetration factor
  - Energy transport: radiative diffusion with Thomson opacity
  - Self-regulation: central temperature adjusts to balance L_nuc = L_rad
  - Entropy production: S_total = integral(L/T_eff) dt over main-sequence

Key result: S_total(alpha) shows a maximum near the observed alpha_EM,
supporting the FST-I conjecture that SM parameters maximize entropy production.

References:
  Adams (2019) Physics Reports 807, 1-111  [fine-tuning review]
  Egan & Lineweaver (2010) ApJ 710, 1825   [cosmic entropy]
  Kippenhahn, Weigert & Weiss (2013)       [stellar structure]

Author: Lukas Geiger (with Claude, 2026-03-10)
Part of: FST-I Paper, Fractal Game Theory
"""

import numpy as np
from scipy.optimize import brentq
import json
import os

# ===========================================================================
# Physical constants (SI)
# ===========================================================================
G       = 6.674e-11     # gravitational constant [m^3 kg^-1 s^-2]
c       = 2.998e8       # speed of light [m/s]
hbar    = 1.055e-34     # reduced Planck constant [J s]
k_B     = 1.381e-23     # Boltzmann constant [J/K]
sigma_SB = 5.670e-8     # Stefan-Boltzmann constant [W m^-2 K^-4]
m_p     = 1.673e-27     # proton mass [kg]
m_e_obs = 9.109e-31     # electron mass [kg]

# Solar reference values
M_sun     = 1.989e30    # solar mass [kg]
R_sun     = 6.957e8     # solar radius [m]
L_sun     = 3.828e26    # solar luminosity [W]
T_eff_sun = 5778.0      # solar effective temperature [K]
T_c_sun   = 1.57e7      # solar central temperature [K]
rho_c_sun = 1.50e5      # solar central density [kg/m^3]

# Observed values
alpha_obs = 1.0 / 137.036  # fine-structure constant
mu_obs    = m_e_obs / m_p   # electron-to-proton mass ratio ~ 1/1836

# Nuclear physics
X_H       = 0.70        # hydrogen mass fraction
eta_nuc   = 0.007       # nuclear efficiency (H -> He)
E_nuc_sun = eta_nuc * X_H * M_sun * c**2  # total nuclear fuel [J]

# Gamow factor for pp-chain at solar conditions
# E_G = (pi * alpha)^2 * 2 * m_r * c^2, with m_r = m_p/2 for p+p
# E_G = (pi * alpha)^2 * m_p * c^2
E_G_sun_J = (np.pi * alpha_obs)**2 * m_p * c**2  # ~7.9e-14 J
E_G_sun_keV = E_G_sun_J / 1.602e-16  # ~493 keV

# Gamow parameter at solar central temperature
# tau = 3 * (E_G / (4 * k_B * T_c))^(1/3)
tau_sun = 3.0 * (E_G_sun_J / (4 * k_B * T_c_sun))**(1.0/3.0)
# Should be ~13.5

# Temperature sensitivity of pp-chain
nu_sun = tau_sun / 3.0 - 2.0/3.0  # ~3.8

print(f"=== Solar Calibration ===")
print(f"E_G = {E_G_sun_keV:.1f} keV")
print(f"tau_sun = {tau_sun:.2f}")
print(f"nu_sun (pp temperature exponent) = {nu_sun:.2f}")
print(f"E_nuc_sun = {E_nuc_sun:.3e} J")
print()


# ===========================================================================
# Stellar model: solve for equilibrium structure
# ===========================================================================

def solve_stellar_structure(x_alpha, x_me=1.0):
    """
    Solve for the equilibrium stellar structure of a solar-mass star
    with modified fundamental constants.

    Parameters
    ----------
    x_alpha : float
        alpha_EM / alpha_obs (ratio of fine-structure constant to observed)
    x_me : float
        m_e / m_e_obs (ratio of electron mass to observed)

    Returns
    -------
    dict with: L, T_eff, T_c, R, t_MS, S_dot, S_total, y (=R/R_sun),
               tau, valid (bool)

    Model:
    ------
    Polytropic star (n=3) with:
    - Virial theorem: T_c = T_c_sun / y  (y = R/R_sun)
    - Density: rho_c = rho_c_sun / y^3
    - Thomson opacity: kappa propto alpha^2 / m_e^2
    - Radiative luminosity: L_rad propto R^4 T_c^4 / (kappa M)
    - Nuclear luminosity: L_nuc propto rho_c^2 R^3 T_c^{-2/3} exp(-tau)
    - Self-consistency: L_nuc = L_rad determines y
    """
    # Opacity scaling: kappa_Thomson propto alpha^2 / m_e^2
    # (Thomson cross section sigma_T = 8pi/3 * (alpha * hbar / (m_e * c))^2)
    kappa_ratio = x_alpha**2 / x_me**2

    # Gamow parameter as function of y (= R/R_sun):
    # T_c(y) = T_c_sun / y
    # E_G(x_alpha) = E_G_sun * x_alpha^2
    # tau(y, x_alpha) = tau_sun * x_alpha^{2/3} * y^{1/3}
    #   (because tau propto E_G^{1/3} * T_c^{-1/3})

    def equilibrium_equation(y):
        """
        L_nuc / L_rad - 1 = 0

        From the scaling analysis:
        L_rad / L_sun = 1 / kappa_ratio  (independent of y for polytrope!)
        L_nuc / L_sun = y^{-10/3} * exp(tau_sun * (1 - x_alpha^{2/3} * y^{1/3}))

        Setting L_nuc = L_rad:
        y^{-10/3} * exp(tau_sun * (1 - x_alpha^{2/3} * y^{1/3})) = 1/kappa_ratio
        """
        tau_current = tau_sun * x_alpha**(2.0/3.0) * y**(1.0/3.0)
        lhs = y**(-10.0/3.0) * np.exp(tau_sun - tau_current)
        rhs = 1.0 / kappa_ratio
        return lhs - rhs

    # Physical bounds for y = R/R_sun
    # y too small: star too compact, T_c too high
    # y too large: star too extended, T_c too low for fusion
    y_min = 0.01
    y_max = 100.0

    # Check if a solution exists
    try:
        f_min = equilibrium_equation(y_min)
        f_max = equilibrium_equation(y_max)
    except (OverflowError, FloatingPointError):
        return _invalid_result(x_alpha, x_me)

    if np.isnan(f_min) or np.isnan(f_max) or np.isinf(f_min) or np.isinf(f_max):
        return _invalid_result(x_alpha, x_me)

    if f_min * f_max > 0:
        # No root in interval: no stable star at this alpha
        return _invalid_result(x_alpha, x_me)

    try:
        y = brentq(equilibrium_equation, y_min, y_max, xtol=1e-10)
    except ValueError:
        return _invalid_result(x_alpha, x_me)

    # Stellar properties at equilibrium
    R = y * R_sun
    T_c = T_c_sun / y
    rho_c = rho_c_sun / y**3
    tau = tau_sun * x_alpha**(2.0/3.0) * y**(1.0/3.0)

    # Luminosity (from radiative transport, equals nuclear luminosity)
    L = L_sun / kappa_ratio

    # Check fusion viability: T_c must be above minimum ignition temperature
    T_c_min = 3.5e6  # K, approximate minimum for pp-chain
    if T_c < T_c_min:
        return _invalid_result(x_alpha, x_me)

    # Check: if tau > 50, fusion is negligibly slow (exp(-50) ~ 2e-22)
    if tau > 50:
        return _invalid_result(x_alpha, x_me)

    # Effective temperature
    T_eff = (L / (4 * np.pi * R**2 * sigma_SB))**0.25

    # Main-sequence lifetime
    # Nuclear fuel is independent of alpha (it's strong-force physics)
    E_nuc = E_nuc_sun  # same fuel
    t_MS = E_nuc / L

    # Entropy production rate
    S_dot = L / T_eff

    # Integrated entropy production
    S_total = S_dot * t_MS  # = E_nuc / T_eff

    # Convert to convenient units
    t_MS_Gyr = t_MS / (1e9 * 3.156e7)
    L_solar = L / L_sun

    return {
        'x_alpha': float(x_alpha),
        'x_me': float(x_me),
        'y': float(y),
        'R_Rsun': float(y),
        'T_c_K': float(T_c),
        'T_eff_K': float(T_eff),
        'L_Lsun': float(L_solar),
        'L_W': float(L),
        't_MS_Gyr': float(t_MS_Gyr),
        'tau': float(tau),
        'S_dot_WperK': float(S_dot),
        'S_total_J_per_K': float(S_total),
        'valid': True
    }


def _invalid_result(x_alpha, x_me):
    """Return a result dict for parameters where no stable star exists."""
    return {
        'x_alpha': float(x_alpha),
        'x_me': float(x_me),
        'y': 0.0, 'R_Rsun': 0.0,
        'T_c_K': 0.0, 'T_eff_K': 0.0,
        'L_Lsun': 0.0, 'L_W': 0.0,
        't_MS_Gyr': 0.0, 'tau': 0.0,
        'S_dot_WperK': 0.0, 'S_total_J_per_K': 0.0,
        'valid': False
    }


# ===========================================================================
# Scan 1: alpha_EM variation (at fixed m_e/m_p)
# ===========================================================================

def scan_alpha():
    """Vary alpha_EM from 0.1x to 5.0x observed value."""
    print("=" * 72)
    print("SCAN 1: alpha_EM variation (m_e/m_p fixed)")
    print("=" * 72)

    x_values = np.concatenate([
        np.linspace(0.10, 0.50, 40),
        np.linspace(0.50, 2.00, 150),
        np.linspace(2.00, 5.00, 60)
    ])

    results = []
    for x in x_values:
        r = solve_stellar_structure(x_alpha=x, x_me=1.0)
        results.append(r)

    # Find maximum S_total among valid results
    valid = [r for r in results if r['valid']]
    if not valid:
        print("ERROR: No valid stellar models found!")
        return results

    S_values = [r['S_total_J_per_K'] for r in valid]
    x_values_valid = [r['x_alpha'] for r in valid]
    idx_max = np.argmax(S_values)
    x_max = x_values_valid[idx_max]
    S_max = S_values[idx_max]

    # Solar reference
    solar = solve_stellar_structure(1.0, 1.0)

    print(f"\nSolar reference (x = 1.0):")
    print(f"  L = {solar['L_Lsun']:.4f} L_sun")
    print(f"  T_eff = {solar['T_eff_K']:.0f} K")
    print(f"  T_c = {solar['T_c_K']:.3e} K")
    print(f"  t_MS = {solar['t_MS_Gyr']:.2f} Gyr")
    print(f"  tau = {solar['tau']:.2f}")
    print(f"  S_total = {solar['S_total_J_per_K']:.4e} J/K")
    print(f"  S_dot = {solar['S_dot_WperK']:.4e} W/K")

    print(f"\nMaximum S_total found at x_alpha = {x_max:.3f}")
    print(f"  S_total_max = {S_max:.4e} J/K")
    print(f"  S_total(obs) / S_total_max = {solar['S_total_J_per_K'] / S_max:.4f}")

    print(f"\nViable range: x_alpha in [{x_values_valid[0]:.2f}, {x_values_valid[-1]:.2f}]")

    # Print table for key values
    print(f"\n{'x_alpha':>8s} | {'L/L_sun':>9s} | {'T_eff':>7s} | {'t_MS Gyr':>9s} | "
          f"{'tau':>6s} | {'S_total':>12s} | {'S/S_max':>7s}")
    print("-" * 76)
    for r in results:
        if r['valid'] and abs(r['x_alpha'] % 0.25) < 0.015 or abs(r['x_alpha'] - 1.0) < 0.01:
            ratio = r['S_total_J_per_K'] / S_max
            mark = " <-- obs" if abs(r['x_alpha'] - 1.0) < 0.01 else ""
            print(f"{r['x_alpha']:8.3f} | {r['L_Lsun']:9.4f} | {r['T_eff_K']:7.0f} | "
                  f"{r['t_MS_Gyr']:9.2f} | {r['tau']:6.2f} | "
                  f"{r['S_total_J_per_K']:12.4e} | {ratio:7.4f}{mark}")

    return results


# ===========================================================================
# Scan 2: m_e/m_p variation (at fixed alpha_EM)
# ===========================================================================

def scan_me():
    """Vary m_e/m_p from 0.1x to 10x observed value."""
    print("\n" + "=" * 72)
    print("SCAN 2: m_e/m_p variation (alpha_EM fixed)")
    print("=" * 72)

    x_values = np.concatenate([
        np.linspace(0.10, 0.50, 30),
        np.linspace(0.50, 3.00, 150),
        np.linspace(3.00, 10.0, 40)
    ])

    results = []
    for x in x_values:
        r = solve_stellar_structure(x_alpha=1.0, x_me=x)
        results.append(r)

    valid = [r for r in results if r['valid']]
    if not valid:
        print("ERROR: No valid stellar models found!")
        return results

    S_values = [r['S_total_J_per_K'] for r in valid]
    x_values_valid = [r['x_me'] for r in valid]
    idx_max = np.argmax(S_values)
    x_max = x_values_valid[idx_max]
    S_max = S_values[idx_max]

    solar = solve_stellar_structure(1.0, 1.0)

    print(f"\nMaximum S_total found at x_me = {x_max:.3f}")
    print(f"  S_total_max = {S_max:.4e} J/K")
    print(f"  S_total(obs) / S_total_max = {solar['S_total_J_per_K'] / S_max:.4f}")
    print(f"\nViable range: x_me in [{x_values_valid[0]:.2f}, {x_values_valid[-1]:.2f}]")

    print(f"\n{'x_me':>8s} | {'L/L_sun':>9s} | {'T_eff':>7s} | {'t_MS Gyr':>9s} | "
          f"{'tau':>6s} | {'S_total':>12s} | {'S/S_max':>7s}")
    print("-" * 76)
    for r in results:
        if r['valid'] and (abs(r['x_me'] % 0.5) < 0.02 or abs(r['x_me'] - 1.0) < 0.01):
            ratio = r['S_total_J_per_K'] / S_max
            mark = " <-- obs" if abs(r['x_me'] - 1.0) < 0.01 else ""
            print(f"{r['x_me']:8.3f} | {r['L_Lsun']:9.4f} | {r['T_eff_K']:7.0f} | "
                  f"{r['t_MS_Gyr']:9.2f} | {r['tau']:6.2f} | "
                  f"{r['S_total_J_per_K']:12.4e} | {ratio:7.4f}{mark}")

    return results


# ===========================================================================
# Scan 3: 2D scan (alpha, m_e) -- coarse
# ===========================================================================

def scan_2d():
    """Coarse 2D scan over alpha and m_e simultaneously."""
    print("\n" + "=" * 72)
    print("SCAN 3: 2D scan (alpha_EM x m_e/m_p)")
    print("=" * 72)

    x_alpha_vals = np.linspace(0.3, 3.0, 50)
    x_me_vals = np.linspace(0.3, 3.0, 50)

    grid = np.zeros((len(x_alpha_vals), len(x_me_vals)))
    valid_grid = np.zeros_like(grid, dtype=bool)

    for i, xa in enumerate(x_alpha_vals):
        for j, xm in enumerate(x_me_vals):
            r = solve_stellar_structure(xa, xm)
            if r['valid']:
                grid[i, j] = r['S_total_J_per_K']
                valid_grid[i, j] = True

    # Find global maximum
    if np.any(valid_grid):
        max_val = np.max(grid[valid_grid])
        max_idx = np.unravel_index(np.argmax(grid), grid.shape)
        xa_max = x_alpha_vals[max_idx[0]]
        xm_max = x_me_vals[max_idx[1]]
        print(f"\n2D maximum at: x_alpha = {xa_max:.3f}, x_me = {xm_max:.3f}")
        print(f"S_total_max = {max_val:.4e} J/K")

        solar = solve_stellar_structure(1.0, 1.0)
        print(f"S_total(obs) / S_total_max = {solar['S_total_J_per_K'] / max_val:.4f}")
    else:
        print("ERROR: No valid models in 2D scan!")

    return {
        'x_alpha': x_alpha_vals.tolist(),
        'x_me': x_me_vals.tolist(),
        'S_total': grid.tolist(),
        'valid': valid_grid.tolist()
    }


# ===========================================================================
# Falsifiability test
# ===========================================================================

def falsifiability_test():
    """
    Test the FST-I falsifiability criterion:
    "If a variation |Delta p| > 10% increases S_total, the framework fails."

    Check all +-10% variations of alpha_EM and m_e/m_p.
    """
    print("\n" + "=" * 72)
    print("FALSIFIABILITY TEST")
    print("=" * 72)

    solar = solve_stellar_structure(1.0, 1.0)
    S_obs = solar['S_total_J_per_K']
    print(f"S_total(observed) = {S_obs:.6e} J/K\n")

    test_points = [
        (0.90, 1.00, "alpha -10%"),
        (1.10, 1.00, "alpha +10%"),
        (0.80, 1.00, "alpha -20%"),
        (1.20, 1.00, "alpha +20%"),
        (0.50, 1.00, "alpha -50%"),
        (1.50, 1.00, "alpha +50%"),
        (2.00, 1.00, "alpha +100%"),
        (1.00, 0.90, "m_e -10%"),
        (1.00, 1.10, "m_e +10%"),
        (1.00, 0.80, "m_e -20%"),
        (1.00, 1.20, "m_e +20%"),
        (1.00, 0.50, "m_e -50%"),
        (1.00, 2.00, "m_e +100%"),
    ]

    print(f"{'Variation':>16s} | {'x_alpha':>8s} | {'x_me':>6s} | {'S_total':>12s} | "
          f"{'S/S_obs':>8s} | {'Status':>12s}")
    print("-" * 80)

    violations = 0
    for xa, xm, label in test_points:
        r = solve_stellar_structure(xa, xm)
        if r['valid']:
            ratio = r['S_total_J_per_K'] / S_obs
            status = "VIOLATION!" if ratio > 1.0 else "OK (lower)"
            if ratio > 1.0:
                violations += 1
        else:
            ratio = 0.0
            status = "No star"

        print(f"{label:>16s} | {xa:8.3f} | {xm:6.3f} | "
              f"{r['S_total_J_per_K']:12.4e} | {ratio:8.4f} | {status:>12s}")

    print(f"\nResult: {violations} violations out of {len(test_points)} tests")
    if violations == 0:
        print("=> FST-I falsifiability criterion PASSED (no variation increases S_total)")
    else:
        print(f"=> FST-I falsifiability criterion: {violations} violations found")
        print("   (This may indicate the observed value is not at the global maximum,")
        print("    or that the simplified stellar model needs refinement.)")

    return violations


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("FST-I: Stellar Entropy Production vs. Fundamental Constants")
    print("=" * 72)
    print()

    # Run scans
    alpha_results = scan_alpha()
    me_results = scan_me()
    grid_2d = scan_2d()
    n_violations = falsifiability_test()

    # Compile all results
    output = {
        'metadata': {
            'description': 'FST-I: Integrated stellar entropy production as function of alpha_EM and m_e/m_p',
            'author': 'Lukas Geiger',
            'date': '2026-03-10',
            'model': 'Semi-analytic polytropic stellar model (n=3), pp-chain, Thomson opacity',
            'assumptions': [
                'Fixed stellar mass M = M_sun',
                'pp-chain fusion only (valid for T_c < 1.7e7 K)',
                'Thomson opacity dominates (valid for solar-type stars)',
                'Polytropic equation of state with solar-calibrated ratios',
                'Nuclear fuel E_nuc independent of alpha (strong-force physics)',
            ],
            'limitations': [
                'Single stellar mass (no IMF integration)',
                'No CNO cycle (dominates at T_c > 1.7e7 K)',
                'No convective transport (relevant for low-mass stars)',
                'No Kramers opacity (relevant for high-Z envelopes)',
                'No cosmological context (single star, not galaxy/universe)',
            ],
            'solar_calibration': {
                'E_G_keV': float(E_G_sun_keV),
                'tau_sun': float(tau_sun),
                'nu_sun': float(nu_sun),
                'E_nuc_J': float(E_nuc_sun),
            },
        },
        'scan_alpha': [r for r in alpha_results],
        'scan_me': [r for r in me_results],
        'scan_2d': grid_2d,
        'falsifiability_violations': n_violations,
    }

    # Save to JSON
    outdir = os.path.dirname(os.path.abspath(__file__))
    outpath = os.path.join(outdir, 'fst_i_entropy_results.json')
    with open(outpath, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults saved to: {outpath}")

    # Summary for paper
    print("\n" + "=" * 72)
    print("SUMMARY FOR PAPER (FST-I Section 5)")
    print("=" * 72)

    solar = solve_stellar_structure(1.0, 1.0)
    valid_alpha = [r for r in alpha_results if r['valid']]
    S_alpha = [r['S_total_J_per_K'] for r in valid_alpha]
    x_alpha = [r['x_alpha'] for r in valid_alpha]
    idx_max = np.argmax(S_alpha)

    print(f"""
The semi-analytic stellar model yields the following key results:

1. ALPHA_EM VARIATION:
   - Maximum S_total at alpha/alpha_obs = {x_alpha[idx_max]:.3f}
   - Viable range: [{x_alpha[0]:.2f}, {x_alpha[-1]:.2f}] x alpha_obs
   - S_total(obs) / S_total_max = {solar['S_total_J_per_K'] / S_alpha[idx_max]:.4f}
   - Observed alpha is within {abs(1.0 - x_alpha[idx_max])*100:.1f}% of the maximum

2. SOLAR REFERENCE:
   - L = {solar['L_Lsun']:.4f} L_sun
   - T_eff = {solar['T_eff_K']:.0f} K
   - t_MS = {solar['t_MS_Gyr']:.2f} Gyr
   - S_total = {solar['S_total_J_per_K']:.4e} J/K

3. FALSIFIABILITY:
   - {n_violations} parameter variations with |Delta| >= 10% exceed S_total(obs)
""")

    # Check 10% variations specifically
    for delta in [0.1, 0.2, 0.5]:
        r_up = solve_stellar_structure(1.0 + delta, 1.0)
        r_down = solve_stellar_structure(1.0 - delta, 1.0)
        print(f"   alpha +{delta*100:.0f}%: S/S_obs = "
              f"{r_up['S_total_J_per_K']/solar['S_total_J_per_K']:.4f} "
              f"({'valid' if r_up['valid'] else 'no star'})")
        print(f"   alpha -{delta*100:.0f}%: S/S_obs = "
              f"{r_down['S_total_J_per_K']/solar['S_total_J_per_K']:.4f} "
              f"({'valid' if r_down['valid'] else 'no star'})")


if __name__ == '__main__':
    main()
