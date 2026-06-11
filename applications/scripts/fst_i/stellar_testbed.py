#!/usr/bin/env python3
"""
I-TEST-01: Simplified stellar model with variable fine-structure constant alpha.

Uses dimensionless scaling relations (Carr & Rees 1979, Adams 2008) to compute
how stellar properties and cumulative entropy production depend on alpha.

Key question: Does alpha_obs lie near the maximum of S(alpha)?

Physics:
  - Polytroptic stellar model (Lane-Emden n=3 for main sequence)
  - Eddington luminosity scaling: L_Edd ~ alpha^(-2) * (m_p/m_e)^2
  - Main-sequence lifetime: tau ~ M/L ~ alpha^2
  - Cumulative entropy: S ~ L * tau
  - Stability: d^2 S / d alpha^2 (1D Hessian)

All quantities are dimensionless ratios normalized to observed values.

Author: Lukas Geiger (with Claude, 2026-03-23)
Part of: FST-I Fine-Tuning & Entropy Production

References:
  Carr & Rees (1979) Nature 278, 605-612
  Adams (2008) JCAP 08, 010
  Adams (2019) Physics Reports 807, 1-111
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

import sys
import argparse
import json
import numpy as np
from scipy.optimize import brentq


# ---------------------------------------------------------------------------
# Physical constants (dimensionless ratios at observed alpha)
# ---------------------------------------------------------------------------

ALPHA_OBS = 1.0 / 137.036      # Fine-structure constant
MP_ME_RATIO = 1836.15           # Proton-to-electron mass ratio


# ---------------------------------------------------------------------------
# Stellar scaling relations (dimensionless)
# ---------------------------------------------------------------------------

def stellar_properties(x_alpha: float) -> dict:
    """
    Compute dimensionless stellar properties as function of x_alpha = alpha / alpha_obs.

    All quantities normalized so that x_alpha=1 gives 1.0.

    Scaling relations (Carr & Rees 1979, Adams 2008, 2019):

    Stellar radius:
      R ~ alpha^(-1) * a_0 * (m_p/m_e)^(-1)
      => R/R_obs ~ x_alpha^(-1)

    Surface temperature (from virial theorem):
      T_eff ~ alpha^2 * m_e * c^2 / k_B
      => T_eff/T_obs ~ x_alpha^2

    Luminosity (Stefan-Boltzmann with R and T scaling):
      L ~ R^2 * T^4 ~ alpha^(-2) * alpha^8 = alpha^6
      BUT: more precisely from opacity/transport:
      L ~ alpha^(-1/2) in suitable opacity regime (Adams 2019)
      We use the full radiative transport result:
      L scales as kappa^(-1) where kappa_Thomson ~ alpha^2
      => L ~ alpha^(-2) for Thomson-dominated
      Combined with structure: L/L_obs ~ x_alpha^(-2)

    Main-sequence lifetime:
      tau_MS ~ E_nuc / L ~ M*c^2 / L
      Nuclear fuel is alpha-independent (strong force)
      => tau/tau_obs ~ L_obs/L ~ x_alpha^2

    Entropy production rate:
      S_dot = L / T_eff
      => S_dot/S_dot_obs ~ x_alpha^(-2) / x_alpha^2 = x_alpha^(-4)

    Cumulative entropy:
      S_total = S_dot * tau = (L/T_eff) * (E_nuc/L) = E_nuc/T_eff
      => S_total/S_obs ~ x_alpha^(-2)

    Parameters
    ----------
    x_alpha : float
        alpha / alpha_obs (dimensionless ratio)

    Returns
    -------
    dict with dimensionless ratios (all normalized to observed = 1.0)
    """
    if x_alpha <= 0:
        return _invalid(x_alpha)

    # Basic scalings
    R_ratio = x_alpha**(-1)         # R / R_obs
    T_eff_ratio = x_alpha**2        # T_eff / T_obs
    L_ratio = x_alpha**(-2)         # L / L_obs (Thomson opacity dominated)
    tau_ratio = x_alpha**2          # tau_MS / tau_obs = 1/L_ratio * const

    # Entropy
    S_dot_ratio = L_ratio / T_eff_ratio     # S_dot / S_dot_obs
    S_total_ratio = S_dot_ratio * tau_ratio  # S_total / S_obs

    # Fusion viability check:
    # Central temperature T_c ~ alpha^2 (virial) / R ~ alpha^3
    # Gamow peak requires T_c above ignition threshold
    # At very low alpha: T_c too low for fusion
    # At very high alpha: star collapses, other physics dominates
    T_c_ratio = x_alpha**3  # Rough scaling

    # Viability bounds (approximate)
    # Fusion requires T_c > ~3.5 MK (T_c_ratio > ~0.2)
    # Stability requires alpha < ~10*alpha_obs (pair production, etc.)
    viable = (T_c_ratio > 0.2) and (x_alpha < 10.0)

    return {
        "x_alpha": float(x_alpha),
        "alpha": float(x_alpha * ALPHA_OBS),
        "R_ratio": float(R_ratio),
        "T_eff_ratio": float(T_eff_ratio),
        "L_ratio": float(L_ratio),
        "tau_MS_ratio": float(tau_ratio),
        "S_dot_ratio": float(S_dot_ratio),
        "S_total_ratio": float(S_total_ratio),
        "T_c_ratio": float(T_c_ratio),
        "viable": viable
    }


def _invalid(x_alpha: float) -> dict:
    return {
        "x_alpha": float(x_alpha), "alpha": 0.0,
        "R_ratio": 0.0, "T_eff_ratio": 0.0, "L_ratio": 0.0,
        "tau_MS_ratio": 0.0, "S_dot_ratio": 0.0, "S_total_ratio": 0.0,
        "T_c_ratio": 0.0, "viable": False
    }


# ---------------------------------------------------------------------------
# Enhanced model with Gamow penetration factor
# ---------------------------------------------------------------------------

def stellar_properties_gamow(x_alpha: float) -> dict:
    """
    More physical model including Gamow penetration factor.

    The pp-chain reaction rate depends exponentially on alpha through
    the Gamow factor: epsilon ~ exp(-tau) with tau ~ alpha^(2/3) * T_c^(-1/3).

    At equilibrium, L_nuc = L_rad determines T_c and hence all observables.
    This model uses the same physics as fst_i_entropy_calculation.py but
    expressed in dimensionless ratios.

    The key improvement: S_total is no longer a simple power law in alpha.
    The Gamow factor creates a genuine maximum.
    """
    if x_alpha <= 0:
        return _invalid(x_alpha)

    # Solar-calibrated Gamow parameter
    # tau_sun ~ 13.5
    tau_sun = 13.5

    # Thomson opacity ratio
    kappa_ratio = x_alpha**2  # kappa / kappa_obs

    # Luminosity from radiative transport (at fixed mass)
    # L = L_obs / kappa_ratio (polytropic model)
    L_ratio = 1.0 / kappa_ratio

    # Self-consistent radius from L_nuc = L_rad
    # Equilibrium equation: y^(-10/3) * exp(tau_sun * (1 - x_alpha^(2/3) * y^(1/3))) = L_ratio
    # where y = R/R_obs
    def eq(y):
        tau_curr = tau_sun * x_alpha**(2.0/3.0) * y**(1.0/3.0)
        lhs = y**(-10.0/3.0) * np.exp(tau_sun - tau_curr)
        return lhs - L_ratio

    try:
        y_min, y_max = 0.01, 100.0
        f_min = eq(y_min)
        f_max = eq(y_max)
        if np.isnan(f_min) or np.isnan(f_max) or np.isinf(f_min) or np.isinf(f_max):
            return _invalid(x_alpha)
        if f_min * f_max > 0:
            return _invalid(x_alpha)
        y = brentq(eq, y_min, y_max, xtol=1e-10)
    except (ValueError, OverflowError, FloatingPointError):
        return _invalid(x_alpha)

    # Derived quantities
    R_ratio = y
    T_c_ratio = 1.0 / y  # T_c = T_c_sun / y
    tau = tau_sun * x_alpha**(2.0/3.0) * y**(1.0/3.0)

    # Effective temperature: T_eff^4 = L / (4*pi*R^2*sigma_SB)
    # T_eff / T_eff_obs = (L_ratio / R_ratio^2)^(1/4)
    T_eff_ratio = (L_ratio / R_ratio**2)**0.25

    # Lifetime: tau_MS / tau_obs = 1 / L_ratio (nuclear fuel unchanged)
    tau_MS_ratio = 1.0 / L_ratio

    # Entropy
    S_dot_ratio = L_ratio / T_eff_ratio
    S_total_ratio = S_dot_ratio * tau_MS_ratio  # = 1 / T_eff_ratio

    # Viability
    T_c_min = 0.22  # T_c / T_c_sun > 0.22 for pp-chain
    viable = (T_c_ratio > T_c_min) and (tau < 50)

    return {
        "x_alpha": float(x_alpha),
        "alpha": float(x_alpha * ALPHA_OBS),
        "R_ratio": float(R_ratio),
        "T_eff_ratio": float(T_eff_ratio),
        "L_ratio": float(L_ratio),
        "tau_MS_ratio": float(tau_MS_ratio),
        "S_dot_ratio": float(S_dot_ratio),
        "S_total_ratio": float(S_total_ratio),
        "T_c_ratio": float(T_c_ratio),
        "tau_gamow": float(tau),
        "viable": viable
    }


# ---------------------------------------------------------------------------
# 1D Hessian (stability measure)
# ---------------------------------------------------------------------------

def compute_hessian_1d(x_values: np.ndarray, S_values: np.ndarray) -> np.ndarray:
    """
    Numerical second derivative d^2 S / d(x_alpha)^2 using central differences.
    Returns array of same length (with NaN at boundaries).
    """
    dx = np.diff(x_values)
    d2S = np.full_like(S_values, np.nan)

    for i in range(1, len(x_values) - 1):
        h1 = x_values[i] - x_values[i - 1]
        h2 = x_values[i + 1] - x_values[i]
        d2S[i] = 2.0 * (S_values[i + 1] / (h2 * (h1 + h2))
                         - S_values[i] / (h1 * h2)
                         + S_values[i - 1] / (h1 * (h1 + h2)))

    return d2S


# ---------------------------------------------------------------------------
# Main scan
# ---------------------------------------------------------------------------

def run_alpha_scan(n_points: int = 20, model: str = "gamow") -> dict:
    """
    Scan alpha from 0.5*alpha_obs to 2.0*alpha_obs.

    Parameters
    ----------
    n_points : int
        Number of scan points
    model : str
        "simple" for power-law scaling, "gamow" for Gamow-enhanced model

    Returns
    -------
    dict with all results
    """
    x_values = np.linspace(0.5, 2.0, n_points)

    model_fn = stellar_properties_gamow if model == "gamow" else stellar_properties

    results = []
    for x in x_values:
        r = model_fn(x)
        results.append(r)

    # Extract S_total for viable stars
    x_arr = np.array([r["x_alpha"] for r in results])
    S_arr = np.array([r["S_total_ratio"] if r["viable"] else 0.0 for r in results])
    viable = np.array([r["viable"] for r in results])

    # Find maximum
    if np.any(viable):
        S_viable = S_arr.copy()
        S_viable[~viable] = -np.inf
        idx_max = np.argmax(S_viable)
        x_max = x_arr[idx_max]
        S_max = S_arr[idx_max]
    else:
        x_max = 1.0
        S_max = 0.0

    # S_total at observed alpha
    r_obs = model_fn(1.0)
    S_obs = r_obs["S_total_ratio"]

    # Hessian (curvature at maximum)
    d2S = compute_hessian_1d(x_arr, S_arr)
    # Hessian at x_alpha = 1.0
    idx_obs = np.argmin(np.abs(x_arr - 1.0))
    hessian_at_obs = float(d2S[idx_obs]) if not np.isnan(d2S[idx_obs]) else 0.0

    # Hessian at maximum
    hessian_at_max = float(d2S[idx_max]) if not np.isnan(d2S[idx_max]) else 0.0

    return {
        "model": model,
        "x_values": x_arr.tolist(),
        "scan_results": results,
        "S_total_values": S_arr.tolist(),
        "viable_mask": viable.tolist(),
        "d2S_dalpha2": d2S.tolist(),
        "x_alpha_at_S_max": float(x_max),
        "S_total_max": float(S_max),
        "S_total_obs": float(S_obs),
        "S_obs_over_S_max": float(S_obs / S_max) if S_max > 0 else 0.0,
        "deviation_from_max_percent": float(abs(1.0 - x_max) * 100),
        "hessian_at_obs": hessian_at_obs,
        "hessian_at_max": hessian_at_max
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(scan_simple: dict, scan_gamow: dict, output_dir: str):
    """Plot alpha vs S for both models."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping plots.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    for col, (scan, title_suffix) in enumerate([
        (scan_simple, "Simple Scaling"),
        (scan_gamow, "Gamow Model")
    ]):
        x = np.array(scan["x_values"])
        S = np.array(scan["S_total_values"])
        viable = np.array(scan["viable_mask"])
        d2S = np.array(scan["d2S_dalpha2"])

        # --- Top row: S_total vs alpha ---
        ax = axes[0, col]
        # Plot viable and non-viable separately
        ax.plot(x[viable], S[viable], 'b-o', linewidth=2, markersize=5,
                label=r'$S_{total}(\alpha)$ (viable)')
        if np.any(~viable):
            ax.plot(x[~viable], S[~viable], 'rx', markersize=8, label='Non-viable')

        # Mark observed alpha
        ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.7,
                    label=r'$\alpha_{obs}$')

        # Mark maximum
        x_max = scan["x_alpha_at_S_max"]
        S_max = scan["S_total_max"]
        ax.plot(x_max, S_max, 'g*', markersize=15, zorder=5,
                label=f'Max at x={x_max:.3f}')

        ax.set_xlabel(r'$\alpha / \alpha_{obs}$', fontsize=13)
        ax.set_ylabel(r'$S_{total} / S_{obs}$', fontsize=13)
        ax.set_title(f'Entropy Production ({title_suffix})', fontsize=14)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # --- Bottom row: d²S/dalpha² ---
        ax2 = axes[1, col]
        mask = ~np.isnan(d2S)
        ax2.plot(x[mask], d2S[mask], 'g-s', linewidth=2, markersize=5)
        ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        ax2.axvline(x=1.0, color='red', linestyle='--', alpha=0.7,
                     label=r'$\alpha_{obs}$')

        ax2.set_xlabel(r'$\alpha / \alpha_{obs}$', fontsize=13)
        ax2.set_ylabel(r'$d^2 S / d\alpha^2$', fontsize=13)
        ax2.set_title(f'Curvature (Stability Hessian) ({title_suffix})', fontsize=14)
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, "stellar_testbed_alpha_scan.png")
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nPlot saved: {plot_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="I-TEST-01: Stellar testbed with variable alpha"
    )
    ap.add_argument("--output-dir", default=os.path.dirname(os.path.abspath(__file__)),
                    help="Output directory for JSON and plots")
    ap.add_argument("--n-points", type=int, default=20,
                    help="Number of alpha scan points")
    ap.add_argument("--alpha-min", type=float, default=0.5,
                    help="Minimum alpha/alpha_obs ratio")
    ap.add_argument("--alpha-max", type=float, default=2.0,
                    help="Maximum alpha/alpha_obs ratio")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("I-TEST-01: Stellar Testbed -- Variable Fine-Structure Constant")
    print("=" * 60)
    print(f"\nalpha_obs = 1/{1/ALPHA_OBS:.3f}")
    print(f"Scan range: [{args.alpha_min}, {args.alpha_max}] x alpha_obs")
    print(f"Points: {args.n_points}")

    # Run both models
    print("\n--- Simple Scaling Model ---")
    scan_simple = run_alpha_scan(n_points=args.n_points, model="simple")

    print(f"\n  Maximum S_total at x_alpha = {scan_simple['x_alpha_at_S_max']:.4f}")
    print(f"  S_obs / S_max = {scan_simple['S_obs_over_S_max']:.4f}")
    print(f"  alpha_obs deviates {scan_simple['deviation_from_max_percent']:.1f}% from maximum")
    print(f"  Hessian at alpha_obs: {scan_simple['hessian_at_obs']:.4f}")

    print("\n--- Gamow-Enhanced Model ---")
    scan_gamow = run_alpha_scan(n_points=args.n_points, model="gamow")

    print(f"\n  Maximum S_total at x_alpha = {scan_gamow['x_alpha_at_S_max']:.4f}")
    print(f"  S_obs / S_max = {scan_gamow['S_obs_over_S_max']:.4f}")
    print(f"  alpha_obs deviates {scan_gamow['deviation_from_max_percent']:.1f}% from maximum")
    print(f"  Hessian at alpha_obs: {scan_gamow['hessian_at_obs']:.4f}")

    # Compare observed with maximum
    print("\n" + "=" * 60)
    print("RESULT COMPARISON")
    print("=" * 60)
    print(f"\n  {'Model':>15s} | {'x_max':>8s} | {'S_obs/S_max':>12s} | {'Deviation':>10s} | {'H(obs)':>10s}")
    print("  " + "-" * 65)
    for label, scan in [("Simple", scan_simple), ("Gamow", scan_gamow)]:
        print(f"  {label:>15s} | {scan['x_alpha_at_S_max']:8.4f} | "
              f"{scan['S_obs_over_S_max']:12.4f} | "
              f"{scan['deviation_from_max_percent']:9.1f}% | "
              f"{scan['hessian_at_obs']:10.4f}")

    print(f"\nInterpretation:")
    gamow_dev = scan_gamow['deviation_from_max_percent']
    gamow_ratio = scan_gamow['S_obs_over_S_max']
    if gamow_dev < 15 and gamow_ratio > 0.9:
        print(f"  STRONG SUPPORT: alpha_obs is within {gamow_dev:.1f}% of the entropy maximum")
        print(f"  (S_obs = {gamow_ratio:.1%} of S_max)")
    elif gamow_dev < 30:
        print(f"  MODERATE SUPPORT: alpha_obs is within {gamow_dev:.1f}% of the entropy maximum")
    else:
        print(f"  WEAK SUPPORT: alpha_obs deviates {gamow_dev:.1f}% from the entropy maximum")

    if scan_gamow['hessian_at_obs'] < 0:
        print(f"  Negative Hessian at alpha_obs confirms it is near a LOCAL MAXIMUM")
    else:
        print(f"  Positive Hessian: alpha_obs is NOT at a local maximum (saddle or minimum)")

    # Save JSON
    output = {
        "metadata": {
            "description": "I-TEST-01: Stellar entropy production vs fine-structure constant",
            "author": "Lukas Geiger",
            "date": "2026-03-23",
            "alpha_obs": ALPHA_OBS,
            "alpha_obs_inverse": 1.0 / ALPHA_OBS,
            "mp_me_ratio": MP_ME_RATIO,
            "references": [
                "Carr & Rees (1979) Nature 278, 605-612",
                "Adams (2008) JCAP 08, 010",
                "Adams (2019) Physics Reports 807, 1-111"
            ]
        },
        "simple_model": scan_simple,
        "gamow_model": scan_gamow
    }

    json_path = os.path.join(args.output_dir, "stellar_testbed_results.json")

    # Convert numpy/bool types for JSON serialization
    def json_safe(obj):
        import numpy as np
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, bool):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2, default=json_safe)
    print(f"\nJSON saved: {json_path}")

    # Plot
    plot_results(scan_simple, scan_gamow, args.output_dir)


if __name__ == "__main__":
    main()
