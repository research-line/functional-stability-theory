#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
II-CALC-01: 3-Spezies-Autokatalyse-Modell
==========================================

Reaktionsschema:
    A + B -> 2B   (Rate k1)
    B + C -> 2C   (Rate k2)
    C -> A        (Rate k3)

Totalmasse konserviert: x_A + x_B + x_C = 1 (Simplex).

Berechnet:
    - Fixpunkt x* via ODE-Integration
    - Nash-Gleichgewicht-Pruefung
    - Stabilitaet (Jacobian, volle Eigenwerte, Tangentialraum-Eigenwerte)
    - Gibbs-Funktion G(x)
    - Entropieproduktion S_dot
    - Korrelationsanalyse rho(J) vs. S_dot ueber k1-Variation

Beweisstück fuer FST-II (Functional Stability Theory).

Autor: Lukas Geiger (ORCID: 0009-0005-7296-1534)
Datum: 2026-03-23
"""

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import pearsonr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# 1. Modell-Parameter
# ---------------------------------------------------------------------------

DEFAULT_K1 = 1.0
DEFAULT_K2 = 0.8
DEFAULT_K3 = 0.5

MU0 = np.array([0.0, 0.2, 0.4])  # chemische Standard-Potentiale

# Verhaeltnis k_fwd / k_rev fuer Entropieproduktion (weit vom GG)
K_FWD_REV_RATIO = 100.0

# Anfangsbedingung auf dem Simplex
X0 = np.array([0.6, 0.3, 0.1])

# ODE-Integrationsparameter
T_MAX = 500.0
T_EVAL_POINTS = 5000
CONVERGENCE_RTOL = 1e-10
CONVERGENCE_ATOL = 1e-12

# Korrelationsanalyse
K1_RANGE = np.linspace(0.1, 5.0, 20)


# ---------------------------------------------------------------------------
# 2. ODEs (Replikator-Form)
# ---------------------------------------------------------------------------

def ode_rhs(t, x, k1, k2, k3):
    """Rechte Seite des ODE-Systems auf dem Simplex."""
    x_a, x_b, x_c = np.clip(x, 1e-15, None)  # numerischer Schutz

    dx_a = -k1 * x_a * x_b + k3 * x_c
    dx_b = k1 * x_a * x_b - k2 * x_b * x_c
    dx_c = k2 * x_b * x_c - k3 * x_c

    return [dx_a, dx_b, dx_c]


def normalize_simplex(x):
    """Projiziere x auf den Simplex (x_i >= 0, sum = 1)."""
    x = np.clip(x, 1e-15, None)
    return x / np.sum(x)


# ---------------------------------------------------------------------------
# 3. Fixpunkt via ODE-Integration
# ---------------------------------------------------------------------------

def find_fixed_point(k1, k2, k3, x0=None, t_max=T_MAX):
    """Integriere ODE bis zur Konvergenz und gib den Fixpunkt zurueck."""
    if x0 is None:
        x0 = X0.copy()

    def rhs(t, x):
        return ode_rhs(t, x, k1, k2, k3)

    sol = solve_ivp(
        rhs,
        [0, t_max],
        x0,
        method="RK45",
        rtol=CONVERGENCE_RTOL,
        atol=CONVERGENCE_ATOL,
        dense_output=True,
        max_step=0.5,
    )

    if not sol.success:
        raise RuntimeError(f"ODE-Integration fehlgeschlagen: {sol.message}")

    x_final = normalize_simplex(sol.y[:, -1])
    return x_final, sol


# ---------------------------------------------------------------------------
# 4. Gibbs-Funktion
# ---------------------------------------------------------------------------

def gibbs(x, mu0=MU0):
    """Berechne G = sum_i x_i * (mu0_i + ln(x_i))."""
    x = np.clip(x, 1e-15, None)
    return np.sum(x * (mu0 + np.log(x)))


# ---------------------------------------------------------------------------
# 5. Nash-Gleichgewicht-Pruefung
# ---------------------------------------------------------------------------

def check_nash(x_star, k1, k2, k3, delta=1e-4):
    """
    Pruefe ob x* ein Nash-GG ist.

    Fuer jeden Spieler i: teste ob unilaterale Abweichung
    (x_i + delta) die Auszahlung (= -Gibbs bei Konvergenz) verbessert.

    In diesem kontinuierlichen Replikator-System entspricht ein
    stabiler Fixpunkt einem Nash-GG.
    """
    g_star = gibbs(x_star)
    is_nash = True
    deviations = []

    for i in range(3):
        # Unilaterale Abweichung: Spieler i erhoehe seinen Anteil
        x_pert = x_star.copy()
        x_pert[i] += delta
        x_pert = normalize_simplex(x_pert)

        # Integriere ab der perturbierten Bedingung
        x_conv, _ = find_fixed_point(k1, k2, k3, x0=x_pert, t_max=200.0)
        g_conv = gibbs(x_conv)

        # Pruefen ob Abweichung profitabel ist (niedrigeres G = besser)
        improved = g_conv < g_star - 1e-8
        deviations.append({
            "player": ["A", "B", "C"][i],
            "G_star": float(g_star),
            "G_deviated": float(g_conv),
            "improved": bool(improved),
        })
        if improved:
            is_nash = False

    return is_nash, deviations


# ---------------------------------------------------------------------------
# 6. Jacobian und Stabilitaetsanalyse
# ---------------------------------------------------------------------------

def jacobian_numerical(x_star, k1, k2, k3, eps=1e-7):
    """Numerischer Jacobian df/dx am Fixpunkt (3x3)."""
    n = len(x_star)
    J = np.zeros((n, n))
    f0 = np.array(ode_rhs(0, x_star, k1, k2, k3))

    for j in range(n):
        x_pert = x_star.copy()
        x_pert[j] += eps
        f_pert = np.array(ode_rhs(0, x_pert, k1, k2, k3))
        J[:, j] = (f_pert - f0) / eps

    return J


def jacobian_analytical(x_star, k1, k2, k3):
    """Analytischer Jacobian df/dx am Fixpunkt."""
    x_a, x_b, x_c = x_star

    J = np.array([
        [-k1 * x_b,       -k1 * x_a,       k3        ],
        [ k1 * x_b,        k1 * x_a - k2 * x_c, -k2 * x_b],
        [ 0,                k2 * x_c,        k2 * x_b - k3],
    ])
    return J


def tangent_restriction_simplex(J):
    """
    Restrict a 3x3 Jacobian to the tangent space sum_i dx_i = 0.

    The full concentration dynamics conserves total mass, so a neutral
    conservation mode is expected in the unreduced Jacobian. Physical
    asymptotic stability is tested on the simplex tangent space.
    """
    basis = np.array([
        [1.0, 0.0],
        [0.0, 1.0],
        [-1.0, -1.0],
    ])
    gram_inv = np.linalg.inv(basis.T @ basis)
    return gram_inv @ basis.T @ J @ basis


def stability_analysis(x_star, k1, k2, k3):
    """Eigenwerte, Spektralradius, Tangentialraum-Stabilitaet."""
    J = jacobian_analytical(x_star, k1, k2, k3)
    eigenvalues = np.linalg.eigvals(J)
    J_tangent = tangent_restriction_simplex(J)
    tangent_eigenvalues = np.linalg.eigvals(J_tangent)
    spectral_radius = np.max(np.abs(eigenvalues))
    tangent_spectral_radius = np.max(np.abs(tangent_eigenvalues))
    max_real = np.max(np.real(eigenvalues))
    tangent_max_real = np.max(np.real(tangent_eigenvalues))
    is_full_stable = max_real < -1e-10
    is_tangent_stable = tangent_max_real < -1e-10

    return {
        "jacobian": J,
        "eigenvalues": eigenvalues,
        "tangent_jacobian": J_tangent,
        "tangent_eigenvalues": tangent_eigenvalues,
        "spectral_radius": float(spectral_radius),
        "tangent_spectral_radius": float(tangent_spectral_radius),
        "max_real_eigenvalue": float(max_real),
        "tangent_max_real_eigenvalue": float(tangent_max_real),
        "is_full_stable": bool(is_full_stable),
        "is_tangent_stable": bool(is_tangent_stable),
        "is_stable": bool(is_tangent_stable),
    }


# ---------------------------------------------------------------------------
# 7. Entropieproduktion
# ---------------------------------------------------------------------------

def entropy_production(x, k1, k2, k3, ratio=K_FWD_REV_RATIO):
    """
    Entropieproduktion S_dot = sum_r J_r * ln(k_fwd / k_rev).

    Reaktionen:
        R1: A + B -> 2B   J1 = k1 * x_A * x_B
        R2: B + C -> 2C   J2 = k2 * x_B * x_C
        R3: C -> A         J3 = k3 * x_C

    Vereinfachung (weit vom GG): S_dot = sum J_r * ln(ratio)
    """
    x_a, x_b, x_c = np.clip(x, 1e-15, None)

    J1 = k1 * x_a * x_b
    J2 = k2 * x_b * x_c
    J3 = k3 * x_c

    ln_ratio = np.log(ratio)
    s_dot = (J1 + J2 + J3) * ln_ratio

    return float(s_dot), {"J1": float(J1), "J2": float(J2), "J3": float(J3)}


# ---------------------------------------------------------------------------
# 8. Korrelationsanalyse: rho(J) vs. S_dot ueber k1-Variation
# ---------------------------------------------------------------------------

def correlation_scan(k1_values, k2, k3):
    """Variiere k1, berechne Fixpunkt, rho(J), S_dot fuer jeden Wert."""
    results = []

    for k1 in k1_values:
        try:
            x_star, _ = find_fixed_point(k1, k2, k3)
            stab = stability_analysis(x_star, k1, k2, k3)
            s_dot, fluxes = entropy_production(x_star, k1, k2, k3)

            results.append({
                "k1": float(k1),
                "x_star": x_star.tolist(),
                "spectral_radius": stab["spectral_radius"],
                "tangent_spectral_radius": stab["tangent_spectral_radius"],
                "max_real_eigenvalue": stab["max_real_eigenvalue"],
                "tangent_max_real_eigenvalue": stab["tangent_max_real_eigenvalue"],
                "is_full_stable": stab["is_full_stable"],
                "is_tangent_stable": stab["is_tangent_stable"],
                "is_stable": stab["is_stable"],
                "S_dot": s_dot,
                "fluxes": fluxes,
            })
        except Exception as e:
            print(f"  [WARNUNG] k1={k1:.2f}: {e}")
            results.append({
                "k1": float(k1),
                "error": str(e),
            })

    # Pearson-Korrelation (nur erfolgreiche Punkte)
    valid = [r for r in results if "spectral_radius" in r]
    rho_vals = np.array([r["spectral_radius"] for r in valid])
    s_dot_vals = np.array([r["S_dot"] for r in valid])

    if len(valid) >= 3:
        corr, p_value = pearsonr(rho_vals, s_dot_vals)
    else:
        corr, p_value = float("nan"), float("nan")

    return results, float(corr), float(p_value), rho_vals, s_dot_vals


# ---------------------------------------------------------------------------
# 9. Plotting
# ---------------------------------------------------------------------------

def create_plot(rho_vals, s_dot_vals, k1_values_valid, corr, p_value, output_path):
    """Erstelle rho(J) vs. S_dot Scatterplot."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    scatter = ax.scatter(
        s_dot_vals, rho_vals,
        c=k1_values_valid, cmap="viridis", s=60, edgecolors="black", linewidth=0.5
    )
    cbar = plt.colorbar(scatter, ax=ax, label=r"$k_1$")

    ax.set_xlabel(r"Entropy Production $\dot{S}$", fontsize=13)
    ax.set_ylabel(r"Spectral Radius $\rho(J)$", fontsize=13)
    ax.set_title(
        f"II-CALC-01: 3-Species Autocatalysis\n"
        f"Pearson r = {corr:.4f}, p = {p_value:.2e}",
        fontsize=14,
    )
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot gespeichert: {output_path}")


# ---------------------------------------------------------------------------
# 10. Hauptprogramm
# ---------------------------------------------------------------------------

def main():
    script_dir = Path(__file__).parent
    output_dir = script_dir / "output"
    output_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("  II-CALC-01: 3-Spezies-Autokatalyse-Modell")
    print("=" * 70)

    k1, k2, k3 = DEFAULT_K1, DEFAULT_K2, DEFAULT_K3
    print(f"\n  Parameter: k1={k1}, k2={k2}, k3={k3}")
    print(f"  Anfangsbedingung: x0 = {X0.tolist()}")

    # --- Fixpunkt ---
    print("\n--- Fixpunkt-Berechnung ---")
    x_star, sol = find_fixed_point(k1, k2, k3)
    print(f"  x* = [{x_star[0]:.8f}, {x_star[1]:.8f}, {x_star[2]:.8f}]")
    print(f"  Summe: {np.sum(x_star):.12f}")

    # --- Gibbs ---
    g_star = gibbs(x_star)
    print(f"  G(x*) = {g_star:.8f}")

    # --- Residuum am Fixpunkt ---
    residuum = np.array(ode_rhs(0, x_star, k1, k2, k3))
    print(f"  |f(x*)| = {np.linalg.norm(residuum):.2e}")

    # --- Nash-GG ---
    print("\n--- Nash-Gleichgewicht-Pruefung ---")
    is_nash, deviations = check_nash(x_star, k1, k2, k3)
    for d in deviations:
        status = "NEIN (profitabel)" if d["improved"] else "JA (nicht profitabel)"
        print(f"  Spieler {d['player']}: Abweichung {status}")
        print(f"    G* = {d['G_star']:.8f}, G_dev = {d['G_deviated']:.8f}")
    print(f"  => Nash-GG: {'JA' if is_nash else 'NEIN'}")

    # --- Stabilitaet ---
    print("\n--- Stabilitaetsanalyse ---")
    stab = stability_analysis(x_star, k1, k2, k3)
    print(f"  Jacobian:\n{stab['jacobian']}")
    print(f"  Eigenwerte: {stab['eigenvalues']}")
    print(f"  Tangentialraum-Jacobian:\n{stab['tangent_jacobian']}")
    print(f"  Tangentialraum-Eigenwerte: {stab['tangent_eigenvalues']}")
    print(f"  Spektralradius rho(J) = {stab['spectral_radius']:.8f}")
    print(f"  Max. Realteil = {stab['max_real_eigenvalue']:.8f}")
    print(f"  Tangentiale Spektralabszisse = {stab['tangent_max_real_eigenvalue']:.8f}")
    print(f"  => Stabil im vollen Raum: {'JA' if stab['is_full_stable'] else 'NEIN (neutraler Erhaltungsmodus)'}")
    print(f"  => Stabil im Simplex-Tangentialraum: {'JA' if stab['is_tangent_stable'] else 'NEIN'}")

    # --- Entropieproduktion ---
    print("\n--- Entropieproduktion ---")
    s_dot, fluxes = entropy_production(x_star, k1, k2, k3)
    print(f"  J1 (A+B->2B) = {fluxes['J1']:.8f}")
    print(f"  J2 (B+C->2C) = {fluxes['J2']:.8f}")
    print(f"  J3 (C->A)    = {fluxes['J3']:.8f}")
    print(f"  S_dot = {s_dot:.8f}")

    # --- Korrelationsanalyse ---
    print("\n--- Korrelationsanalyse: rho(J) vs. S_dot ---")
    print(f"  k1-Bereich: {K1_RANGE[0]:.1f} bis {K1_RANGE[-1]:.1f} ({len(K1_RANGE)} Punkte)")
    scan_results, corr, p_value, rho_vals, s_dot_vals = correlation_scan(K1_RANGE, k2, k3)

    valid_results = [r for r in scan_results if "spectral_radius" in r]
    print(f"  Erfolgreiche Punkte: {len(valid_results)}/{len(K1_RANGE)}")
    print(f"  Pearson r = {corr:.6f}")
    print(f"  p-Wert   = {p_value:.2e}")

    if abs(corr) > 0.7:
        corr_desc = "STARK"
    elif abs(corr) > 0.4:
        corr_desc = "MITTEL"
    else:
        corr_desc = "SCHWACH"
    print(f"  => Korrelation: {corr_desc} ({'positiv' if corr > 0 else 'negativ'})")

    # --- Plot ---
    print("\n--- Plot ---")
    k1_valid = np.array([r["k1"] for r in valid_results])
    plot_path = str(output_dir / "ii_calc_01_rho_vs_sdot.png")
    create_plot(rho_vals, s_dot_vals, k1_valid, corr, p_value, plot_path)

    # --- JSON-Export ---
    print("\n--- JSON-Export ---")
    json_output = {
        "model": "II-CALC-01: 3-Species Autocatalysis",
        "reactions": [
            "A + B -> 2B (k1)",
            "B + C -> 2C (k2)",
            "C -> A (k3)",
        ],
        "default_parameters": {"k1": k1, "k2": k2, "k3": k3},
        "mu0": MU0.tolist(),
        "initial_condition": X0.tolist(),
        "fixed_point": {
            "x_star": x_star.tolist(),
            "sum_check": float(np.sum(x_star)),
            "residuum_norm": float(np.linalg.norm(residuum)),
            "G_star": float(g_star),
        },
        "nash_equilibrium": {
            "is_nash": is_nash,
            "deviations": deviations,
        },
        "stability": {
            "eigenvalues_real": np.real(stab["eigenvalues"]).tolist(),
            "eigenvalues_imag": np.imag(stab["eigenvalues"]).tolist(),
            "tangent_eigenvalues_real": np.real(stab["tangent_eigenvalues"]).tolist(),
            "tangent_eigenvalues_imag": np.imag(stab["tangent_eigenvalues"]).tolist(),
            "spectral_radius": stab["spectral_radius"],
            "tangent_spectral_radius": stab["tangent_spectral_radius"],
            "max_real_eigenvalue": stab["max_real_eigenvalue"],
            "tangent_max_real_eigenvalue": stab["tangent_max_real_eigenvalue"],
            "is_full_stable": stab["is_full_stable"],
            "is_tangent_stable": stab["is_tangent_stable"],
            "is_stable": stab["is_stable"],
            "stability_interpretation": (
                "Full-space stability is false because total concentration is conserved; "
                "physical stability is tested on the simplex tangent space."
            ),
        },
        "entropy_production": {
            "S_dot": s_dot,
            "fluxes": fluxes,
            "k_fwd_rev_ratio": K_FWD_REV_RATIO,
        },
        "correlation_analysis": {
            "k1_range": K1_RANGE.tolist(),
            "pearson_r": corr,
            "p_value": p_value,
            "n_valid_points": len(valid_results),
            "scan_results": scan_results,
        },
    }

    json_path = str(output_dir / "ii_calc_01_results.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_output, f, indent=2, ensure_ascii=False)
    print(f"  JSON gespeichert: {json_path}")

    # --- Zusammenfassung ---
    print("\n" + "=" * 70)
    print("  ZUSAMMENFASSUNG")
    print("=" * 70)
    print(f"  Fixpunkt x*:        [{x_star[0]:.6f}, {x_star[1]:.6f}, {x_star[2]:.6f}]")
    print(f"  Nash-GG:            {'JA' if is_nash else 'NEIN'}")
    print(f"  Stabil (Tangential): {'JA' if stab['is_tangent_stable'] else 'NEIN'}")
    print(f"  Spektralradius:     {stab['spectral_radius']:.6f}")
    print(f"  Spektralabszisse T: {stab['tangent_max_real_eigenvalue']:.6f}")
    print(f"  Entropieproduktion: {s_dot:.6f}")
    print(f"  Korrelation r:      {corr:.6f} ({corr_desc})")
    print(f"  p-Wert:             {p_value:.2e}")
    print("=" * 70)

    return json_output


if __name__ == "__main__":
    main()
