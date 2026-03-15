#!/usr/bin/env python3
"""
FST-II: Replicator Dynamics and Entropy Production for Autocatalytic Networks
=============================================================================

Implements the quantitative test P3 from FST-II paper:
Three competing autocatalytic networks (Azoarcus-type ribozymes) modeled as
replicator dynamics. Tests whether the observed winner (cooperative network)
maximizes entropy production rate.

Physics:
  - Replicator equation: dx_i/dt = x_i * (f_i(x) - <f>(x))
  - Fitness from catalytic rate constants (Vaidya et al. 2012)
  - Entropy production: sigma = sum_i x_i * f_i * ln(f_i / <f>)
  - Nash equilibrium: interior ESS where no species can improve unilaterally

Key predictions:
  1. Cooperative network reaches Nash equilibrium (ESS)
  2. At ESS, entropy production rate is maximized among stable fixed points
  3. Selfish replicators are dynamically unstable (higher eigenvalues)

References:
  Vaidya et al. (2012) Nature 491, 72-77 [Azoarcus cooperative RNA]
  Hofbauer & Sigmund (1998) Evolutionary Games and Population Dynamics
  England (2013) J. Chem. Phys. 139, 121923 [dissipative adaptation]

Author: Lukas Geiger (with Claude, 2026-03-15)
Part of: FST-II Paper, Fractal Game Theory
"""

import json
import numpy as np
from scipy.integrate import solve_ivp
import os
import argparse

# ===========================================================================
# Autocatalytic network parameters (from Vaidya et al. 2012, Table 1)
# ===========================================================================

# Three competing network types:
# 1. Cooperative (C): Cross-catalytic, Azoarcus-type
# 2. Selfish (S): Self-catalytic replicator
# 3. Parasitic (P): Exploits catalysis without reciprocating

# Fitness matrix A[i,j] = catalytic rate of species j on species i
# Based on Vaidya et al. (2012) measurements:
# - Cooperative fragments: ~10x higher catalytic rate in network
# - Selfish: ~3x rate advantage initially but no cooperation bonus
# - Parasitic: benefits from others but contributes ~0.1x

FITNESS_MATRICES = {
    "azoarcus_3species": {
        "description": "3-species: Cooperative, Selfish, Parasitic",
        "names": ["Cooperative", "Selfish", "Parasitic"],
        "A": np.array([
            [1.0, 0.1, 0.05],   # Cooperative: good self, low from others
            [0.3, 3.0, 0.1],    # Selfish: strong self-catalysis
            [0.8, 0.5, 0.1],    # Parasitic: exploits cooperative
        ]),
        "cooperation_bonus": np.array([
            [5.0, 0.0, 0.0],   # C gets 5x bonus from C-C interaction
            [0.0, 0.0, 0.0],   # S gets no cooperation bonus
            [2.0, 0.0, 0.0],   # P gets 2x bonus from C (parasitism)
        ]),
    },
    "azoarcus_4species": {
        "description": "4-species: Cooperative A, Cooperative B, Selfish, Parasitic",
        "names": ["Coop-A", "Coop-B", "Selfish", "Parasitic"],
        "A": np.array([
            [0.5, 4.0, 0.1, 0.05],  # Coop-A: needs Coop-B
            [4.0, 0.5, 0.1, 0.05],  # Coop-B: needs Coop-A
            [0.3, 0.3, 3.0, 0.1],   # Selfish: self-catalytic
            [1.5, 1.5, 0.5, 0.1],   # Parasitic: exploits coops
        ]),
        "cooperation_bonus": np.array([
            [0.0, 3.0, 0.0, 0.0],
            [3.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0, 0.0],
        ]),
    },
}


def effective_fitness(x, A, coop_bonus):
    """Compute effective fitness vector f(x) including cooperation effects."""
    f_base = A @ x
    f_coop = coop_bonus @ x
    return f_base + f_coop * x


def replicator_rhs(t, x, A, coop_bonus):
    """Replicator equation dx_i/dt = x_i * (f_i - <f>)."""
    x = np.maximum(x, 1e-15)
    x = x / x.sum()
    f = effective_fitness(x, A, coop_bonus)
    f_mean = x @ f
    return x * (f - f_mean)


def entropy_production_rate(x, A, coop_bonus):
    """
    Entropy production rate sigma = sum_i x_i * f_i * ln(f_i / <f>).
    KL divergence of fitness from mean, weighted by population and fitness.
    """
    x = np.maximum(x, 1e-15)
    f = effective_fitness(x, A, coop_bonus)
    f = np.maximum(f, 1e-15)
    f_mean = max(x @ f, 1e-15)
    return float(np.sum(x * f * np.log(f / f_mean)))


def compute_jacobian(x, A, coop_bonus):
    """Jacobian of the replicator equation at state x."""
    n = len(x)
    f = effective_fitness(x, A, coop_bonus)
    f_mean = x @ f
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                dfi_dxi = A[i, i] + 2 * coop_bonus[i, i] * x[i]
                df_mean_dxi = f[i] + sum(x[k] * (A[k, i] + coop_bonus[k, i] * x[k]) for k in range(n))
                J[i, j] = (f[i] - f_mean) + x[i] * (dfi_dxi - df_mean_dxi)
            else:
                dfi_dxj = A[i, j] + coop_bonus[i, j] * x[i]
                df_mean_dxj = f[j] + sum(x[k] * (A[k, j] + coop_bonus[k, j] * x[k]) for k in range(n))
                J[i, j] = x[i] * (dfi_dxj - df_mean_dxj)
    return J


def find_fixed_points(A, coop_bonus, n_random=500):
    """Find fixed points of replicator dynamics by simulation from random starts."""
    n = A.shape[0]
    fixed_points = []

    for trial in range(n_random):
        x0 = np.random.dirichlet(np.ones(n))
        try:
            sol = solve_ivp(
                replicator_rhs, [0, 500], x0,
                args=(A, coop_bonus),
                method="RK45", rtol=1e-10, atol=1e-12,
                max_step=1.0
            )
            x_final = sol.y[:, -1]
            x_final = np.maximum(x_final, 0)
            x_final = x_final / x_final.sum()

            rhs = replicator_rhs(0, x_final, A, coop_bonus)
            if np.max(np.abs(rhs)) < 1e-8:
                is_new = True
                for fp in fixed_points:
                    if np.max(np.abs(np.array(fp["x"]) - x_final)) < 1e-4:
                        is_new = False
                        fp["count"] += 1
                        break
                if is_new:
                    fixed_points.append({"x": x_final.tolist(), "count": 1})
        except Exception:
            continue

    # Analyze each fixed point
    for fp in fixed_points:
        x = np.array(fp["x"])
        f = effective_fitness(x, A, coop_bonus)
        f_mean = float(x @ f)

        fp["sigma"] = entropy_production_rate(x, A, coop_bonus)
        fp["f_mean"] = f_mean
        fp["f"] = f.tolist()

        J = compute_jacobian(x, A, coop_bonus)
        eigs = np.linalg.eigvals(J)
        eigs_sorted = sorted(eigs.real)
        fp["eigenvalues"] = [float(e) for e in eigs_sorted]
        fp["max_eigenvalue"] = float(max(eigs.real))
        fp["stable"] = all(e.real < 1e-6 for e in eigs)
        fp["n_positive_eigs"] = sum(1 for e in eigs if e.real > 1e-6)
        fp["dominant"] = int(np.argmax(x))

    return fixed_points


def run_dynamics(A, coop_bonus, names, x0=None, t_max=200):
    """Run replicator dynamics and track entropy production."""
    n = A.shape[0]
    if x0 is None:
        x0 = np.ones(n) / n

    t_eval = np.linspace(0, t_max, 2000)
    sol = solve_ivp(
        replicator_rhs, [0, t_max], x0,
        args=(A, coop_bonus),
        method="RK45", t_eval=t_eval,
        rtol=1e-10, atol=1e-12,
        max_step=0.5
    )

    sigmas = []
    for i in range(len(sol.t)):
        x = sol.y[:, i]
        x = np.maximum(x, 1e-15)
        x = x / x.sum()
        sigmas.append(entropy_production_rate(x, A, coop_bonus))

    return {
        "t": sol.t.tolist(),
        "x": {names[i]: sol.y[i].tolist() for i in range(n)},
        "sigma": sigmas,
        "x_final": sol.y[:, -1].tolist(),
        "sigma_final": sigmas[-1],
    }


def main():
    parser = argparse.ArgumentParser(description="FST-II Replicator Entropy Test")
    parser.add_argument("--scenario", default="azoarcus_3species",
                        choices=list(FITNESS_MATRICES.keys()),
                        help="Which network scenario to simulate")
    parser.add_argument("--out", default=None, help="Output JSON file")
    parser.add_argument("--n-random", type=int, default=500,
                        help="Number of random starts for fixed point search")
    args = parser.parse_args()

    scenario = FITNESS_MATRICES[args.scenario]
    A = scenario["A"]
    coop = scenario["cooperation_bonus"]
    names = scenario["names"]
    n = len(names)

    print("FST-II Replicator Entropy Test: %s" % scenario["description"])
    print("Species: %s" % names)
    print("\nFitness matrix A:")
    for i in range(n):
        print("  %12s: %s" % (names[i], A[i]))
    print("\nCooperation bonus:")
    for i in range(n):
        print("  %12s: %s" % (names[i], coop[i]))

    # Phase 1: Find fixed points
    print("\n" + "=" * 60)
    print("PHASE 1: Finding fixed points (%d random starts)" % args.n_random)
    print("=" * 60)

    fps = find_fixed_points(A, coop, n_random=args.n_random)
    fps.sort(key=lambda fp: -fp["sigma"])

    print("\nFound %d distinct fixed points:" % len(fps))
    for i, fp in enumerate(fps):
        dom = names[fp["dominant"]]
        stab = "STABLE" if fp["stable"] else "UNSTABLE"
        print("\n  FP-%d (%s, found %dx):" % (i + 1, stab, fp["count"]))
        for j in range(n):
            print("    %12s: x=%.6f, f=%.4f" % (names[j], fp["x"][j], fp["f"][j]))
        print("    Entropy production: sigma = %.6f" % fp["sigma"])
        print("    Mean fitness: <f> = %.4f" % fp["f_mean"])
        print("    Dominant: %s" % dom)
        print("    Eigenvalues: %s" % [round(e, 4) for e in fp["eigenvalues"]])
        print("    Positive eigenvalues: %d" % fp["n_positive_eigs"])

    # Phase 2: Dynamics from uniform start
    print("\n" + "=" * 60)
    print("PHASE 2: Dynamics from uniform start")
    print("=" * 60)

    dyn = run_dynamics(A, coop, names)
    print("\nFinal state (t=%.0f):" % dyn["t"][-1])
    for name in names:
        print("  %12s: %.6f" % (name, dyn["x"][name][-1]))
    print("  Entropy production: %.6f" % dyn["sigma_final"])

    # Phase 3: P3 Test
    print("\n" + "=" * 60)
    print("PHASE 3: P3 TEST -- Does the ESS maximize entropy production?")
    print("=" * 60)

    stable_fps = [fp for fp in fps if fp["stable"]]
    unstable_fps = [fp for fp in fps if not fp["stable"]]

    if stable_fps:
        max_sigma_stable = max(fp["sigma"] for fp in stable_fps)
        best_stable = max(stable_fps, key=lambda fp: fp["sigma"])
        print("\nBest stable FP: %s (sigma = %.6f)" % (names[best_stable["dominant"]], best_stable["sigma"]))

        if unstable_fps:
            best_unstable = max(unstable_fps, key=lambda fp: fp["sigma"])
            print("Best unstable FP: %s (sigma = %.6f)" % (names[best_unstable["dominant"]], best_unstable["sigma"]))

        coop_fps = [fp for fp in stable_fps if fp["dominant"] == 0]
        if coop_fps:
            coop_sigma = max(fp["sigma"] for fp in coop_fps)
            print("\nP3 Result: Cooperative ESS sigma = %.6f" % coop_sigma)
            print("           Max stable sigma      = %.6f" % max_sigma_stable)
            if abs(coop_sigma - max_sigma_stable) < 1e-6:
                print("  => P3 CONFIRMED: Cooperative ESS maximizes entropy production among stable states!")
            else:
                print("  => P3 FAILED: Another stable state has higher entropy production")
        else:
            print("\nNo cooperative-dominated stable FP found.")
    else:
        print("\nNo stable fixed points found.")

    # Phase 4: Nash equilibrium analysis
    print("\n" + "=" * 60)
    print("PHASE 4: Nash Equilibrium Analysis")
    print("=" * 60)

    for fp in stable_fps:
        x = np.array(fp["x"])
        f = np.array(fp["f"])
        f_mean = fp["f_mean"]
        active = [i for i in range(n) if x[i] > 0.01]
        deviations = [abs(f[i] - f_mean) for i in active]
        is_nash = all(d < 0.1 for d in deviations)
        print("\nFP (dominant=%s): Nash=%s" % (names[fp["dominant"]], "YES" if is_nash else "NO"))
        for i in active:
            print("  %s: f=%.4f, <f>=%.4f, deviation=%.4f" % (names[i], f[i], f_mean, abs(f[i] - f_mean)))

    # Save results
    results = {
        "scenario": args.scenario,
        "description": scenario["description"],
        "species": names,
        "fixed_points": fps,
        "dynamics": {
            "x_final": dyn["x_final"],
            "sigma_final": dyn["sigma_final"],
            "sigma_trajectory": dyn["sigma"][::20],
            "t_subsample": dyn["t"][::20],
        },
        "p3_test": {
            "n_stable": len(stable_fps),
            "n_unstable": len(unstable_fps),
            "cooperative_dominates": any(fp["dominant"] == 0 for fp in stable_fps),
            "max_sigma_stable": max(fp["sigma"] for fp in stable_fps) if stable_fps else None,
        }
    }

    out_path = args.out or "/opt/fst_calculations/results/fst_ii/replicator_%s.json" % args.scenario
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)
    print("\nResults saved to: %s" % out_path)


if __name__ == "__main__":
    main()
