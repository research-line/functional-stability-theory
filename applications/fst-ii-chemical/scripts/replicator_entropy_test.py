#!/usr/bin/env python3
"""
FST-II: Replicator Dynamics Toy Proxy for Autocatalytic Networks
================================================================

Legacy exploratory model for the P3 idea in the FST-II paper:
three competing autocatalytic-network types are modeled as replicator dynamics.
This script does not measure thermodynamic entropy production. It reports a
fitness-dispersion / throughput proxy and must not be cited as P3 evidence.

Physics:
  - Replicator equation: dx_i/dt = x_i * (f_i(x) - <f>(x))
  - Fitness from catalytic rate constants (Vaidya et al. 2012)
  - Toy proxy: sigma_proxy = sum_i x_i * f_i * ln(f_i / <f>)
  - Nash equilibrium: interior ESS where no species can improve unilaterally

Guardrail:
  P3 requires independent observables for game status and dissipation
  (for example sequencing plus calorimetry). This script only illustrates
  candidate game dynamics and cannot confirm or falsify thermodynamic P3.

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
import argparse
from pathlib import Path

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


def entropy_proxy_rate(x, A, coop_bonus):
    """
    Fitness-dispersion proxy sigma_proxy = sum_i x_i * f_i * ln(f_i / <f>).

    This is a game/throughput proxy, not a thermodynamic entropy-production
    rate. It can be useful for toy-model ranking, but P3 requires an
    independent physical dissipation observable.
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
    support_tol = 1e-6
    stability_tol = 1e-6

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

        fp["sigma_proxy"] = entropy_proxy_rate(x, A, coop_bonus)
        fp["sigma"] = fp["sigma_proxy"]
        fp["f_mean"] = f_mean
        fp["f"] = f.tolist()

        J = compute_jacobian(x, A, coop_bonus)
        eigs = np.linalg.eigvals(J)
        eigs_sorted = sorted(eigs.real)
        max_real = float(max(eigs.real))
        support = x > support_tol
        invasion_channels = []
        neutral_invasion_channels = []
        for idx, gap in enumerate(f - f_mean):
            if support[idx]:
                continue
            gap = float(gap)
            if gap > stability_tol:
                invasion_channels.append({"species": idx, "fitness_gap": gap})
            elif abs(gap) <= stability_tol:
                neutral_invasion_channels.append({"species": idx, "fitness_gap": gap})

        fp["eigenvalues"] = [float(e) for e in eigs_sorted]
        fp["max_eigenvalue"] = max_real
        fp["n_positive_eigs"] = sum(1 for e in eigs if e.real > stability_tol)
        fp["invasion_channels"] = invasion_channels
        fp["neutral_invasion_channels"] = neutral_invasion_channels
        if max_real > stability_tol or invasion_channels:
            fp["stability_status"] = "unstable"
        elif abs(max_real) <= stability_tol or neutral_invasion_channels:
            fp["stability_status"] = "neutral_boundary"
        else:
            fp["stability_status"] = "asymptotically_stable"
        fp["stable"] = fp["stability_status"] == "asymptotically_stable"
        fp["dominant"] = int(np.argmax(x))

    return fixed_points


def run_dynamics(A, coop_bonus, names, x0=None, t_max=200):
    """Run replicator dynamics and track the fitness/throughput proxy."""
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
        sigmas.append(entropy_proxy_rate(x, A, coop_bonus))

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
    parser.add_argument("--seed", type=int, default=20260522,
                        help="Random seed for reproducible exploratory scans")
    args = parser.parse_args()
    np.random.seed(args.seed)

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
        status_labels = {
            "asymptotically_stable": "ASYMPTOTICALLY STABLE",
            "neutral_boundary": "NEUTRAL BOUNDARY",
            "unstable": "UNSTABLE",
        }
        stab = status_labels.get(fp["stability_status"], "UNCLASSIFIED")
        print("\n  FP-%d (%s, found %dx):" % (i + 1, stab, fp["count"]))
        for j in range(n):
            print("    %12s: x=%.6f, f=%.4f" % (names[j], fp["x"][j], fp["f"][j]))
        print("    Fitness/throughput proxy: sigma_proxy = %.6f" % fp["sigma_proxy"])
        print("    Mean fitness: <f> = %.4f" % fp["f_mean"])
        print("    Dominant: %s" % dom)
        print("    Eigenvalues: %s" % [round(e, 4) for e in fp["eigenvalues"]])
        print("    Positive eigenvalues: %d" % fp["n_positive_eigs"])
        if fp["neutral_invasion_channels"]:
            neutral_names = [names[item["species"]] for item in fp["neutral_invasion_channels"]]
            print("    Neutral boundary channels: %s" % neutral_names)
        if fp["invasion_channels"]:
            invasion_names = [names[item["species"]] for item in fp["invasion_channels"]]
            print("    Invasion channels: %s" % invasion_names)

    # Phase 2: Dynamics from uniform start
    print("\n" + "=" * 60)
    print("PHASE 2: Dynamics from uniform start")
    print("=" * 60)

    dyn = run_dynamics(A, coop, names)
    print("\nFinal state (t=%.0f):" % dyn["t"][-1])
    for name in names:
        print("  %12s: %.6f" % (name, dyn["x"][name][-1]))
    print("  Fitness/throughput proxy: %.6f" % dyn["sigma_final"])

    # Phase 3: P3 Test
    print("\n" + "=" * 60)
    print("PHASE 3: P3 GUARDRAIL -- toy proxy only, no thermodynamic confirmation")
    print("=" * 60)

    stable_fps = [fp for fp in fps if fp["stable"]]
    neutral_fps = [fp for fp in fps if fp.get("stability_status") == "neutral_boundary"]
    unstable_fps = [fp for fp in fps if fp.get("stability_status") == "unstable"]
    nonstable_fps = neutral_fps + unstable_fps

    if stable_fps:
        max_sigma_stable = max(fp["sigma"] for fp in stable_fps)
        best_stable = max(stable_fps, key=lambda fp: fp["sigma"])
        print("\nBest stable FP by proxy: %s (sigma_proxy = %.6f)" % (names[best_stable["dominant"]], best_stable["sigma_proxy"]))

        if nonstable_fps:
            best_nonstable = max(nonstable_fps, key=lambda fp: fp["sigma"])
            print("Best non-asymptotically-stable FP by proxy: %s (sigma_proxy = %.6f)" % (names[best_nonstable["dominant"]], best_nonstable["sigma_proxy"]))

        coop_fps = [fp for fp in stable_fps if fp["dominant"] == 0]
        if coop_fps:
            coop_sigma = max(fp["sigma"] for fp in coop_fps)
            print("\nToy proxy result: Cooperative ESS sigma_proxy = %.6f" % coop_sigma)
            print("                  Max stable sigma_proxy      = %.6f" % max_sigma_stable)
            print("  => P3 STATUS: NOT TESTED. This proxy is not independent calorimetric entropy production.")
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
        "random_seed": args.seed,
        "guardrail": (
            "sigma_proxy is a fitness/throughput proxy, not thermodynamic entropy production; "
            "this script is not P3 evidence."
        ),
        "species": names,
        "fixed_points": fps,
        "dynamics": {
            "x_final": dyn["x_final"],
            "sigma_final": dyn["sigma_final"],
            "sigma_trajectory": dyn["sigma"][::20],
            "t_subsample": dyn["t"][::20],
        },
        "p3_test": {
            "claim_status": "not_tested_toy_proxy_only",
            "n_stable": len(stable_fps),
            "n_neutral_boundary": sum(1 for fp in fps if fp.get("stability_status") == "neutral_boundary"),
            "n_unstable": len(unstable_fps),
            "cooperative_dominates": any(fp["dominant"] == 0 for fp in stable_fps),
            "max_sigma_stable": max(fp["sigma"] for fp in stable_fps) if stable_fps else None,
        }
    }

    base_dir = Path(__file__).resolve().parents[1]
    out_path = Path(args.out) if args.out else base_dir / "results" / f"replicator_{args.scenario}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    print("\nResults saved to: %s" % out_path)


if __name__ == "__main__":
    main()
