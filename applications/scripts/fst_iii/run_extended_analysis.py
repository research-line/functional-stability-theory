#!/usr/bin/env python3
"""
III-CALC-01: Extended Nash-Frustration Analysis for all 3 proteins.

Batch analysis wrapper for protein_fold_nash_pdb.py:
  - Self-consistency fit (train on itself)
  - Hessian, spectral radius rho(J), frustration map
  - Best-response dynamics: convergence check
  - Nash stability validation
  - Comparison table (console + LaTeX)

Proteins: HP35 (1YRF), Protein G (1PGA), 2XWR

Author: Lukas Geiger (with Claude, 2026-03-23)
Part of: FST-III Biological Game Theory
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

import sys
import argparse
import json
import numpy as np
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from protein_fold_nash_pdb import (
    AA_TO_IDX,
    load_ca_coords_and_seq, extract_phi_psi, build_contacts,
    fit_params_from_structure, energy_and_grad, numerical_hessian,
    frustration_map, analyze_frustration,
    best_response_dynamics, validate_nash_stability
)


# ---------------------------------------------------------------------------
# Analysis for a single protein
# ---------------------------------------------------------------------------

def analyze_protein(pdb_path: str, chain_id: str, label: str,
                    lr: float = 0.015, rcut: float = 8.0,
                    sigma: float = 1.0, br_sweeps: int = 300,
                    verbose: bool = True) -> dict:
    """
    Full Nash analysis for a single protein.

    Steps:
    1. Load structure
    2. Fit Nash potential (self-consistency)
    3. Compute Hessian, rho(J), frustration map
    4. Run best-response dynamics from native
    5. Validate Nash stability (perturbed start)
    """
    print(f"\n{'='*70}")
    print(f"PROTEIN: {label}")
    print(f"PDB: {pdb_path}, Chain: {chain_id}")
    print(f"{'='*70}")

    # 1. Load structure
    ca, seq = load_ca_coords_and_seq(pdb_path, chain_id)
    phi, psi = extract_phi_psi(pdb_path, chain_id)
    edges, r = build_contacts(ca, r_cut=rcut)
    N = len(seq)
    n_contacts = len(edges)

    print(f"  Residues: {N}")
    print(f"  Sequence: {seq}")
    print(f"  Contacts: {n_contacts}")

    # 2. Fit Nash potential (self-consistency)
    print("\n  [Step 2] Fitting Nash potential...")
    params = fit_params_from_structure(seq, phi, psi, edges, r, sigma=sigma)

    aa_idx = np.array([AA_TO_IDX.get(a, 0) for a in seq], dtype=int)
    F, gphi, gpsi = energy_and_grad(
        phi, psi, aa_idx, edges, r,
        params.mu_phi, params.mu_psi,
        params.k_phi, params.k_psi,
        params.w0, params.w_rbf,
        params.r_centers, sigma
    )
    grad_norm = float(np.linalg.norm(np.concatenate([gphi, gpsi])))
    print(f"  F(theta*) = {F:.6f}")
    print(f"  |grad F| = {grad_norm:.6e}")

    # 3. Hessian and spectral analysis
    print("\n  [Step 3] Computing Hessian...")
    H = numerical_hessian(phi, psi, aa_idx, edges, r, params, sigma=sigma)

    # Eigenvalues of Hessian (for convexity)
    evals_H = np.linalg.eigvalsh(H)
    n_positive = int(np.sum(evals_H > 1e-8))
    frac_positive = float(n_positive / len(evals_H))

    # Jacobian J = I - lr * H
    J = np.eye(H.shape[0]) - lr * H
    evals_J = np.linalg.eigvals(J)
    rho_J = float(np.max(np.abs(evals_J)))
    n_unstable = int(np.sum(np.abs(evals_J) >= 1.0))

    print(f"  Hessian: {n_positive}/{len(evals_H)} positive eigenvalues ({frac_positive:.1%})")
    print(f"  rho(J) = {rho_J:.4f} ({'STABLE' if rho_J < 1.0 else 'UNSTABLE'})")
    print(f"  Unstable modes: {n_unstable}/{2*N}")

    # 4. Frustration analysis
    print("\n  [Step 4] Frustration analysis...")
    frust_result = analyze_frustration(seq, H, lr=lr, verbose=verbose)
    frust_scores = np.array(frust_result["frustration_scores"])
    mean_frust = float(np.mean(frust_scores))
    max_frust = float(np.max(frust_scores))

    # 5. Best-response dynamics from native (convergence check)
    print("\n  [Step 5] Best-response dynamics...")
    br_result = best_response_dynamics(
        seq, edges, r, params,
        phi, psi,
        sigma=sigma, lr=0.05,
        sweeps=br_sweeps, verbose=verbose
    )
    br_converged = br_result["converged"]
    br_final_F = br_result["F"]
    br_sweeps_used = br_result["sweeps"]

    print(f"  Converged: {br_converged} (after {br_sweeps_used} sweeps)")
    print(f"  Final F: {br_final_F:.6f}")

    # 6. Nash stability validation
    print("\n  [Step 6] Nash stability validation...")
    nash_result = validate_nash_stability(
        seq, edges, r, params,
        phi, psi,
        sigma=sigma,
        perturbation=0.1,
        sweeps=br_sweeps,
        verbose=verbose
    )

    # Compile results
    result = {
        "protein": label,
        "pdb_path": pdb_path,
        "chain": chain_id,
        "sequence": seq,
        "n_residues": N,
        "n_contacts": n_contacts,

        # Fit quality
        "F_at_native": float(F),
        "grad_norm_at_native": grad_norm,

        # Hessian / spectral
        "hessian_positive_frac": frac_positive,
        "hessian_eigenvalue_range": {
            "min": float(np.min(evals_H)),
            "median": float(np.median(evals_H)),
            "max": float(np.max(evals_H))
        },
        "rho_J": rho_J,
        "n_unstable_modes": n_unstable,
        "spectral_stable": rho_J < 1.0,

        # Frustration
        "mean_frustration": mean_frust,
        "max_frustration": max_frust,
        "frustration_scores": frust_scores.tolist(),
        "high_frustration_residues": frust_result["high_frustration_residues"],

        # Best-response dynamics
        "br_converged": br_converged,
        "br_sweeps_used": br_sweeps_used,
        "br_final_F": br_final_F,

        # Nash validation
        "nash_stable": nash_result.get("is_nash_stable", False),
        "nash_phi_rms": nash_result.get("phi_rms_vs_native", None),
        "nash_psi_rms": nash_result.get("psi_rms_vs_native", None),
        "nash_converged": nash_result.get("converged", False),
    }

    return result


# ---------------------------------------------------------------------------
# LaTeX table generation
# ---------------------------------------------------------------------------

def generate_latex_table(results: list) -> str:
    """Generate a LaTeX comparison table."""
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Extended Nash-Frustration Analysis -- Protein Comparison}")
    lines.append(r"\label{tab:nash_extended}")
    lines.append(r"\begin{tabular}{lcccccc}")
    lines.append(r"\toprule")
    lines.append(r"Protein & $N_{res}$ & $\rho(J)$ & $\bar{f}$ & $f_{max}$ & BR conv. & Nash stable \\")
    lines.append(r"\midrule")

    for r in results:
        stable_sym = r"\checkmark" if r["spectral_stable"] else r"\times"
        br_sym = r"\checkmark" if r["br_converged"] else r"\times"
        nash_sym = r"\checkmark" if r["nash_stable"] else r"\times"
        lines.append(
            f"  {r['protein']} & {r['n_residues']} & "
            f"{r['rho_J']:.4f} & {r['mean_frustration']:.4f} & "
            f"{r['max_frustration']:.4f} & {br_sym} & {nash_sym} \\\\"
        )

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


def generate_extended_latex(results: list) -> str:
    """Generate extended LaTeX table with all metrics."""
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Full Nash Stability Metrics}")
    lines.append(r"\label{tab:nash_full}")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{lcccccccc}")
    lines.append(r"\toprule")
    lines.append(r"Protein & $N$ & Contacts & $|{\nabla F}|$ & $\rho(J)$ & "
                 r"Unstable & $\bar{f}$ & $\phi_{rms}$ & $\psi_{rms}$ \\")
    lines.append(r"\midrule")

    for r in results:
        phi_rms = f"{r['nash_phi_rms']:.4f}" if r['nash_phi_rms'] is not None else "---"
        psi_rms = f"{r['nash_psi_rms']:.4f}" if r['nash_psi_rms'] is not None else "---"
        lines.append(
            f"  {r['protein']} & {r['n_residues']} & {r['n_contacts']} & "
            f"{r['grad_norm_at_native']:.2e} & {r['rho_J']:.4f} & "
            f"{r['n_unstable_modes']} & {r['mean_frustration']:.4f} & "
            f"{phi_rms} & {psi_rms} \\\\"
        )

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="III-CALC-01: Extended Nash-Frustration Analysis"
    )
    ap.add_argument("--data-dir", required=True,
                    help="Directory containing PDB files (1YRF.pdb, 1PGA.pdb, 2XWR.pdb)")
    ap.add_argument("--output-dir", default=None,
                    help="Output directory (default: same as data-dir)")
    ap.add_argument("--lr", type=float, default=0.015,
                    help="Learning rate for Jacobian")
    ap.add_argument("--rcut", type=float, default=8.0,
                    help="Contact distance cutoff")
    ap.add_argument("--sigma", type=float, default=1.0,
                    help="RBF sigma")
    ap.add_argument("--br-sweeps", type=int, default=300,
                    help="Max sweeps for best-response dynamics")
    ap.add_argument("--quiet", action="store_true",
                    help="Reduce per-residue output")
    args = ap.parse_args()

    if args.output_dir is None:
        args.output_dir = args.data_dir
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 70)
    print("III-CALC-01: Extended Nash-Frustration Analysis")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 70)

    # Protein definitions
    proteins = [
        ("1YRF.pdb", "A", "HP35 (1YRF)"),
        ("1PGA.pdb", "A", "Protein G (1PGA)"),
        ("2XWR.pdb", "A", "2XWR"),
    ]

    all_results = []
    for pdb_file, chain, label in proteins:
        pdb_path = os.path.join(args.data_dir, pdb_file)
        if not os.path.exists(pdb_path):
            print(f"\nWARNING: {pdb_path} not found, skipping {label}")
            continue

        result = analyze_protein(
            pdb_path, chain, label,
            lr=args.lr, rcut=args.rcut, sigma=args.sigma,
            br_sweeps=args.br_sweeps,
            verbose=not args.quiet
        )
        all_results.append(result)

    if not all_results:
        print("ERROR: No PDB files found. Check --data-dir.")
        sys.exit(1)

    # ---------------------------------------------------------------------------
    # Comparison table (console)
    # ---------------------------------------------------------------------------
    print("\n\n" + "=" * 70)
    print("COMPARISON TABLE")
    print("=" * 70)

    header = (f"  {'Protein':>18s} | {'N_res':>5s} | {'Contacts':>8s} | {'rho(J)':>8s} | "
              f"{'mean_f':>7s} | {'max_f':>7s} | {'BR conv':>8s} | {'Nash':>6s}")
    print(header)
    print("  " + "-" * len(header))

    for r in all_results:
        br = "YES" if r["br_converged"] else "NO"
        nash = "YES" if r["nash_stable"] else "NO"
        stable = "S" if r["spectral_stable"] else "U"
        print(f"  {r['protein']:>18s} | {r['n_residues']:5d} | {r['n_contacts']:8d} | "
              f"{r['rho_J']:8.4f}{stable} | {r['mean_frustration']:7.4f} | "
              f"{r['max_frustration']:7.4f} | {br:>8s} | {nash:>6s}")

    # ---------------------------------------------------------------------------
    # LaTeX output
    # ---------------------------------------------------------------------------
    latex_simple = generate_latex_table(all_results)
    latex_full = generate_extended_latex(all_results)

    print("\n\n" + "=" * 70)
    print("LATEX TABLE (compact)")
    print("=" * 70)
    print(latex_simple)

    print("\n" + "=" * 70)
    print("LATEX TABLE (extended)")
    print("=" * 70)
    print(latex_full)

    # Save LaTeX
    latex_path = os.path.join(args.output_dir, "nash_comparison_table.tex")
    with open(latex_path, 'w', encoding='utf-8') as f:
        f.write("% Auto-generated by run_extended_analysis.py\n")
        f.write(f"% Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        f.write("% Compact table\n")
        f.write(latex_simple)
        f.write("\n\n% Extended table\n")
        f.write(latex_full)

    print(f"\nLaTeX saved: {latex_path}")

    # ---------------------------------------------------------------------------
    # Save JSON
    # ---------------------------------------------------------------------------
    output = {
        "metadata": {
            "description": "III-CALC-01: Extended Nash-Frustration Analysis",
            "author": "Lukas Geiger",
            "date": datetime.now().strftime("%Y-%m-%d"),
            "parameters": {
                "lr": args.lr,
                "rcut": args.rcut,
                "sigma": args.sigma,
                "br_sweeps": args.br_sweeps
            }
        },
        "proteins": all_results,
        "latex_compact": latex_simple,
        "latex_extended": latex_full
    }

    json_path = os.path.join(args.output_dir, "extended_analysis_results.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2)
    print(f"JSON saved: {json_path}")

    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    all_stable = all(r["spectral_stable"] for r in all_results)
    all_nash = all(r["nash_stable"] for r in all_results)
    all_br = all(r["br_converged"] for r in all_results)

    n_total = len(all_results)
    n_stable = sum(1 for r in all_results if r["spectral_stable"])
    n_nash = sum(1 for r in all_results if r["nash_stable"])
    n_br = sum(1 for r in all_results if r["br_converged"])

    print(f"\n  Spectral stability (rho(J) < 1): {n_stable}/{n_total}")
    print(f"  Best-response convergence:       {n_br}/{n_total}")
    print(f"  Nash stability (validated):      {n_nash}/{n_total}")

    if all_stable and all_nash:
        print("\n  => ALL proteins are Nash-stable. Consistent with FST-III prediction P1.")
    elif all_stable:
        print("\n  => All spectrally stable, but Nash validation incomplete for some.")
    else:
        unstable = [r["protein"] for r in all_results if not r["spectral_stable"]]
        print(f"\n  => Spectrally unstable: {', '.join(unstable)}")
        print("     These show structural frustration (game-theoretic coordination breakdown).")

    rho_values = [r["rho_J"] for r in all_results]
    print(f"\n  rho(J) range: [{min(rho_values):.4f}, {max(rho_values):.4f}]")
    print(f"  Mean frustration range: [{min(r['mean_frustration'] for r in all_results):.4f}, "
          f"{max(r['mean_frustration'] for r in all_results):.4f}]")


if __name__ == "__main__":
    main()
