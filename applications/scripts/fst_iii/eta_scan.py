#!/usr/bin/env python3
"""
QUANT-III-2: Systematic eta-scan for Nash frustration.

For each PDB structure (1YRF, 1PGA, 2XWR) and a range of eta values,
compute the best-response Jacobian J = I - eta*H, spectral radius rho(J),
and frustration map. Output JSON results and matplotlib plots.

Author: Lukas Geiger (with Claude, 2026-03-23)
Part of: FST-III Biological Game Theory
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

import sys
import argparse
import json
import numpy as np

# Import from main module in same directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from protein_fold_nash_pdb import (
    AA_TO_IDX,
    load_ca_coords_and_seq, extract_phi_psi, build_contacts,
    fit_params_from_structure, numerical_hessian,
    frustration_map, analyze_frustration
)


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def eta_scan_single_protein(pdb_path: str, chain_id: str, pdb_label: str,
                            eta_values: list, rcut: float = 8.0,
                            sigma: float = 1.0) -> dict:
    """
    Run eta-scan for a single protein structure.

    For each eta:
      - Compute J = I - eta * H
      - Spectral radius rho(J)
      - Frustration map from unstable modes

    Returns dict with all results for this protein.
    """
    print(f"\n{'='*60}")
    print(f"Protein: {pdb_label} ({pdb_path}, chain {chain_id})")
    print(f"{'='*60}")

    # Load structure
    ca, seq = load_ca_coords_and_seq(pdb_path, chain_id)
    phi, psi = extract_phi_psi(pdb_path, chain_id)
    edges, r = build_contacts(ca, r_cut=rcut)
    N = len(seq)

    print(f"  Residues: {N}, Contacts: {len(edges)}")

    # Fit parameters (self-consistency)
    print("  Fitting Nash potential...")
    params = fit_params_from_structure(seq, phi, psi, edges, r, sigma=sigma)

    # Compute Hessian once
    print("  Computing Hessian...")
    aa_idx = np.array([AA_TO_IDX.get(a, 0) for a in seq], dtype=int)
    H = numerical_hessian(phi, psi, aa_idx, edges, r, params, sigma=sigma)

    # Scan eta values
    results = []
    print(f"\n  {'eta':>8s} | {'rho(J)':>10s} | {'mean_frust':>11s} | {'max_frust':>10s} | {'n_unstable':>11s} | high_frust_residues")
    print("  " + "-" * 90)

    for eta in eta_values:
        J = np.eye(H.shape[0]) - eta * H
        eigvals = np.linalg.eigvals(J)
        rho = float(np.max(np.abs(eigvals)))
        n_unstable = int(np.sum(np.abs(eigvals) >= 1.0))

        # Frustration map
        frust = frustration_map(H, eta, N)
        mean_frust = float(np.mean(frust))
        max_frust = float(np.max(frust))

        # High-frustration residues (> 0.5)
        high_frust_res = [(int(i), seq[i], float(frust[i]))
                          for i in range(N) if frust[i] > 0.5]

        result_entry = {
            "eta": float(eta),
            "rho_J": rho,
            "mean_frustration": mean_frust,
            "max_frustration": max_frust,
            "n_unstable_modes": n_unstable,
            "high_frustration_residues": [
                {"index": idx, "aa": aa, "score": sc}
                for idx, aa, sc in high_frust_res
            ],
            "frustration_scores": frust.tolist()
        }
        results.append(result_entry)

        high_res_str = ", ".join(f"{aa}{idx}" for idx, aa, _ in high_frust_res[:5])
        if len(high_frust_res) > 5:
            high_res_str += f" (+{len(high_frust_res)-5} more)"
        print(f"  {eta:8.4f} | {rho:10.4f} | {mean_frust:11.4f} | {max_frust:10.4f} | {n_unstable:11d} | {high_res_str}")

    return {
        "pdb_label": pdb_label,
        "pdb_path": pdb_path,
        "chain": chain_id,
        "sequence": seq,
        "n_residues": N,
        "n_contacts": int(len(edges)),
        "eta_scan": results
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_eta_scan(all_results: list, output_dir: str):
    """Plot eta vs rho(J) for all proteins."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping plot.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    # --- Plot 1: eta vs rho(J) ---
    ax1 = axes[0]
    for i, prot in enumerate(all_results):
        etas = [r["eta"] for r in prot["eta_scan"]]
        rhos = [r["rho_J"] for r in prot["eta_scan"]]
        ax1.plot(etas, rhos, 'o-', color=colors[i % len(colors)],
                 label=prot["pdb_label"], linewidth=2, markersize=6)

    ax1.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label=r'$\rho(J)=1$ (stability boundary)')
    ax1.set_xlabel(r'$\eta$ (learning rate)', fontsize=13)
    ax1.set_ylabel(r'$\rho(J)$ (spectral radius)', fontsize=13)
    ax1.set_title('Nash Stability vs. Learning Rate', fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xscale('log')

    # --- Plot 2: eta vs mean frustration ---
    ax2 = axes[1]
    for i, prot in enumerate(all_results):
        etas = [r["eta"] for r in prot["eta_scan"]]
        mean_frust = [r["mean_frustration"] for r in prot["eta_scan"]]
        ax2.plot(etas, mean_frust, 's-', color=colors[i % len(colors)],
                 label=prot["pdb_label"], linewidth=2, markersize=6)

    ax2.set_xlabel(r'$\eta$ (learning rate)', fontsize=13)
    ax2.set_ylabel('Mean Frustration', fontsize=13)
    ax2.set_title('Nash Frustration vs. Learning Rate', fontsize=14)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log')

    plt.tight_layout()
    plot_path = os.path.join(output_dir, "eta_scan_results.png")
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nPlot saved: {plot_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="QUANT-III-2: Systematic eta-scan for Nash frustration"
    )
    ap.add_argument("--data-dir", required=True,
                    help="Directory containing PDB files (1YRF.pdb, 1PGA.pdb, 2XWR.pdb)")
    ap.add_argument("--output-dir", default=None,
                    help="Output directory for JSON and plots (default: same as data-dir)")
    ap.add_argument("--rcut", type=float, default=8.0,
                    help="Contact distance cutoff in Angstrom")
    ap.add_argument("--sigma", type=float, default=1.0,
                    help="RBF sigma for distance weighting")
    args = ap.parse_args()

    if args.output_dir is None:
        args.output_dir = args.data_dir
    os.makedirs(args.output_dir, exist_ok=True)

    # PDB files to process
    proteins = [
        ("1YRF.pdb", "A", "HP35 (1YRF)"),
        ("1PGA.pdb", "A", "Protein G (1PGA)"),
        ("2XWR.pdb", "A", "2XWR"),
    ]

    # eta values to scan
    eta_values = [0.001, 0.005, 0.01, 0.015, 0.02, 0.05, 0.1]

    all_results = []
    for pdb_file, chain, label in proteins:
        pdb_path = os.path.join(args.data_dir, pdb_file)
        if not os.path.exists(pdb_path):
            print(f"WARNING: {pdb_path} not found, skipping {label}")
            continue

        result = eta_scan_single_protein(
            pdb_path, chain, label,
            eta_values=eta_values,
            rcut=args.rcut,
            sigma=args.sigma
        )
        all_results.append(result)

    if not all_results:
        print("ERROR: No PDB files found. Check --data-dir.")
        sys.exit(1)

    # Save JSON
    output = {
        "metadata": {
            "description": "QUANT-III-2: Systematic eta-scan for Nash frustration",
            "author": "Lukas Geiger",
            "date": "2026-03-23",
            "eta_values": eta_values,
            "parameters": {
                "rcut": args.rcut,
                "sigma": args.sigma
            }
        },
        "proteins": all_results
    }

    json_path = os.path.join(args.output_dir, "eta_scan_results.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2)
    print(f"\nJSON saved: {json_path}")

    # Plot
    plot_eta_scan(all_results, args.output_dir)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for prot in all_results:
        print(f"\n{prot['pdb_label']} ({prot['n_residues']} residues, {prot['n_contacts']} contacts):")
        for r in prot["eta_scan"]:
            stable = "STABLE" if r["rho_J"] < 1.0 else "UNSTABLE"
            print(f"  eta={r['eta']:.3f}: rho(J)={r['rho_J']:.4f} [{stable}], "
                  f"mean_frust={r['mean_frustration']:.4f}, n_unstable={r['n_unstable_modes']}")


if __name__ == "__main__":
    main()
