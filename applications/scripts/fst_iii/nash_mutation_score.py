#!/usr/bin/env python3
"""
nash_mutation_score.py
======================
Score pathogenicity of mutations using Nash stability analysis.

Pipeline:
1. Load WT structure (PDB/AlphaFold)
2. Fit Nash potential to WT
3. For each mutation: perturb local region, re-equilibrate, measure Δρ(J)
4. Output scores for ROC/PR-AUC analysis

Score candidates:
- Δρ(J): Change in spectral radius (stability loss)
- ΔN_unstable: Increase in unstable modes
- ΔFrust_max: Hotspot intensification
- ΔFrust_sum: Global coordination burden

Reference: FST-III Biological Game Theory (Geiger, 2026)

Usage:
    python nash_mutation_score.py --pdb structure.pdb --chain A --mutations variants.csv --out scores.json
"""

import argparse
import json
import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import csv

# Import from main module
from protein_fold_nash_pdb import (
    AA20, AA_TO_IDX,
    load_ca_coords_and_seq, extract_phi_psi, build_contacts,
    fit_params_from_structure, energy_and_grad, numerical_hessian,
    frustration_map, best_response_dynamics, LearnedParams
)


@dataclass
class MutationScore:
    """Score for a single mutation."""
    mutation: str           # e.g., "A15V"
    position: int           # 0-indexed
    wt_aa: str
    mut_aa: str
    # Primary clinical scores
    theta_drift_patch: float    # Δθ-Drift in local patch (main score)
    theta_drift_site: float     # Δθ-Drift at mutation site
    delta_frust_patch: float    # ΔFrustration in patch
    delta_frust_site: float     # ΔFrustration at mutation site
    # Secondary scores
    delta_rho: float            # Δρ(J) - for reference
    delta_unstable: int         # ΔN_unstable
    wt_rho: float
    mut_rho: float
    converged: bool
    patch_size: int             # Number of residues in patch


def parse_mutation(mutation: str, pdb_offset: int = 0) -> Tuple[str, int, str]:
    """
    Parse mutation string like 'G105R' -> ('G', 14, 'R') when pdb_offset=91

    Args:
        mutation: Mutation string (e.g., 'G105R')
        pdb_offset: First residue number in PDB (e.g., 91 for p53 DBD)
                    Position 105 becomes index 105-91=14

    Returns:
        wt_aa: Wild-type amino acid
        pos: 0-indexed position in the loaded sequence
        mut_aa: Mutant amino acid
    """
    wt_aa = mutation[0]
    mut_aa = mutation[-1]
    pdb_pos = int(mutation[1:-1])
    pos = pdb_pos - pdb_offset  # 0-indexed: pos 105 with offset 91 = index 14
    return wt_aa, pos, mut_aa


def apply_mutation_to_sequence(seq: str, pos: int, new_aa: str) -> str:
    """Apply mutation to sequence."""
    seq_list = list(seq)
    seq_list[pos] = new_aa
    return "".join(seq_list)


def compute_stability_metrics(phi: np.ndarray, psi: np.ndarray, seq: str,
                               edges: np.ndarray, r: np.ndarray,
                               params: LearnedParams, sigma: float = 1.0,
                               lr: float = 0.015) -> Dict:
    """
    Compute Nash stability metrics for a configuration.
    """
    aa_idx = np.array([AA_TO_IDX.get(a, 0) for a in seq], dtype=int)
    N = len(seq)

    # Hessian
    H = numerical_hessian(phi, psi, aa_idx, edges, r, params, sigma=sigma)

    # Jacobian and spectral radius
    J = np.eye(H.shape[0]) - lr * H
    eigvals = np.linalg.eigvals(J)
    rho = float(np.max(np.abs(eigvals)))
    n_unstable = int(np.sum(np.abs(eigvals) >= 1.0))

    # Frustration
    frust = frustration_map(H, lr, N)

    return {
        "rho": rho,
        "n_unstable": n_unstable,
        "frust_max": float(np.max(frust)),
        "frust_mean": float(np.mean(frust)),
        "frustration": frust
    }


def local_patch_relax(phi: np.ndarray, psi: np.ndarray, seq: str,
                      edges: np.ndarray, r: np.ndarray, params: LearnedParams,
                      center: int, patch_radius: int = 15,
                      sigma: float = 1.0, lr: float = 0.05,
                      sweeps: int = 50) -> Tuple[np.ndarray, np.ndarray, bool]:
    """
    Relax only a local patch around center position.
    Rest of the protein stays fixed at input configuration.

    This is the clinically correct approach: mutations affect local coordination,
    not global structure.
    """
    N = len(seq)
    aa_idx = np.array([AA_TO_IDX.get(a, 0) for a in seq], dtype=int)

    phi_out = phi.copy()
    psi_out = psi.copy()

    # Define patch
    patch_start = max(0, center - patch_radius)
    patch_end = min(N, center + patch_radius + 1)
    patch_indices = set(range(patch_start, patch_end))

    converged = False
    tol = 1e-4

    for sweep in range(sweeps):
        max_update = 0.0

        # Only update residues in patch
        for i in range(patch_start, patch_end):
            _, gphi, gpsi = energy_and_grad(
                phi_out, psi_out, aa_idx, edges, r,
                params.mu_phi, params.mu_psi,
                params.k_phi, params.k_psi,
                params.w0, params.w_rbf,
                params.r_centers, sigma
            )

            # Update only this residue
            dphi = -lr * gphi[i]
            dpsi = -lr * gpsi[i]
            phi_out[i] += dphi
            psi_out[i] += dpsi

            max_update = max(max_update, abs(dphi), abs(dpsi))

        # Wrap angles
        phi_out = (phi_out + np.pi) % (2 * np.pi) - np.pi
        psi_out = (psi_out + np.pi) % (2 * np.pi) - np.pi

        if max_update < tol:
            converged = True
            break

    return phi_out, psi_out, converged


def angular_distance(a: float, b: float) -> float:
    """Compute angular distance accounting for periodicity."""
    d = a - b
    return abs(np.arctan2(np.sin(d), np.cos(d)))


def compute_theta_drift(phi1: np.ndarray, psi1: np.ndarray,
                        phi2: np.ndarray, psi2: np.ndarray,
                        indices: Optional[List[int]] = None) -> float:
    """
    Compute mean θ-drift: sqrt(mean(Δφ² + Δψ²)) over specified indices.
    This is the primary clinical score.
    """
    if indices is None:
        indices = list(range(len(phi1)))

    drift_sq = 0.0
    for i in indices:
        dphi = angular_distance(phi1[i], phi2[i])
        dpsi = angular_distance(psi1[i], psi2[i])
        drift_sq += dphi**2 + dpsi**2

    return float(np.sqrt(drift_sq / len(indices)))


def score_mutation(mutation: str, wt_seq: str, wt_phi: np.ndarray, wt_psi: np.ndarray,
                   edges: np.ndarray, r: np.ndarray, params: LearnedParams,
                   wt_metrics: Dict, sigma: float = 1.0, lr: float = 0.015,
                   patch_radius: int = 15, relax_sweeps: int = 50,
                   pdb_offset: int = 0) -> MutationScore:
    """
    Score a single mutation using local patch relaxation and θ-drift.

    Clinical scoring approach:
    1. Apply mutation to sequence
    2. Relax only local patch (rest fixed at WT)
    3. Measure θ-drift from WT configuration
    4. Measure frustration shift

    Higher θ-drift = mutation disrupts local coordination = likely pathogenic

    Args:
        mutation: Mutation string (e.g., 'A15V')
        wt_seq: Wild-type sequence
        wt_phi, wt_psi: Wild-type angles (reference)
        edges, r: Contact graph
        params: Learned potential parameters
        wt_metrics: Pre-computed WT metrics
        patch_radius: Residues around mutation (±radius)
        relax_sweeps: Number of local relaxation sweeps
    """
    wt_aa, pos, mut_aa = parse_mutation(mutation, pdb_offset)
    N = len(wt_seq)

    # Validate
    if pos < 0 or pos >= N:
        raise ValueError(f"Position {pos+1} out of range for sequence length {N}")
    if wt_seq[pos] != wt_aa:
        raise ValueError(f"WT amino acid mismatch: expected {wt_aa}, found {wt_seq[pos]} at position {pos+1}")

    # Apply mutation
    mut_seq = apply_mutation_to_sequence(wt_seq, pos, mut_aa)

    # Start from WT configuration
    mut_phi = wt_phi.copy()
    mut_psi = wt_psi.copy()

    # Small perturbation at mutation site to break symmetry
    np.random.seed(hash(mutation) % (2**31))
    mut_phi[pos] += 0.05 * np.random.randn()
    mut_psi[pos] += 0.05 * np.random.randn()

    # Local patch relaxation (rest stays fixed)
    mut_phi, mut_psi, converged = local_patch_relax(
        mut_phi, mut_psi, mut_seq, edges, r, params,
        center=pos, patch_radius=patch_radius,
        sigma=sigma, lr=0.05, sweeps=relax_sweeps
    )

    # Define patch indices
    patch_start = max(0, pos - patch_radius)
    patch_end = min(N, pos + patch_radius + 1)
    patch_indices = list(range(patch_start, patch_end))

    # Compute θ-drift scores
    theta_drift_patch = compute_theta_drift(wt_phi, wt_psi, mut_phi, mut_psi, patch_indices)
    theta_drift_site = compute_theta_drift(wt_phi, wt_psi, mut_phi, mut_psi, [pos])

    # Compute mutant metrics for frustration
    mut_metrics = compute_stability_metrics(
        mut_phi, mut_psi, mut_seq, edges, r, params, sigma, lr
    )

    # Frustration shifts
    wt_frust = np.array(wt_metrics["frustration"])
    mut_frust = np.array(mut_metrics["frustration"])

    delta_frust_patch = float(np.sum(mut_frust[patch_start:patch_end]) -
                               np.sum(wt_frust[patch_start:patch_end]))
    delta_frust_site = float(mut_frust[pos] - wt_frust[pos])

    # Stability metrics (for reference)
    delta_rho = mut_metrics["rho"] - wt_metrics["rho"]
    delta_unstable = mut_metrics["n_unstable"] - wt_metrics["n_unstable"]

    return MutationScore(
        mutation=mutation,
        position=pos,
        wt_aa=wt_aa,
        mut_aa=mut_aa,
        theta_drift_patch=theta_drift_patch,
        theta_drift_site=theta_drift_site,
        delta_frust_patch=delta_frust_patch,
        delta_frust_site=delta_frust_site,
        delta_rho=delta_rho,
        delta_unstable=delta_unstable,
        wt_rho=wt_metrics["rho"],
        mut_rho=mut_metrics["rho"],
        converged=converged,
        patch_size=len(patch_indices)
    )


def load_mutations_from_csv(csv_path: str) -> List[Dict]:
    """
    Load mutations from CSV file.
    Expected columns: mutation, label (optional: pathogenic/benign)
    """
    mutations = []
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mutations.append({
                "mutation": row.get("mutation", row.get("Mutation", "")),
                "label": row.get("label", row.get("Label", "unknown"))
            })
    return mutations


def main():
    ap = argparse.ArgumentParser(
        description="Score mutations using Nash stability analysis"
    )
    ap.add_argument("--pdb", required=True, help="PDB file for WT structure")
    ap.add_argument("--chain", default="A", help="Chain ID")
    ap.add_argument("--mutations", required=True,
                    help="CSV file with mutations (column: mutation) or comma-separated list")
    ap.add_argument("--out", default="mutation_scores.json", help="Output JSON file")
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--lr", type=float, default=0.015)
    ap.add_argument("--rcut", type=float, default=8.0)
    ap.add_argument("--patch-radius", type=int, default=15,
                    help="Radius of local patch for relaxation (±residues)")
    ap.add_argument("--relax-sweeps", type=int, default=50,
                    help="Number of local relaxation sweeps")
    ap.add_argument("--pdb-offset", type=int, default=0,
                    help="PDB residue numbering offset (e.g., 91 for p53 DBD)")
    args = ap.parse_args()

    # Load WT structure
    print(f"Loading WT structure: {args.pdb}, chain {args.chain}")
    ca, seq = load_ca_coords_and_seq(args.pdb, args.chain)
    phi, psi = extract_phi_psi(args.pdb, args.chain)
    edges, r = build_contacts(ca, r_cut=args.rcut)

    print(f"  Sequence length: {len(seq)}")
    print(f"  Contacts: {len(edges)}")

    # Fit parameters to WT
    print("Fitting Nash potential to WT structure...")
    params = fit_params_from_structure(seq, phi, psi, edges, r)

    # Compute WT metrics
    print("Computing WT stability metrics...")
    wt_metrics = compute_stability_metrics(phi, psi, seq, edges, r, params,
                                            args.sigma, args.lr)
    print(f"  WT ρ(J) = {wt_metrics['rho']:.4f}")
    print(f"  WT unstable modes: {wt_metrics['n_unstable']}")

    # Load mutations
    if args.mutations.endswith('.csv'):
        mutations = load_mutations_from_csv(args.mutations)
    else:
        # Comma-separated list
        mutations = [{"mutation": m.strip(), "label": "unknown"}
                     for m in args.mutations.split(",")]

    print(f"\nScoring {len(mutations)} mutations...")

    # Score each mutation
    results = {
        "wt": {
            "pdb": args.pdb,
            "chain": args.chain,
            "sequence": seq,
            "rho": wt_metrics["rho"],
            "n_unstable": wt_metrics["n_unstable"],
            "frust_max": wt_metrics["frust_max"],
        },
        "mutations": []
    }

    for i, mut_info in enumerate(mutations):
        mutation = mut_info["mutation"]
        label = mut_info["label"]

        try:
            score = score_mutation(
                mutation, seq, phi, psi, edges, r, params, wt_metrics,
                args.sigma, args.lr,
                patch_radius=args.patch_radius,
                relax_sweeps=args.relax_sweeps,
                pdb_offset=args.pdb_offset
            )

            result = {
                "mutation": mutation,
                "label": label,
                "position": score.position + 1,  # 1-indexed for output
                "theta_drift_patch": score.theta_drift_patch,
                "theta_drift_site": score.theta_drift_site,
                "delta_frust_patch": score.delta_frust_patch,
                "delta_frust_site": score.delta_frust_site,
                "delta_rho": score.delta_rho,
                "delta_unstable": score.delta_unstable,
                "patch_size": score.patch_size,
                "converged": score.converged,
            }
            results["mutations"].append(result)

            # Progress - θ-drift is now the primary score
            print(f"  [{i+1}/{len(mutations)}] {mutation}: θ-drift={score.theta_drift_patch:.4f} ({label})")

        except Exception as e:
            print(f"  [{i+1}/{len(mutations)}] {mutation}: ERROR - {e}")
            results["mutations"].append({
                "mutation": mutation,
                "label": label,
                "error": str(e)
            })

    # Save results
    with open(args.out, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)

    print(f"\nWrote: {args.out}")

    # Summary statistics
    valid_scores = [m for m in results["mutations"] if "theta_drift_patch" in m]
    if valid_scores:
        pathogenic = [m for m in valid_scores if m["label"].lower() in ["pathogenic", "pathogen", "p", "1"]]
        benign = [m for m in valid_scores if m["label"].lower() in ["benign", "b", "0"]]

        if pathogenic and benign:
            path_drift = np.mean([m["theta_drift_patch"] for m in pathogenic])
            ben_drift = np.mean([m["theta_drift_patch"] for m in benign])
            print(f"\n=== Summary (θ-drift as primary score) ===")
            print(f"Pathogenic (n={len(pathogenic)}): mean θ-drift = {path_drift:.4f}")
            print(f"Benign (n={len(benign)}): mean θ-drift = {ben_drift:.4f}")
            print(f"Separation: {path_drift - ben_drift:.4f}")

            # Simple AUC estimate (Mann-Whitney U)
            # Higher drift = pathogenic, so we test path > ben
            from scipy import stats
            path_scores = [m["theta_drift_patch"] for m in pathogenic]
            ben_scores = [m["theta_drift_patch"] for m in benign]
            U, p = stats.mannwhitneyu(path_scores, ben_scores, alternative='greater')
            auc = U / (len(path_scores) * len(ben_scores))
            print(f"Estimated AUC (θ-drift): {auc:.3f} (p={p:.4f})")
        else:
            # Show overall statistics
            drifts = [m["theta_drift_patch"] for m in valid_scores]
            print(f"\n=== θ-drift Statistics ===")
            print(f"Mean: {np.mean(drifts):.4f}")
            print(f"Std: {np.std(drifts):.4f}")
            print(f"Min: {np.min(drifts):.4f}")
            print(f"Max: {np.max(drifts):.4f}")


if __name__ == "__main__":
    main()
