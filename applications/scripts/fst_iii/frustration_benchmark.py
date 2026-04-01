#!/usr/bin/env python3
"""
III-BENCH-01: ROC/AUPRC Benchmark -- Nash Frustration vs. Known Functional Residues.

Evaluates whether Nash frustration scores can distinguish functionally
important residues from structurally passive ones.

Ground truth:
  - HP35 (1YRF): Helix-cap residues are functionally critical
    Functional: 1-3 (N-cap), 10-12 (helix junction), 20-23 (turn), 32-35 (C-cap)
  - Protein G (1PGA): Beta-sheet core vs. loops
    Functional (core): 1-8 (beta1), 41-56 (beta3+beta4)

Outputs: JSON with AUC, AUPRC, optimal threshold; ROC/PR plots as PNG.

Author: Lukas Geiger (with Claude, 2026-03-23)
Part of: FST-III Biological Game Theory
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

import sys
import argparse
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from protein_fold_nash_pdb import (
    AA_TO_IDX,
    load_ca_coords_and_seq, extract_phi_psi, build_contacts,
    fit_params_from_structure, numerical_hessian,
    frustration_map
)


# ---------------------------------------------------------------------------
# Ground truth definitions
# ---------------------------------------------------------------------------

def get_ground_truth_1yrf(n_residues: int) -> np.ndarray:
    """
    HP35 villin headpiece: functional residues based on literature.

    Helix-cap residues are critical for folding nucleation and stability:
    - Residues 1-3: N-terminal cap (helix 1 initiation)
    - Residues 10-12: Helix 1/2 junction (inter-helix contacts)
    - Residues 20-23: Turn region (helix 2/3 linker, hydrophobic core)
    - Residues 32-35: C-terminal cap (helix 3 termination)

    References:
      McKnight et al. (1997) Nat Struct Biol 4, 180-184
      Kubelka et al. (2003) JACS 125, 12397-12407
    """
    labels = np.zeros(n_residues, dtype=int)
    functional_ranges = [(0, 3), (9, 12), (19, 23), (31, 35)]
    for start, end in functional_ranges:
        for i in range(start, min(end, n_residues)):
            labels[i] = 1
    return labels


def get_ground_truth_1pga(n_residues: int) -> np.ndarray:
    """
    Protein G (B1 domain): Beta-sheet core residues as functional ground truth.

    The beta-sheet core is essential for structural stability:
    - Residues 1-8: Beta-strand 1
    - Residues 41-56: Beta-strands 3+4 (C-terminal sheet)

    Loop regions are structurally less constrained.

    References:
      Gronenborn et al. (1991) Science 253, 657-661
      Alexander et al. (2009) Biochemistry 48, 3473-3482
    """
    labels = np.zeros(n_residues, dtype=int)
    functional_ranges = [(0, 8), (40, 56)]
    for start, end in functional_ranges:
        for i in range(start, min(end, n_residues)):
            labels[i] = 1
    return labels


# ---------------------------------------------------------------------------
# ROC / AUPRC computation (no sklearn dependency)
# ---------------------------------------------------------------------------

def compute_roc(labels: np.ndarray, scores: np.ndarray):
    """
    Compute ROC curve and AUC from labels (0/1) and continuous scores.

    Returns:
        fpr: array of false positive rates
        tpr: array of true positive rates
        auc: area under ROC curve
        thresholds: corresponding thresholds
    """
    # Sort by descending score
    order = np.argsort(-scores)
    labels_sorted = labels[order]
    scores_sorted = scores[order]

    n_pos = np.sum(labels == 1)
    n_neg = np.sum(labels == 0)

    if n_pos == 0 or n_neg == 0:
        return np.array([0, 1]), np.array([0, 1]), 0.5, np.array([1.0, 0.0])

    tpr_list = [0.0]
    fpr_list = [0.0]
    thresholds = [scores_sorted[0] + 1.0]

    tp = 0
    fp = 0
    for i in range(len(labels_sorted)):
        if labels_sorted[i] == 1:
            tp += 1
        else:
            fp += 1
        tpr_list.append(tp / n_pos)
        fpr_list.append(fp / n_neg)
        thresholds.append(scores_sorted[i])

    fpr = np.array(fpr_list)
    tpr = np.array(tpr_list)
    thresholds = np.array(thresholds)

    # AUC via trapezoidal rule
    auc = float(np.trapezoid(tpr, fpr))

    return fpr, tpr, auc, thresholds


def compute_precision_recall(labels: np.ndarray, scores: np.ndarray):
    """
    Compute Precision-Recall curve and AUPRC.

    Returns:
        precision: array
        recall: array
        auprc: area under precision-recall curve
    """
    order = np.argsort(-scores)
    labels_sorted = labels[order]

    n_pos = np.sum(labels == 1)
    if n_pos == 0:
        return np.array([0]), np.array([0]), 0.0

    precision_list = []
    recall_list = []
    tp = 0

    for i in range(len(labels_sorted)):
        if labels_sorted[i] == 1:
            tp += 1
        precision_list.append(tp / (i + 1))
        recall_list.append(tp / n_pos)

    precision = np.array(precision_list)
    recall = np.array(recall_list)

    # AUPRC via trapezoidal rule (step function approximation)
    auprc = float(np.trapezoid(precision, recall))

    return precision, recall, auprc


def find_optimal_threshold(labels: np.ndarray, scores: np.ndarray) -> dict:
    """Find threshold that maximizes Youden's J (TPR - FPR)."""
    fpr, tpr, auc, thresholds = compute_roc(labels, scores)
    j_stat = tpr - fpr
    idx_best = np.argmax(j_stat)
    return {
        "threshold": float(thresholds[idx_best]),
        "youden_j": float(j_stat[idx_best]),
        "tpr_at_best": float(tpr[idx_best]),
        "fpr_at_best": float(fpr[idx_best])
    }


# ---------------------------------------------------------------------------
# Benchmark for a single protein
# ---------------------------------------------------------------------------

def benchmark_protein(pdb_path: str, chain_id: str, label: str,
                      ground_truth_fn, lr: float = 0.015,
                      rcut: float = 8.0, sigma: float = 1.0) -> dict:
    """
    Run frustration benchmark for one protein.
    """
    print(f"\n{'='*60}")
    print(f"Benchmark: {label}")
    print(f"{'='*60}")

    # Load structure
    ca, seq = load_ca_coords_and_seq(pdb_path, chain_id)
    phi, psi = extract_phi_psi(pdb_path, chain_id)
    edges, r = build_contacts(ca, r_cut=rcut)
    N = len(seq)

    print(f"  Residues: {N}, Contacts: {len(edges)}")

    # Ground truth
    labels = ground_truth_fn(N)
    n_functional = int(np.sum(labels))
    n_passive = N - n_functional
    print(f"  Functional residues: {n_functional}, Passive: {n_passive}")

    # Fit Nash potential
    print("  Fitting Nash potential...")
    params = fit_params_from_structure(seq, phi, psi, edges, r, sigma=sigma)

    # Compute Hessian and frustration
    print("  Computing Hessian and frustration...")
    aa_idx = np.array([AA_TO_IDX.get(a, 0) for a in seq], dtype=int)
    H = numerical_hessian(phi, psi, aa_idx, edges, r, params, sigma=sigma)

    frust = frustration_map(H, lr, N)

    # ROC analysis
    fpr, tpr, auc, thresholds_roc = compute_roc(labels, frust)
    precision, recall, auprc = compute_precision_recall(labels, frust)
    opt = find_optimal_threshold(labels, frust)

    print(f"\n  Results:")
    print(f"    AUC (ROC):  {auc:.4f}")
    print(f"    AUPRC:      {auprc:.4f}")
    print(f"    Optimal threshold: {opt['threshold']:.4f} (Youden J={opt['youden_j']:.4f})")
    print(f"    TPR at optimal: {opt['tpr_at_best']:.4f}, FPR: {opt['fpr_at_best']:.4f}")

    # Per-residue breakdown
    print(f"\n  Per-residue frustration (functional marked with *):")
    print(f"    {'Idx':>4s} {'AA':>3s} {'Frust':>7s} {'Label':>6s}")
    print(f"    " + "-" * 30)
    for i in range(N):
        marker = " *" if labels[i] == 1 else "  "
        print(f"    {i:4d} {seq[i]:>3s} {frust[i]:7.3f} {marker}")

    # Mean frustration comparison
    func_frust = frust[labels == 1]
    pass_frust = frust[labels == 0]
    print(f"\n  Mean frustration (functional): {np.mean(func_frust):.4f} +/- {np.std(func_frust):.4f}")
    print(f"  Mean frustration (passive):    {np.mean(pass_frust):.4f} +/- {np.std(pass_frust):.4f}")

    return {
        "protein": label,
        "pdb_path": pdb_path,
        "chain": chain_id,
        "n_residues": N,
        "n_functional": n_functional,
        "n_passive": n_passive,
        "auc_roc": auc,
        "auprc": auprc,
        "optimal_threshold": opt,
        "mean_frust_functional": float(np.mean(func_frust)),
        "mean_frust_passive": float(np.mean(pass_frust)),
        "frustration_scores": frust.tolist(),
        "labels": labels.tolist(),
        "roc": {"fpr": fpr.tolist(), "tpr": tpr.tolist()},
        "pr": {"precision": precision.tolist(), "recall": recall.tolist()},
        "sequence": seq
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_benchmark(results: list, output_dir: str):
    """Plot ROC and PR curves for all benchmarked proteins."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping plots.")
        return

    n_prot = len(results)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

    # --- ROC curves ---
    ax1 = axes[0]
    for i, res in enumerate(results):
        fpr = res["roc"]["fpr"]
        tpr = res["roc"]["tpr"]
        auc = res["auc_roc"]
        ax1.plot(fpr, tpr, color=colors[i % len(colors)], linewidth=2,
                 label=f'{res["protein"]} (AUC={auc:.3f})')

    ax1.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random (AUC=0.5)')
    ax1.set_xlabel('False Positive Rate', fontsize=13)
    ax1.set_ylabel('True Positive Rate', fontsize=13)
    ax1.set_title('ROC: Nash Frustration vs. Functional Residues', fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(-0.02, 1.02)
    ax1.set_ylim(-0.02, 1.02)

    # --- Precision-Recall curves ---
    ax2 = axes[1]
    for i, res in enumerate(results):
        precision = res["pr"]["precision"]
        recall = res["pr"]["recall"]
        auprc = res["auprc"]
        prevalence = res["n_functional"] / res["n_residues"]
        ax2.plot(recall, precision, color=colors[i % len(colors)], linewidth=2,
                 label=f'{res["protein"]} (AUPRC={auprc:.3f})')
        # Baseline = prevalence
        ax2.axhline(y=prevalence, color=colors[i % len(colors)],
                     linestyle=':', alpha=0.4)

    ax2.set_xlabel('Recall', fontsize=13)
    ax2.set_ylabel('Precision', fontsize=13)
    ax2.set_title('Precision-Recall: Nash Frustration', fontsize=14)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(-0.02, 1.02)
    ax2.set_ylim(-0.02, 1.05)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, "frustration_benchmark_roc.png")
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nPlot saved: {plot_path}")

    # --- Frustration bar plots ---
    fig2, axes2 = plt.subplots(n_prot, 1, figsize=(12, 4 * n_prot))
    if n_prot == 1:
        axes2 = [axes2]

    for idx, res in enumerate(results):
        ax = axes2[idx]
        frust = np.array(res["frustration_scores"])
        labels = np.array(res["labels"])
        N = len(frust)
        x = np.arange(N)

        bar_colors = ['#d62728' if l == 1 else '#1f77b4' for l in labels]
        ax.bar(x, frust, color=bar_colors, width=1.0, edgecolor='none')
        ax.set_xlabel('Residue Index', fontsize=12)
        ax.set_ylabel('Frustration', fontsize=12)
        ax.set_title(f'{res["protein"]}: Frustration by Residue '
                     f'(red = functional, blue = passive)', fontsize=13)
        ax.grid(True, alpha=0.3, axis='y')

        # Threshold line
        thresh = res["optimal_threshold"]["threshold"]
        ax.axhline(y=thresh, color='green', linestyle='--', alpha=0.7,
                    label=f'Optimal threshold = {thresh:.3f}')
        ax.legend(fontsize=10)

    plt.tight_layout()
    bar_path = os.path.join(output_dir, "frustration_benchmark_bars.png")
    plt.savefig(bar_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Bar plot saved: {bar_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="III-BENCH-01: ROC/AUPRC Benchmark for Nash Frustration"
    )
    ap.add_argument("--data-dir", required=True,
                    help="Directory containing PDB files")
    ap.add_argument("--output-dir", default=None,
                    help="Output directory (default: same as data-dir)")
    ap.add_argument("--lr", type=float, default=0.015,
                    help="Learning rate for Jacobian J = I - lr*H")
    ap.add_argument("--rcut", type=float, default=8.0,
                    help="Contact distance cutoff")
    ap.add_argument("--sigma", type=float, default=1.0,
                    help="RBF sigma")
    args = ap.parse_args()

    if args.output_dir is None:
        args.output_dir = args.data_dir
    os.makedirs(args.output_dir, exist_ok=True)

    # Benchmark definitions: (pdb_file, chain, label, ground_truth_function)
    benchmarks = [
        ("1YRF.pdb", "A", "HP35 (1YRF)", get_ground_truth_1yrf),
        ("1PGA.pdb", "A", "Protein G (1PGA)", get_ground_truth_1pga),
    ]

    all_results = []
    for pdb_file, chain, label, gt_fn in benchmarks:
        pdb_path = os.path.join(args.data_dir, pdb_file)
        if not os.path.exists(pdb_path):
            print(f"WARNING: {pdb_path} not found, skipping {label}")
            continue

        result = benchmark_protein(
            pdb_path, chain, label, gt_fn,
            lr=args.lr, rcut=args.rcut, sigma=args.sigma
        )
        all_results.append(result)

    if not all_results:
        print("ERROR: No PDB files found. Check --data-dir.")
        sys.exit(1)

    # Save JSON
    output = {
        "metadata": {
            "description": "III-BENCH-01: ROC/AUPRC Benchmark Nash Frustration vs Functional Residues",
            "author": "Lukas Geiger",
            "date": "2026-03-23",
            "parameters": {
                "lr": args.lr,
                "rcut": args.rcut,
                "sigma": args.sigma
            },
            "ground_truth_sources": {
                "1YRF": "Helix-cap residues (McKnight 1997, Kubelka 2003)",
                "1PGA": "Beta-sheet core (Gronenborn 1991, Alexander 2009)"
            }
        },
        "benchmarks": all_results
    }

    # Remove large arrays for compact JSON (keep in separate key)
    output_compact = json.loads(json.dumps(output))
    for bench in output_compact["benchmarks"]:
        bench.pop("roc", None)
        bench.pop("pr", None)

    json_path = os.path.join(args.output_dir, "frustration_benchmark_results.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2)
    print(f"\nJSON saved: {json_path}")

    # Plot
    plot_benchmark(all_results, args.output_dir)

    # Summary
    print("\n" + "=" * 60)
    print("BENCHMARK SUMMARY")
    print("=" * 60)
    print(f"\n  {'Protein':>20s} | {'AUC':>6s} | {'AUPRC':>6s} | {'Threshold':>10s} | {'Youden J':>9s}")
    print("  " + "-" * 65)
    for res in all_results:
        opt = res["optimal_threshold"]
        print(f"  {res['protein']:>20s} | {res['auc_roc']:6.4f} | {res['auprc']:6.4f} | "
              f"{opt['threshold']:10.4f} | {opt['youden_j']:9.4f}")

    print(f"\nInterpretation:")
    for res in all_results:
        auc = res['auc_roc']
        if auc > 0.7:
            print(f"  {res['protein']}: GOOD discrimination (AUC={auc:.3f})")
        elif auc > 0.55:
            print(f"  {res['protein']}: MODERATE discrimination (AUC={auc:.3f})")
        else:
            print(f"  {res['protein']}: WEAK discrimination (AUC={auc:.3f})")


if __name__ == "__main__":
    main()
