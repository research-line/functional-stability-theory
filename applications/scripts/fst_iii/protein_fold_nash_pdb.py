#!/usr/bin/env python3
"""
Reverse-engineer a continuous potential-game (Nash) model from Protein A (known structure),
then predict folding (phi,psi) for Protein B using the learned parameters.

State per residue i: theta_i = (phi_i, psi_i)

Potential (sum of two coupled angle-games):
  F = F_phi + F_psi

  F_phi = sum_i k_phi[a_i]*(1 - cos(phi_i - mu_phi[a_i]))
        + sum_(i,j in contacts) w(r_ij)*(1 - cos(phi_i - phi_j))

  F_psi = sum_i k_psi[a_i]*(1 - cos(psi_i - mu_psi[a_i]))
        + sum_(i,j in contacts) w(r_ij)*(1 - cos(psi_i - psi_j))

Reverse engineering:
  Fit mu_phi, mu_psi, k_phi, k_psi, and w(r) (RBF basis) such that grad F(theta*) ~ 0.

Stability:
  Compute numerical Hessian at theta* and report eigenvalue stats.

Prediction:
  Multi-start gradient descent to minimize F for Protein B.

Best-Response Dynamics:
  True game-theoretic dynamics where each residue (player) optimizes only its own
  strategy (phi_i, psi_i) while all others remain fixed.
  Convergence => Nash equilibrium; Divergence/cycles => structural frustration.

Notes:
- This is a deliberately minimal, invertible prototype (not Rosetta/AlphaFold).
- Contacts are derived from CA distances in the *given* PDB (so for B this is a "structure-informed"
  contact graph; if you want true ab initio, replace contacts by predicted contacts).

Reference: FST-III Biological Game Theory (Geiger, 2026)
Implements Prediction P1: Nash equilibrium stability analysis of protein energy landscapes.

Requirements:
    pip install biopython scipy numpy

Usage:
    python protein_fold_nash_pdb.py --pdbA 1YRF.pdb --chainA A --pdbB 1PGA.pdb --chainB A
"""

import argparse
import json
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.vectors import calc_dihedral
from scipy.optimize import least_squares


AA20 = list("ACDEFGHIKLMNPQRSTVWY")
AA_TO_IDX = {a: i for i, a in enumerate(AA20)}


# -----------------------------
# PDB utilities
# -----------------------------

def load_chain_residues(pdb_path: str, chain_id: str):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", pdb_path)
    model = next(structure.get_models())
    chain = model[chain_id]
    residues = [r for r in chain.get_residues() if r.id[0] == " "]
    # keep only residues with backbone atoms
    residues = [r for r in residues if ("N" in r and "CA" in r and "C" in r)]
    return residues


def load_ca_coords_and_seq(pdb_path: str, chain_id: str) -> Tuple[np.ndarray, str]:
    residues = load_chain_residues(pdb_path, chain_id)

    from Bio.Data.IUPACData import protein_letters_3to1_extended as map3

    ca = []
    seq = []
    for res in residues:
        ca.append(res["CA"].coord.astype(float))
        name3 = res.get_resname().strip().upper()
        aa = map3.get(name3.capitalize(), "A")
        aa = aa if aa in AA_TO_IDX else "A"
        seq.append(aa)

    return np.array(ca, dtype=float), "".join(seq)


def extract_phi_psi(pdb_path: str, chain_id: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    phi(i) uses C(i-1)-N(i)-CA(i)-C(i)  (undefined at i=0 -> 0)
    psi(i) uses N(i)-CA(i)-C(i)-N(i+1)  (undefined at i=N-1 -> 0)
    """
    residues = load_chain_residues(pdb_path, chain_id)
    N = len(residues)
    phi = np.zeros(N, dtype=float)
    psi = np.zeros(N, dtype=float)

    for i in range(1, N):
        prev = residues[i - 1]
        cur = residues[i]
        try:
            phi[i] = float(calc_dihedral(prev["C"].get_vector(), cur["N"].get_vector(),
                                         cur["CA"].get_vector(), cur["C"].get_vector()))
        except Exception:
            phi[i] = 0.0

    for i in range(0, N - 1):
        cur = residues[i]
        nxt = residues[i + 1]
        try:
            psi[i] = float(calc_dihedral(cur["N"].get_vector(), cur["CA"].get_vector(),
                                         cur["C"].get_vector(), nxt["N"].get_vector()))
        except Exception:
            psi[i] = 0.0

    return phi, psi


def build_contacts(ca: np.ndarray, r_cut: float = 8.0, seq_sep: int = 3) -> Tuple[np.ndarray, np.ndarray]:
    """
    Contact edges (i,j) for CA distance <= r_cut and |i-j|>=seq_sep.
    Returns edges (M,2) and distances r (M,).
    """
    N = ca.shape[0]
    edges = []
    dists = []
    for i in range(N):
        for j in range(i + seq_sep, N):
            r = np.linalg.norm(ca[i] - ca[j])
            if r <= r_cut:
                edges.append((i, j))
                dists.append(r)
    if not edges:
        edges = [(0, N - 1)]
        dists = [np.linalg.norm(ca[0] - ca[-1])]
    return np.array(edges, dtype=int), np.array(dists, dtype=float)


# -----------------------------
# Model
# -----------------------------

@dataclass
class LearnedParams:
    mu_phi: np.ndarray     # (20,)
    mu_psi: np.ndarray     # (20,)
    k_phi: np.ndarray      # (20,)
    k_psi: np.ndarray      # (20,)
    r_centers: np.ndarray  # (B,)
    w_rbf: np.ndarray      # (B,)
    w0: float              # scalar bias


def rbf_features(r: np.ndarray, centers: np.ndarray, sigma: float) -> np.ndarray:
    return np.exp(-0.5 * ((r[:, None] - centers[None, :]) / sigma) ** 2)


def unpack_params(x: np.ndarray, B: int):
    # layout:
    # 0:20 mu_phi
    # 20:40 mu_psi
    # 40:60 logk_phi
    # 60:80 logk_psi
    # 80 w0
    # 81:81+B w_rbf
    mu_phi = x[0:20]
    mu_psi = x[20:40]
    k_phi = np.exp(x[40:60])
    k_psi = np.exp(x[60:80])
    w0 = x[80]
    w_rbf = x[81:81 + B]
    return mu_phi, mu_psi, k_phi, k_psi, w0, w_rbf


def energy_and_grad(phi: np.ndarray, psi: np.ndarray, aa_idx: np.ndarray,
                    edges: np.ndarray, r: np.ndarray,
                    mu_phi: np.ndarray, mu_psi: np.ndarray,
                    k_phi: np.ndarray, k_psi: np.ndarray,
                    w0: float, w_rbf: np.ndarray,
                    centers: np.ndarray, sigma: float,
                    lambda_hb: float = 0.0) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Returns F, grad_phi, grad_psi

    Includes hydrogen bond coordination term for helix stabilization.
    """
    # unary terms
    mu_i_phi = mu_phi[aa_idx]
    mu_i_psi = mu_psi[aa_idx]
    k_i_phi = k_phi[aa_idx]
    k_i_psi = k_psi[aa_idx]

    F_u = np.sum(k_i_phi * (1.0 - np.cos(phi - mu_i_phi))) + np.sum(k_i_psi * (1.0 - np.cos(psi - mu_i_psi)))
    g_phi = k_i_phi * np.sin(phi - mu_i_phi)
    g_psi = k_i_psi * np.sin(psi - mu_i_psi)

    # pair weights
    Phi = rbf_features(r, centers, sigma)   # (M,B)
    w = w0 + Phi @ w_rbf                    # (M,)

    i_idx = edges[:, 0]
    j_idx = edges[:, 1]

    dphi = phi[i_idx] - phi[j_idx]
    dpsi = psi[i_idx] - psi[j_idx]

    F_p = np.sum(w * (1.0 - np.cos(dphi))) + np.sum(w * (1.0 - np.cos(dpsi)))

    sphi = w * np.sin(dphi)
    spsi = w * np.sin(dpsi)

    np.add.at(g_phi, i_idx, sphi)
    np.add.at(g_phi, j_idx, -sphi)

    np.add.at(g_psi, i_idx, spsi)
    np.add.at(g_psi, j_idx, -spsi)

    # --- Hydrogen bond coordination term ---
    # Bilateral coordination reward: neighbors should have SIMILAR conformations
    # This is Nash-relevant: couples strategies via conformational coherence
    # More general than fixed helix angles - works for any secondary structure
    sigma_hb = 0.8

    # Reward conformational similarity between contact pairs
    # F_hb = -lambda_hb * sum exp(-|theta_i - theta_j|^2 / 2sigma^2)
    dphi_hb = phi[i_idx] - phi[j_idx]
    dpsi_hb = psi[i_idx] - psi[j_idx]

    # Periodic distance (angles wrap around)
    dphi_hb = np.arctan2(np.sin(dphi_hb), np.cos(dphi_hb))
    dpsi_hb = np.arctan2(np.sin(dpsi_hb), np.cos(dpsi_hb))

    coherence = np.exp(-(dphi_hb**2 + dpsi_hb**2) / (2 * sigma_hb**2))
    F_hb = -lambda_hb * np.sum(coherence)

    # Gradients
    d_dphi = coherence * (-dphi_hb / (sigma_hb**2))
    d_dpsi = coherence * (-dpsi_hb / (sigma_hb**2))

    np.add.at(g_phi, i_idx, -lambda_hb * d_dphi)
    np.add.at(g_phi, j_idx, lambda_hb * d_dphi)
    np.add.at(g_psi, i_idx, -lambda_hb * d_dpsi)
    np.add.at(g_psi, j_idx, lambda_hb * d_dpsi)

    return float(F_u + F_p + F_hb), g_phi, g_psi


# -----------------------------
# Reverse engineering
# -----------------------------

def fit_params_from_structure(seq: str, phi_star: np.ndarray, psi_star: np.ndarray,
                              edges: np.ndarray, r: np.ndarray,
                              B: int = 12, sigma: float = 1.0,
                              lam_reg: float = 1e-2) -> LearnedParams:
    aa_idx = np.array([AA_TO_IDX.get(a, 0) for a in seq], dtype=int)

    r_min, r_max = float(np.min(r)), float(np.max(r))
    centers = np.linspace(r_min, r_max, B)

    x0 = np.zeros(81 + B, dtype=float)
    # logk init 0 -> k=1
    x0[40:60] = 0.0
    x0[60:80] = 0.0
    x0[80] = 0.1  # w0

    def residual(x: np.ndarray) -> np.ndarray:
        mu_phi, mu_psi, k_phi, k_psi, w0, w_rbf = unpack_params(x, B)
        _, gphi, gpsi = energy_and_grad(phi_star, psi_star, aa_idx, edges, r,
                                        mu_phi, mu_psi, k_phi, k_psi,
                                        w0, w_rbf, centers, sigma)
        g = np.concatenate([gphi, gpsi])
        reg = math.sqrt(lam_reg) * x
        return np.concatenate([g, reg])

    res = least_squares(residual, x0, method="trf", max_nfev=3000)

    mu_phi, mu_psi, k_phi, k_psi, w0, w_rbf = unpack_params(res.x, B)
    return LearnedParams(mu_phi=mu_phi, mu_psi=mu_psi, k_phi=k_phi, k_psi=k_psi,
                         r_centers=centers, w_rbf=w_rbf, w0=w0)


# -----------------------------
# Stability diagnostics
# -----------------------------

def numerical_hessian(phi: np.ndarray, psi: np.ndarray, aa_idx: np.ndarray,
                      edges: np.ndarray, r: np.ndarray, params: LearnedParams,
                      sigma: float = 1.0, eps: float = 1e-4) -> np.ndarray:
    """
    Numerical Hessian of F wrt x = [phi, psi] at (phi,psi).
    Central differences on gradient.
    """
    N = len(phi)
    x0 = np.concatenate([phi, psi])
    H = np.zeros((2 * N, 2 * N), dtype=float)

    def grad(x):
        ph = x[:N]
        ps = x[N:]
        _, gph, gps = energy_and_grad(ph, ps, aa_idx, edges, r,
                                      params.mu_phi, params.mu_psi,
                                      params.k_phi, params.k_psi,
                                      params.w0, params.w_rbf,
                                      params.r_centers, sigma)
        return np.concatenate([gph, gps])

    g0 = grad(x0)
    for j in range(2 * N):
        dx = np.zeros_like(x0)
        dx[j] = eps
        gp = grad(x0 + dx)
        gm = grad(x0 - dx)
        H[:, j] = (gp - gm) / (2.0 * eps)

    # symmetrize (numerical noise)
    H = 0.5 * (H + H.T)
    return H


# -----------------------------
# Prediction
# -----------------------------

def fold_protein(seq: str, edges: np.ndarray, r: np.ndarray, params: LearnedParams,
                 sigma: float = 1.0, n_starts: int = 64, steps: int = 6000, lr: float = 0.015,
                 seed: int = 0) -> Dict:
    rng = np.random.default_rng(seed)
    aa_idx = np.array([AA_TO_IDX.get(a, 0) for a in seq], dtype=int)
    N = len(seq)

    best = None
    for s in range(n_starts):
        phi = rng.uniform(-np.pi, np.pi, size=N)
        psi = rng.uniform(-np.pi, np.pi, size=N)

        for _ in range(steps):
            F, gphi, gpsi = energy_and_grad(phi, psi, aa_idx, edges, r,
                                            params.mu_phi, params.mu_psi,
                                            params.k_phi, params.k_psi,
                                            params.w0, params.w_rbf,
                                            params.r_centers, sigma)
            phi -= lr * gphi
            psi -= lr * gpsi

        F, gphi, gpsi = energy_and_grad(phi, psi, aa_idx, edges, r,
                                        params.mu_phi, params.mu_psi,
                                        params.k_phi, params.k_psi,
                                        params.w0, params.w_rbf,
                                        params.r_centers, sigma)

        if best is None or F < best["F"]:
            # wrap to [-pi,pi]
            phi_w = (phi + np.pi) % (2 * np.pi) - np.pi
            psi_w = (psi + np.pi) % (2 * np.pi) - np.pi
            best = {
                "F": float(F),
                "grad_norm": float(np.linalg.norm(np.concatenate([gphi, gpsi]))),
                "phi": phi_w.tolist(),
                "psi": psi_w.tolist(),
            }

    return best


# -----------------------------
# Best-Response Dynamics (True Nash)
# -----------------------------

def best_response_dynamics(seq: str, edges: np.ndarray, r: np.ndarray, params: LearnedParams,
                           phi_init: np.ndarray, psi_init: np.ndarray,
                           sigma: float = 1.0, lr: float = 0.05,
                           sweeps: int = 200, inner_steps: int = 10,
                           tol: float = 1e-4, verbose: bool = False) -> Dict:
    """
    Continuous Nash best-response dynamics.

    Each sweep: every residue i updates (phi_i, psi_i) via local gradient descent
    while ALL other residues remain fixed. This is the game-theoretic definition
    of best-response dynamics.

    Convergence => Nash equilibrium (no player benefits from unilateral deviation)
    Divergence/cycles => structural frustration, no stable Nash

    Args:
        seq: Amino acid sequence
        edges: Contact edges (M, 2)
        r: Contact distances (M,)
        params: Learned potential parameters
        phi_init, psi_init: Initial angles
        sigma: RBF sigma for distance weighting
        lr: Learning rate for local optimization
        sweeps: Number of full passes over all residues
        inner_steps: Steps of local gradient descent per residue per sweep
        tol: Convergence tolerance on max angle update
        verbose: Print progress

    Returns:
        Dict with final angles, convergence info, and Nash stability metrics
    """
    aa_idx = np.array([AA_TO_IDX.get(a, 0) for a in seq], dtype=int)
    N = len(seq)

    phi = phi_init.copy()
    psi = psi_init.copy()

    history = []
    converged = False

    for sweep in range(sweeps):
        max_update = 0.0
        sweep_energy = 0.0

        for i in range(N):
            # Store old values
            phi_old = phi[i]
            psi_old = psi[i]

            # Local optimization: only residue i moves
            for _ in range(inner_steps):
                _, gphi, gpsi = energy_and_grad(
                    phi, psi, aa_idx, edges, r,
                    params.mu_phi, params.mu_psi,
                    params.k_phi, params.k_psi,
                    params.w0, params.w_rbf,
                    params.r_centers, sigma
                )

                # Unilateral update (only residue i)
                phi[i] -= lr * gphi[i]
                psi[i] -= lr * gpsi[i]

            # Track max update
            dphi = abs(phi[i] - phi_old)
            dpsi = abs(psi[i] - psi_old)
            max_update = max(max_update, dphi, dpsi)

        # Wrap angles to [-pi, pi]
        phi = (phi + np.pi) % (2 * np.pi) - np.pi
        psi = (psi + np.pi) % (2 * np.pi) - np.pi

        # Compute energy after sweep
        F, gphi, gpsi = energy_and_grad(
            phi, psi, aa_idx, edges, r,
            params.mu_phi, params.mu_psi,
            params.k_phi, params.k_psi,
            params.w0, params.w_rbf,
            params.r_centers, sigma
        )
        grad_norm = float(np.linalg.norm(np.concatenate([gphi, gpsi])))

        history.append({
            "sweep": sweep,
            "F": float(F),
            "grad_norm": grad_norm,
            "max_update": float(max_update)
        })

        if verbose and sweep % 20 == 0:
            print(f"  Sweep {sweep}: F={F:.4f}, |grad|={grad_norm:.4e}, max_upd={max_update:.4e}")

        if max_update < tol:
            converged = True
            if verbose:
                print(f"  Converged after {sweep+1} sweeps (max_update < {tol})")
            break

    # Final state
    F_final, gphi_final, gpsi_final = energy_and_grad(
        phi, psi, aa_idx, edges, r,
        params.mu_phi, params.mu_psi,
        params.k_phi, params.k_psi,
        params.w0, params.w_rbf,
        params.r_centers, sigma
    )

    return {
        "phi": phi,
        "psi": psi,
        "converged": converged,
        "sweeps": len(history),
        "F": float(F_final),
        "grad_norm": float(np.linalg.norm(np.concatenate([gphi_final, gpsi_final]))),
        "history": history
    }


def validate_nash_stability(seq: str, edges: np.ndarray, r: np.ndarray, params: LearnedParams,
                            phi_native: np.ndarray, psi_native: np.ndarray,
                            sigma: float = 1.0, perturbation: float = 0.1,
                            sweeps: int = 300, seed: int = 42, verbose: bool = True) -> Dict:
    """
    Validate that native structure is Nash-stable under best-response dynamics.

    Perturb native angles slightly, then run best-response dynamics.
    If dynamics converge back to native => Nash stable.

    Args:
        phi_native, psi_native: Native (experimental) angles
        perturbation: Std of Gaussian perturbation to add
        ...

    Returns:
        Dict with RMS deviations and convergence info
    """
    rng = np.random.default_rng(seed)

    # Perturb native structure
    phi_init = phi_native + perturbation * rng.standard_normal(len(phi_native))
    psi_init = psi_native + perturbation * rng.standard_normal(len(psi_native))

    if verbose:
        print(f"  Initial perturbation: phi_rms={np.sqrt(np.mean((phi_init - phi_native)**2)):.4f}, "
              f"psi_rms={np.sqrt(np.mean((psi_init - psi_native)**2)):.4f}")

    # Run best-response dynamics
    result = best_response_dynamics(
        seq, edges, r, params,
        phi_init, psi_init,
        sigma=sigma, sweeps=sweeps, verbose=verbose
    )

    # Compare to native
    phi_br = result["phi"]
    psi_br = result["psi"]

    # Angular difference (accounting for periodicity)
    def angular_diff(a, b):
        d = a - b
        return (d + np.pi) % (2 * np.pi) - np.pi

    dphi = angular_diff(phi_br, phi_native)
    dpsi = angular_diff(psi_br, psi_native)

    phi_rms = float(np.sqrt(np.mean(dphi**2)))
    psi_rms = float(np.sqrt(np.mean(dpsi**2)))

    result["phi_rms_vs_native"] = phi_rms
    result["psi_rms_vs_native"] = psi_rms
    result["is_nash_stable"] = (phi_rms < 0.3) and (psi_rms < 0.3) and result["converged"]

    if verbose:
        print(f"  Final: phi_rms={phi_rms:.4f}, psi_rms={psi_rms:.4f}")
        print(f"  Nash stable: {result['is_nash_stable']}")

    return result


# -----------------------------
# Frustration Analysis (Nash Instability per Residue)
# -----------------------------

def frustration_map(H: np.ndarray, lr: float, N: int) -> np.ndarray:
    """
    Compute per-residue Nash frustration from unstable eigenmodes
    of the best-response Jacobian J = I - lr * H.

    This identifies residues that prevent Nash stability - where
    game-theoretic coordination breaks down.

    Args:
        H: Hessian matrix (2N x 2N) at the fixed point
        lr: Learning rate used in best-response dynamics
        N: Number of residues

    Returns:
        Array of frustration scores (N,), normalized to [0, 1]
    """
    J = np.eye(H.shape[0]) - lr * H
    eigvals, eigvecs = np.linalg.eig(J)

    frust = np.zeros(N, dtype=float)

    for k, lam in enumerate(eigvals):
        if abs(lam) >= 1.0:
            weight = abs(lam) - 1.0
            v = eigvecs[:, k].real  # (2N,)
            for i in range(N):
                # Contributions from phi_i and psi_i components
                frust[i] += weight * (v[i]**2 + v[i + N]**2)

    # Normalize for readability
    if frust.max() > 0:
        frust /= frust.max()

    return frust


def analyze_frustration(seq: str, H: np.ndarray, lr: float, verbose: bool = True) -> Dict:
    """
    Full frustration analysis with statistics and residue breakdown.

    Returns:
        Dict with frustration scores, statistics, and interpretation
    """
    N = len(seq)
    frust = frustration_map(H, lr, N)

    # Jacobian eigenvalue analysis
    J = np.eye(H.shape[0]) - lr * H
    eigvals = np.linalg.eigvals(J)
    unstable_count = np.sum(np.abs(eigvals) >= 1.0)
    spectral_radius = float(np.max(np.abs(eigvals)))

    # Classify residues
    high_frust = [(i, seq[i], f) for i, f in enumerate(frust) if f > 0.5]
    low_frust = [(i, seq[i], f) for i, f in enumerate(frust) if f < 0.1]

    result = {
        "frustration_scores": frust.tolist(),
        "sequence": seq,
        "mean_frustration": float(np.mean(frust)),
        "max_frustration": float(np.max(frust)),
        "unstable_modes": int(unstable_count),
        "spectral_radius": spectral_radius,
        "high_frustration_residues": [(i, aa, float(f)) for i, aa, f in high_frust],
        "low_frustration_residues": [(i, aa, float(f)) for i, aa, f in low_frust],
    }

    if verbose:
        print("\n=== Nash Frustration Analysis ===")
        print(f"Spectral radius of J: {spectral_radius:.4f} (stable if < 1)")
        print(f"Unstable modes: {unstable_count}/{2*N}")
        print(f"Mean frustration: {np.mean(frust):.4f}")
        print(f"\nPer-residue frustration:")
        print("  Idx  AA  Frustration  Bar")
        print("  " + "-"*40)
        for i, (aa, f) in enumerate(zip(seq, frust)):
            bar = "#" * int(f * 20)
            marker = " <-- HIGH" if f > 0.5 else ""
            print(f"  {i:3d}  {aa}   {f:.3f}       |{bar:<20}|{marker}")

        if high_frust:
            print(f"\nHigh-frustration residues (>0.5): {[f'{aa}{i}' for i, aa, _ in high_frust]}")
        if low_frust:
            print(f"Low-frustration anchors (<0.1): {[f'{aa}{i}' for i, aa, _ in low_frust]}")

    return result


def write_frustration_pdb(pdb_path: str, chain_id: str, frustration: np.ndarray,
                          out_path: str, scale: float = 100.0) -> None:
    """
    Write PDB with frustration scores as B-factors for 3D visualization.

    The output PDB can be visualized in PyMOL or ChimeraX with:
    - PyMOL: spectrum b, blue_white_red, minimum=0, maximum=100
    - ChimeraX: color bfactor palette blue:white:red range 0,100

    Args:
        pdb_path: Input PDB file
        chain_id: Chain to modify
        frustration: Per-residue frustration scores (0-1 range)
        out_path: Output PDB filename
        scale: Multiply scores for B-factor (default 100 -> 0-100 range)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", pdb_path)
    model = next(structure.get_models())
    chain = model[chain_id]

    # Get residues with CA atoms (standard residues)
    residues = [r for r in chain.get_residues()
                if r.id[0] == " " and "CA" in r]

    if len(residues) != len(frustration):
        print(f"Warning: {len(residues)} residues vs {len(frustration)} frustration scores")
        # Truncate to minimum
        n = min(len(residues), len(frustration))
        residues = residues[:n]
        frustration = frustration[:n]

    # Set B-factor for all atoms in each residue
    for res, frust in zip(residues, frustration):
        b = float(frust * scale)
        for atom in res.get_atoms():
            atom.set_bfactor(b)

    # Write modified PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path)
    print(f"Wrote frustration PDB: {out_path}")
    print(f"  Visualize in PyMOL: spectrum b, blue_white_red, minimum=0, maximum={int(scale)}")
    print(f"  Visualize in ChimeraX: color bfactor palette blue:white:red range 0,{int(scale)}")


def predict_with_best_response(seq: str, edges: np.ndarray, r: np.ndarray, params: LearnedParams,
                               sigma: float = 1.0, n_starts: int = 16,
                               sweeps: int = 500, seed: int = 0,
                               verbose: bool = True) -> Dict:
    """
    Predict protein structure using best-response dynamics from random starts.

    Args:
        n_starts: Number of random initializations
        sweeps: Max sweeps per run
        ...

    Returns:
        Best result (lowest energy) among all runs
    """
    rng = np.random.default_rng(seed)
    N = len(seq)

    best = None
    all_results = []

    for s in range(n_starts):
        phi_init = rng.uniform(-np.pi, np.pi, size=N)
        psi_init = rng.uniform(-np.pi, np.pi, size=N)

        if verbose:
            print(f"  Start {s+1}/{n_starts}...")

        result = best_response_dynamics(
            seq, edges, r, params,
            phi_init, psi_init,
            sigma=sigma, sweeps=sweeps, verbose=False
        )

        all_results.append({
            "start": s,
            "F": result["F"],
            "converged": result["converged"],
            "sweeps": result["sweeps"],
            "grad_norm": result["grad_norm"]
        })

        if best is None or result["F"] < best["F"]:
            best = result.copy()
            best["start_idx"] = s

    if verbose:
        print(f"  Best: start={best['start_idx']}, F={best['F']:.4f}, converged={best['converged']}")

    best["all_starts"] = all_results
    return best


# -----------------------------
# CLI
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Protein Folding as Nash Equilibrium with Best-Response Dynamics"
    )
    ap.add_argument("--pdbA", required=True, help="PDB file for protein A (training), e.g. 1YRF.pdb")
    ap.add_argument("--chainA", required=True, help="Chain ID for protein A, e.g. A")
    ap.add_argument("--pdbB", required=True, help="PDB file for protein B (prediction), e.g. 1PGA.pdb")
    ap.add_argument("--chainB", required=True, help="Chain ID for protein B, e.g. A")
    ap.add_argument("--rcut", type=float, default=8.0)
    ap.add_argument("--B", type=int, default=12)
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--lam", type=float, default=1e-2)
    ap.add_argument("--starts", type=int, default=64)
    ap.add_argument("--steps", type=int, default=6000)
    ap.add_argument("--lr", type=float, default=0.015)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--hess_eps", type=float, default=1e-4)
    ap.add_argument("--out", default="nash_phi_psi_fit_and_fold.json")
    # Best-response dynamics options
    ap.add_argument("--best-response", action="store_true",
                    help="Use best-response dynamics instead of global gradient descent")
    ap.add_argument("--br-sweeps", type=int, default=300,
                    help="Number of sweeps for best-response dynamics")
    ap.add_argument("--br-starts", type=int, default=16,
                    help="Number of random starts for best-response prediction")
    ap.add_argument("--validate-nash", action="store_true",
                    help="Validate that native structure A is Nash-stable")
    ap.add_argument("--write-frustration-pdb", action="store_true",
                    help="Write PDB with frustration as B-factor for visualization")
    ap.add_argument("--lambda-hb", type=float, default=0.0,
                    help="H-bond coordination strength (default: 0.0, set >0 to enable)")
    args = ap.parse_args()

    # --- A: reverse engineer ---
    caA, seqA = load_ca_coords_and_seq(args.pdbA, args.chainA)
    phiA, psiA = extract_phi_psi(args.pdbA, args.chainA)
    edgesA, rA = build_contacts(caA, r_cut=args.rcut)

    params = fit_params_from_structure(seqA, phiA, psiA, edgesA, rA,
                                       B=args.B, sigma=args.sigma, lam_reg=args.lam)

    aaA = np.array([AA_TO_IDX.get(a, 0) for a in seqA], dtype=int)
    FA, gphiA, gpsiA = energy_and_grad(phiA, psiA, aaA, edgesA, rA,
                                       params.mu_phi, params.mu_psi,
                                       params.k_phi, params.k_psi,
                                       params.w0, params.w_rbf,
                                       params.r_centers, args.sigma)
    grad_norm_A = float(np.linalg.norm(np.concatenate([gphiA, gpsiA])))

    # Hessian stability at theta*
    H = numerical_hessian(phiA, psiA, aaA, edgesA, rA, params, sigma=args.sigma, eps=args.hess_eps)
    evals = np.linalg.eigvalsh(H)
    pos_frac = float(np.mean(evals > 1e-8))
    evals_min = float(np.min(evals))
    evals_med = float(np.median(evals))
    evals_max = float(np.max(evals))

    # Frustration analysis - identify residues that prevent Nash stability
    print("\n=== Frustration Analysis (Protein A) ===")
    frustration_A = analyze_frustration(seqA, H, lr=args.lr, verbose=True)

    # Write frustration as B-factor PDB for 3D visualization
    if args.write_frustration_pdb:
        frust_scores = np.array(frustration_A["frustration_scores"])
        out_pdb = args.pdbA.replace(".pdb", "_frustration.pdb")
        write_frustration_pdb(args.pdbA, args.chainA, frust_scores, out_pdb)

    # --- Nash stability validation on A ---
    nash_validation_A = None
    if args.validate_nash or args.best_response:
        print("\n=== Nash Stability Validation (Protein A) ===")
        nash_validation_A = validate_nash_stability(
            seqA, edgesA, rA, params,
            phiA, psiA,
            sigma=args.sigma,
            perturbation=0.1,
            sweeps=args.br_sweeps,
            seed=args.seed,
            verbose=True
        )

    # --- B: predict ---
    caB, seqB = load_ca_coords_and_seq(args.pdbB, args.chainB)
    phiB_native, psiB_native = extract_phi_psi(args.pdbB, args.chainB)
    edgesB, rB = build_contacts(caB, r_cut=args.rcut)

    if args.best_response:
        print("\n=== Best-Response Prediction (Protein B) ===")
        predB_br = predict_with_best_response(
            seqB, edgesB, rB, params,
            sigma=args.sigma,
            n_starts=args.br_starts,
            sweeps=args.br_sweeps,
            seed=args.seed,
            verbose=True
        )

        # Compare to native B
        def angular_diff(a, b):
            d = a - b
            return (d + np.pi) % (2 * np.pi) - np.pi

        dphi = angular_diff(predB_br["phi"], phiB_native)
        dpsi = angular_diff(predB_br["psi"], psiB_native)
        predB_br["phi_rms_vs_native"] = float(np.sqrt(np.mean(dphi**2)))
        predB_br["psi_rms_vs_native"] = float(np.sqrt(np.mean(dpsi**2)))

        print(f"  Comparison to native B: phi_rms={predB_br['phi_rms_vs_native']:.4f}, "
              f"psi_rms={predB_br['psi_rms_vs_native']:.4f}")

        # Convert numpy arrays to lists for JSON
        predB = {
            "method": "best_response",
            "F": predB_br["F"],
            "grad_norm": predB_br["grad_norm"],
            "converged": predB_br["converged"],
            "sweeps": predB_br["sweeps"],
            "phi": predB_br["phi"].tolist(),
            "psi": predB_br["psi"].tolist(),
            "phi_rms_vs_native": predB_br["phi_rms_vs_native"],
            "psi_rms_vs_native": predB_br["psi_rms_vs_native"],
            "n_starts": args.br_starts,
        }
    else:
        predB = fold_protein(seqB, edgesB, rB, params,
                             sigma=args.sigma, n_starts=args.starts,
                             steps=args.steps, lr=args.lr, seed=args.seed)
        predB["method"] = "global_gradient_descent"

    result = {
        "trainA": {
            "pdb": args.pdbA,
            "chain": args.chainA,
            "seq_len": len(seqA),
            "contacts": int(len(edgesA)),
            "F(thetaA)": FA,
            "grad_norm(thetaA)": grad_norm_A,
            "hessian": {
                "eps": args.hess_eps,
                "pos_frac": pos_frac,
                "min": evals_min,
                "median": evals_med,
                "max": evals_max,
            },
            "frustration": frustration_A,
        },
        "learned_params": {
            "AA20": AA20,
            "mu_phi": params.mu_phi.tolist(),
            "mu_psi": params.mu_psi.tolist(),
            "k_phi": params.k_phi.tolist(),
            "k_psi": params.k_psi.tolist(),
            "w0": float(params.w0),
            "r_centers": params.r_centers.tolist(),
            "w_rbf": params.w_rbf.tolist(),
            "sigma": args.sigma,
            "lam": args.lam,
        },
        "predictB": {
            "pdb": args.pdbB,
            "chain": args.chainB,
            "seq_len": len(seqB),
            "contacts": int(len(edgesB)),
            "best": predB,
        },
    }

    # Add Nash validation results if computed
    if nash_validation_A is not None:
        result["nash_validation_A"] = {
            "converged": nash_validation_A["converged"],
            "sweeps": nash_validation_A["sweeps"],
            "phi_rms_vs_native": nash_validation_A["phi_rms_vs_native"],
            "psi_rms_vs_native": nash_validation_A["psi_rms_vs_native"],
            "is_nash_stable": nash_validation_A["is_nash_stable"],
            "F_final": nash_validation_A["F"],
            "grad_norm_final": nash_validation_A["grad_norm"],
        }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    print(f"\nWrote: {args.out}")
    print("\n=== Summary ===")
    print(f"A: grad_norm(theta*) = {grad_norm_A:.6f}")
    print(f"A: Hessian pos_frac = {pos_frac:.2f}, min/med/max = {evals_min:.4f} / {evals_med:.4f} / {evals_max:.4f}")

    if nash_validation_A is not None:
        print(f"A: Nash stable = {nash_validation_A['is_nash_stable']} "
              f"(phi_rms={nash_validation_A['phi_rms_vs_native']:.4f}, "
              f"psi_rms={nash_validation_A['psi_rms_vs_native']:.4f})")

    print(f"B: best F = {predB['F']:.4f}, grad_norm = {predB['grad_norm']:.6f}")
    if "phi_rms_vs_native" in predB:
        print(f"B: vs native: phi_rms={predB['phi_rms_vs_native']:.4f}, psi_rms={predB['psi_rms_vs_native']:.4f}")


if __name__ == "__main__":
    main()
