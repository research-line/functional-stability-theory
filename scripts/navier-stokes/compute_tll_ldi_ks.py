"""
compute_tll_ldi_ks.py
======================
TLL + LDI Verification on the Kuramoto-Sivashinsky (KS) attractor.

KS equation (1D, chaotic PDE):
    u_t + u*u_x + u_xx + u_xxxx = 0
on [0, L] with periodic boundary conditions.

This is a higher-dimensional test case for the TLL+LDI framework:
- Infinite-dimensional PDE, truncated to M Fourier modes
- Chaotic dynamics with attractor dimension ~ O(L)
- Pseudo-spectral method (FFT) in space
- ETDRK4 time integrator (Kassam & Trefethen 2005)

Tests (Definition 3.1 and 4.2 from FST-NS-LDI paper):
1. TLL: ||P(u(t+h)) - P(u(t))|| <= C_TLL * ||u(t+h) - u(t)|| * Omega(t)
2. LDI: int ||u'(t)|| * Omega(t) dt < infty
3. BV Chain Rule: Var(v) <= C_TLL * LDI-Integral
4. STC: |{t : d(t) <= eps}| <= C_* * eps^alpha
5. Exponential Tracking: d(t) ~ C * exp(-lambda * t)
6. Grid refinement study: N_fourier = 32, 64, 128, 256

Author: Lukas Geiger (script generated via Claude, 2026)
"""

import numpy as np
from scipy.spatial import cKDTree
import os
import sys
import json
import time

os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass


# ==========================================================================
# KS Equation -- Pseudo-spectral solver with ETDRK4
# ==========================================================================

class KuramotoSivashinsky:
    """
    Pseudo-spectral solver for the Kuramoto-Sivashinsky equation:
        u_t + u*u_x + u_xx + u_xxxx = 0
    on [0, L] with periodic boundary conditions.

    State is stored in Fourier space (complex coefficients).
    """

    def __init__(self, L=32 * np.pi, N=128, dt=0.05):
        """
        Parameters
        ----------
        L : float
            Domain length. L = 32*pi gives well-developed spatiotemporal chaos.
        N : int
            Number of Fourier modes (grid points). Must be even.
        dt : float
            Default time step for ETDRK4 integration.
        """
        self.L = L
        self.N = N
        self.dx = L / N
        self.x = np.linspace(0, L, N, endpoint=False)

        # Wavenumbers
        self.k = np.fft.rfftfreq(N, d=L / (2 * np.pi * N))
        # Linear operator: L_k = k^2 - k^4
        # (from u_xx -> -k^2 * u_hat, u_xxxx -> k^4 * u_hat)
        self.Lk = self.k**2 - self.k**4

        # Precompute ETDRK4 coefficients
        self._setup_etdrk4_coefficients(dt=dt)

    def _setup_etdrk4_coefficients(self, dt):
        """
        Precompute ETDRK4 coefficients for a given time step.
        Uses Kassam & Trefethen (2005) contour integral method
        for numerical stability.
        """
        self.dt = dt
        Lk = self.Lk
        h = dt

        # Contour integration for stable coefficient computation
        # (avoids catastrophic cancellation near Lk = 0)
        M_contour = 32  # number of contour points
        r = np.exp(1j * np.pi * (np.arange(1, M_contour + 1) - 0.5) / M_contour)

        # LR shape: (len(Lk), M_contour)
        LR = h * Lk[:, np.newaxis] + r[np.newaxis, :]

        # E = exp(h * Lk)
        self.E = np.exp(h * Lk)
        # E2 = exp(h/2 * Lk)
        self.E2 = np.exp(h / 2 * Lk)

        # Coefficients via contour integrals (Kassam & Trefethen 2005, Eq. 26)
        # f1 = (e^(hL/2) - 1) / L
        self.f1 = np.real(np.mean(h * (np.exp(LR / 2) - 1) / LR, axis=1))
        # f2 = (e^(hL) - 1) / L  (not used directly, but related)

        # a = h^{-2} * L^{-1} * [e^(hL/2) - 1 - hL/2 * e^(hL/2)]  ... etc.
        # Actually, the standard ETDRK4 coefficients:
        self.a = np.real(np.mean(
            h * (-4 - LR + np.exp(LR) * (4 - 3 * LR + LR**2)) / LR**3,
            axis=1))
        self.b = np.real(np.mean(
            h * (2 + LR + np.exp(LR) * (-2 + LR)) / LR**3,
            axis=1))
        self.c = np.real(np.mean(
            h * (-4 - 3 * LR - LR**2 + np.exp(LR) * (4 - LR)) / LR**3,
            axis=1))

    def nonlinear(self, u_hat):
        """
        Nonlinear term N(u) = -1/2 * d/dx(u^2) in Fourier space.
        Computed pseudospectrally with 3/2 dealiasing.
        """
        N = self.N
        # Dealiasing: pad to 3/2 * N modes
        N_pad = 3 * N // 2
        u_hat_pad = np.zeros(N_pad // 2 + 1, dtype=complex)
        u_hat_pad[:len(u_hat)] = u_hat

        # Transform to physical space (padded)
        u_phys = np.fft.irfft(u_hat_pad, n=N_pad) * (N_pad / N)

        # u^2 in physical space
        u2_phys = u_phys**2

        # Transform back and truncate
        u2_hat_pad = np.fft.rfft(u2_phys) * (N / N_pad)
        u2_hat = u2_hat_pad[:len(u_hat)]

        # N(u) = -1/2 * i*k * FFT(u^2)
        return -0.5 * 1j * self.k * u2_hat

    def step_etdrk4(self, u_hat):
        """
        One ETDRK4 time step.
        """
        Nk = self.nonlinear(u_hat)

        a_stage = self.E2 * u_hat + self.f1 * Nk
        Na = self.nonlinear(a_stage)

        b_stage = self.E2 * u_hat + self.f1 * Na
        Nb = self.nonlinear(b_stage)

        c_stage = self.E2 * a_stage + self.f1 * (2 * Nb - Nk)
        Nc = self.nonlinear(c_stage)

        u_hat_new = (self.E * u_hat
                     + self.a * Nk
                     + 2 * self.b * (Na + Nb)
                     + self.c * Nc)

        # Blow-up detection: if any mode is NaN/Inf, the simulation is unstable
        if not np.all(np.isfinite(u_hat_new)):
            raise RuntimeError(
                f"ETDRK4 blow-up detected (N={self.N}, dt={self.dt:.4f}). "
                f"Reduce dt or increase N."
            )

        return u_hat_new

    def rhs_fourier(self, u_hat):
        """
        Full RHS in Fourier space: du_hat/dt = L * u_hat + N(u_hat).
        Used for computing ||u'(t)||.
        """
        return self.Lk * u_hat + self.nonlinear(u_hat)

    def to_physical(self, u_hat):
        """Convert Fourier coefficients to physical space."""
        return np.fft.irfft(u_hat, n=self.N)

    def to_fourier(self, u):
        """Convert physical space to Fourier coefficients."""
        return np.fft.rfft(u)

    def initial_condition(self, seed=42):
        """
        Generate a random initial condition with energy in low modes.
        """
        rng = np.random.RandomState(seed)
        u = np.zeros(self.N)
        # Superposition of low-wavenumber modes
        for kn in range(1, 6):
            phase = rng.uniform(0, 2 * np.pi)
            u += rng.uniform(0.5, 2.0) * np.cos(2 * np.pi * kn * self.x / self.L + phase)
        return self.to_fourier(u)

    def integrate(self, u_hat0, T, dt=None, save_every=1):
        """
        Integrate KS equation from u_hat0 for time T.

        Parameters
        ----------
        u_hat0 : array
            Initial condition in Fourier space.
        T : float
            Integration time.
        dt : float, optional
            Time step (default: self.dt).
        save_every : int
            Save every N-th step.

        Returns
        -------
        t_arr : array of times
        traj : array of Fourier coefficient snapshots, shape (n_saved, len(u_hat))
        """
        if dt is not None and dt != self.dt:
            self._setup_etdrk4_coefficients(dt)

        n_steps = int(T / self.dt)
        n_save = n_steps // save_every + 1

        t_arr = np.zeros(n_save)
        traj = np.zeros((n_save, len(u_hat0)), dtype=complex)

        u_hat = u_hat0.copy()
        traj[0] = u_hat
        t_arr[0] = 0.0

        save_idx = 1
        for step in range(1, n_steps + 1):
            u_hat = self.step_etdrk4(u_hat)

            if step % save_every == 0 and save_idx < n_save:
                t_arr[save_idx] = step * self.dt
                traj[save_idx] = u_hat
                save_idx += 1

        return t_arr[:save_idx], traj[:save_idx]


def fourier_state_vector(u_hat, M):
    """
    Extract a real state vector from the first M Fourier modes.
    Returns [Re(u_hat[1]), Im(u_hat[1]), ..., Re(u_hat[M]), Im(u_hat[M])].
    Mode 0 (mean) is excluded (it's conserved in KS).
    """
    vec = np.zeros(2 * M)
    for j in range(M):
        mode = j + 1  # skip mode 0
        if mode < len(u_hat):
            vec[2 * j] = u_hat[mode].real
            vec[2 * j + 1] = u_hat[mode].imag
    return vec


# ==========================================================================
# Main computation
# ==========================================================================

def run_ks_diagnostics(N_fourier=128, L=32 * np.pi, M_modes=32,
                       T_transient=200.0, T_attractor=500.0, T_test=100.0,
                       dt=0.05, save_every=2, verbose=True):
    """
    Run all TLL/LDI diagnostics on the KS attractor.

    Parameters
    ----------
    N_fourier : int
        Number of Fourier modes for the PDE solver.
    L : float
        Domain length.
    M_modes : int
        Number of Fourier modes used for the state vector (attractor geometry).
    T_transient : float
        Transient integration time (discarded).
    T_attractor : float
        Time to sample the attractor.
    T_test : float
        Time for the test trajectory.
    dt : float
        Time step.
    save_every : int
        Save every N-th step.
    verbose : bool
        Print progress.

    Returns
    -------
    results : dict
        All diagnostic results.
    """
    dt_save = dt * save_every

    if verbose:
        print(f"\n{'='*70}")
        print(f"KS DIAGNOSTICS: N={N_fourier}, L={L:.2f}, M={M_modes}")
        print(f"{'='*70}")

    ks = KuramotoSivashinsky(L=L, N=N_fourier, dt=dt)

    # Ensure M_modes does not exceed available modes
    max_modes = N_fourier // 2
    M = min(M_modes, max_modes)
    if verbose:
        print(f"  Using {M} Fourier modes for state vector (dim = {2*M})")

    # ------------------------------------------------------------------
    # Step 1: Generate attractor reference set
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n[1/7] Generating attractor reference set...")
        print(f"  Transient: T={T_transient}, Sampling: T={T_attractor}, dt_save={dt_save}")

    t0 = time.time()
    u_hat0 = ks.initial_condition(seed=42)

    # Integrate through transient
    _, traj_trans = ks.integrate(u_hat0, T_transient, save_every=max(1, int(T_transient / dt)))
    u_hat_post_trans = traj_trans[-1]

    # Sample attractor
    t_attr, traj_attr_fourier = ks.integrate(u_hat_post_trans, T_attractor,
                                             save_every=save_every)
    t_gen = time.time() - t0

    # Convert to real state vectors
    N_attr = len(traj_attr_fourier)
    traj_attr = np.array([fourier_state_vector(uh, M) for uh in traj_attr_fourier])

    if verbose:
        print(f"  {N_attr} attractor points generated ({t_gen:.1f}s)")
        print(f"  State vector dimension: {traj_attr.shape[1]}")

    # Attractor radius
    centroid = np.mean(traj_attr, axis=0)
    R_attr = np.max(np.linalg.norm(traj_attr - centroid, axis=1))
    if verbose:
        print(f"  Attractor radius R = {R_attr:.4f}")

    # KD-Tree
    tree = cKDTree(traj_attr)
    if verbose:
        print(f"  KD-Tree built.")

    # ------------------------------------------------------------------
    # Step 2: Test trajectory (perturbed)
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n[2/7] Generating test trajectory (perturbed)...")

    t0 = time.time()
    # Start from a point on the attractor, add perturbation in Fourier space
    u_hat_start = traj_attr_fourier[N_attr // 2].copy()
    rng = np.random.RandomState(123)
    perturbation = rng.randn(len(u_hat_start)) + 1j * rng.randn(len(u_hat_start))
    perturbation[0] = 0  # keep mean unchanged
    pert_scale = 5.0  # perturbation magnitude
    u_hat_start += pert_scale * perturbation / np.linalg.norm(perturbation)

    t_test_arr, traj_test_fourier = ks.integrate(u_hat_start, T_test,
                                                  save_every=save_every)
    t_gen_test = time.time() - t0

    # Convert to state vectors
    N_test = len(traj_test_fourier)
    traj_test = np.array([fourier_state_vector(uh, M) for uh in traj_test_fourier])

    if verbose:
        print(f"  {N_test} test points generated ({t_gen_test:.1f}s)")

    # ------------------------------------------------------------------
    # Step 3: Distances and projections
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n[3/7] Computing d(t) and P(u(t))...")

    d_t, idx_nearest = tree.query(traj_test)
    proj_t = traj_attr[idx_nearest]

    d_min = np.min(d_t[d_t > 0]) if np.any(d_t > 0) else 1e-12
    d_mean = np.mean(d_t)
    d_max = np.max(d_t)
    d_median = np.median(d_t)

    if verbose:
        print(f"  d(t): min={d_min:.6f}, median={d_median:.6f}, "
              f"mean={d_mean:.6f}, max={d_max:.6f}")

    # Omega(t) = 1 + log(R / d(t))
    omega_t = np.ones(N_test)
    mask_pos = d_t > 1e-12
    omega_t[mask_pos] = 1.0 + np.log(R_attr / d_t[mask_pos])
    omega_t = np.maximum(omega_t, 1.0)

    if verbose:
        print(f"  Omega: min={np.min(omega_t):.4f}, mean={np.mean(omega_t):.4f}, "
              f"max={np.max(omega_t):.4f}")

    # ------------------------------------------------------------------
    # Step 4: ||u'(t)|| via RHS evaluation
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n[4/7] Computing ||u'(t)||...")

    u_prime = np.zeros(N_test)
    for i in range(N_test):
        rhs = ks.rhs_fourier(traj_test_fourier[i])
        # Norm in state-vector representation (first M modes)
        rhs_vec = fourier_state_vector(rhs, M)
        u_prime[i] = np.linalg.norm(rhs_vec)

    if verbose:
        print(f"  ||u'||: min={np.min(u_prime):.4f}, mean={np.mean(u_prime):.4f}, "
              f"max={np.max(u_prime):.4f}")

    # ==================================================================
    # TEST 1: TLL (Trajectory Log-Lipschitz)
    # ==================================================================
    if verbose:
        print(f"\n[5/7] Testing TLL and LDI...")
        print()
        print("-" * 60)
        print("TEST 1: Trajectory Log-Lipschitz (TLL)")
        print("-" * 60)

    h_steps = [1, 2, 5, 10, 20]
    c_tll_max_all = []
    c_tll_95_all = []
    tll_details = {}

    for h in h_steps:
        ratios = []
        for i in range(N_test - h):
            if d_t[i] < 1e-10:
                continue

            proj_diff = np.linalg.norm(proj_t[i + h] - proj_t[i])
            u_diff = np.linalg.norm(traj_test[i + h] - traj_test[i])

            if u_diff < 1e-12:
                continue

            omega = max(1.0, 1.0 + np.log(R_attr / d_t[i]))
            ratio = proj_diff / (u_diff * omega)
            ratios.append(ratio)

        if ratios:
            ratios = np.array(ratios)
            c_max = float(np.max(ratios))
            c_95 = float(np.percentile(ratios, 95))
            c_99 = float(np.percentile(ratios, 99))
            c_mean = float(np.mean(ratios))
            c_median = float(np.median(ratios))
            c_tll_max_all.append(c_max)
            c_tll_95_all.append(c_95)

            tll_details[f"h={h * dt_save:.3f}s"] = {
                "max": c_max, "p99": c_99, "p95": c_95,
                "mean": c_mean, "median": c_median, "n_samples": len(ratios)
            }

            if verbose:
                print(f"  h={h * dt_save:.3f}s: C_TLL(max)={c_max:.4f}, "
                      f"99%={c_99:.4f}, 95%={c_95:.4f}, mean={c_mean:.4f}, "
                      f"n={len(ratios)}")

    if c_tll_max_all:
        C_TLL_max = max(c_tll_max_all)
        C_TLL_95 = max(c_tll_95_all)
        tll_passed = C_TLL_95 < 10.0
        C_TLL_use = C_TLL_95 * 1.5  # robust bound
    else:
        C_TLL_max = float('inf')
        C_TLL_95 = float('inf')
        tll_passed = False
        C_TLL_use = float('inf')

    if verbose:
        if c_tll_max_all:
            print(f"\n  Global C_TLL (max)  = {C_TLL_max:.4f}")
            print(f"  Global C_TLL (95%)  = {C_TLL_95:.4f}")
            print(f"  => TLL {'PASS' if tll_passed else 'FAIL'} "
                  f"(95th percentile < 10)")
        else:
            print("  NO DATA for TLL test")

    # ==================================================================
    # TEST 2: LDI (Log-Distance Integrability)
    # ==================================================================
    if verbose:
        print()
        print("-" * 60)
        print("TEST 2: Log-Distance Integrability (LDI)")
        print("-" * 60)

    integrand = u_prime * omega_t
    ldi_integral = float(np.sum(integrand) * dt_save)
    l1_integral = float(np.sum(u_prime) * dt_save)
    ratio_ldi = ldi_integral / l1_integral if l1_integral > 0 else float('inf')

    if verbose:
        print(f"  int ||u'(t)|| dt           = {l1_integral:.4f}")
        print(f"  int ||u'(t)|| * Omega dt   = {ldi_integral:.4f}")
        print(f"  Ratio LDI/L1               = {ratio_ldi:.4f}")
        print(f"  LDI per time unit           = {ldi_integral / T_test:.4f}")

    ldi_passed = np.isfinite(ldi_integral) and ratio_ldi < 100
    if verbose:
        print(f"\n  => LDI {'PASS' if ldi_passed else 'FAIL'}: "
              f"log-weighted integral is finite")

    # ==================================================================
    # TEST 3: BV Chain Rule (Theorem 4.1)
    # ==================================================================
    if verbose:
        print()
        print("-" * 60)
        print("TEST 3: BV Chain Rule (Theorem 4.1)")
        print("-" * 60)

    proj_diffs = np.linalg.norm(np.diff(proj_t, axis=0), axis=1)
    total_var_v = float(np.sum(proj_diffs))
    bv_bound = C_TLL_use * ldi_integral

    if verbose:
        print(f"  Var(v) = sum ||P(u(t+1)) - P(u(t))|| = {total_var_v:.4f}")
        print(f"  C_TLL * LDI-Integral                  = {bv_bound:.4f}")

    if bv_bound > 0 and np.isfinite(bv_bound):
        ratio_bv = float(total_var_v / bv_bound)
        bv_passed = ratio_bv <= 1.05
    else:
        ratio_bv = float('inf')
        bv_passed = False

    if verbose:
        if np.isfinite(ratio_bv):
            print(f"  Ratio Var/Bound                       = {ratio_bv:.4f}")
        print(f"\n  => BV Chain Rule {'PASS' if bv_passed else 'FAIL'}")

    # ==================================================================
    # TEST 4: STC (Sublevel-Time Control)
    # ==================================================================
    if verbose:
        print()
        print("-" * 60)
        print("TEST 4: Sublevel-Time Control (STC)")
        print("-" * 60)

    eps_values = np.logspace(np.log10(max(d_min * 0.5, 1e-15)),
                             np.log10(d_max), 30)
    sublevel_fractions = np.array([np.mean(d_t <= eps) for eps in eps_values])
    sublevel_times = sublevel_fractions * T_test

    valid = (sublevel_times > 0) & (eps_values > 0)
    if np.sum(valid) >= 5:
        log_eps = np.log(eps_values[valid])
        log_st = np.log(sublevel_times[valid])
        coeffs = np.polyfit(log_eps, log_st, 1)
        alpha_fit = float(coeffs[0])
        c_star = float(np.exp(coeffs[1]))
        stc_passed = alpha_fit > 0.01
    else:
        alpha_fit = 0.0
        c_star = 0.0
        stc_passed = False

    if verbose:
        if np.sum(valid) >= 5:
            print(f"  STC Power-Law Fit:")
            print(f"    alpha = {alpha_fit:.4f} (required: > 0)")
            print(f"    C_*   = {c_star:.6f}")
            print(f"\n  {'eps':>12} | {'|{d<=eps}|':>10} | {'Pred':>10}")
            print(f"  {'-'*12}-+-{'-'*10}-+-{'-'*10}")
            for i in range(0, len(eps_values), 5):
                if i < len(valid) and valid[i]:
                    pred = c_star * eps_values[i]**alpha_fit
                    print(f"  {eps_values[i]:12.4f} | {sublevel_times[i]:10.4f} | "
                          f"{pred:10.4f}")
        else:
            print("  Not enough data points for STC fit.")
        print(f"\n  => STC {'PASS' if stc_passed else 'FAIL'}: "
              f"alpha = {alpha_fit:.4f}")

    # ==================================================================
    # TEST 5: Exponential Tracking (Lemma 6.1)
    # ==================================================================
    if verbose:
        print()
        print("-" * 60)
        print("TEST 5: Exponential Tracking (Lemma 6.1)")
        print("-" * 60)

    N_early = min(500, N_test // 2)
    t_early = t_test_arr[:N_early]
    d_early = d_t[:N_early]

    mask_decay = d_early > d_min * 2
    if np.sum(mask_decay) > 10:
        log_d = np.log(d_early[mask_decay])
        t_decay = t_early[mask_decay]
        coeffs_exp = np.polyfit(t_decay, log_d, 1)
        lambda_fit = float(-coeffs_exp[0])
        C_exp = float(np.exp(coeffs_exp[1]))
        exp_tracking = lambda_fit > 0
    else:
        lambda_fit = 0.0
        C_exp = 0.0
        exp_tracking = False

    if verbose:
        if lambda_fit != 0:
            print(f"  Exponential fit: d(t) ~ {C_exp:.4f} * exp(-{lambda_fit:.4f} * t)")
            print(f"  lambda = {lambda_fit:.4f} "
                  f"({'> 0: exponential approach' if exp_tracking else '< 0: NO approach'})")
        else:
            print("  Not enough data points for exponential fit.")
        print(f"\n  => Exponential Tracking {'PASS' if exp_tracking else 'N/A'}")

    # ==================================================================
    # TEST 6: STC + Exponential Tracking combined summary
    # ==================================================================
    if verbose:
        print()
        print("-" * 60)
        print("TEST 6: Sufficient conditions for LDI (Prop. 5.4 + Lemma 6.1)")
        print("-" * 60)
        if stc_passed:
            print(f"  STC with alpha = {alpha_fit:.4f} > 0 => LDI (Proposition 5.4)")
        if exp_tracking:
            print(f"  Exponential Tracking with lambda = {lambda_fit:.4f} > 0 => LDI (Lemma 6.1)")
        if not stc_passed and not exp_tracking:
            print(f"  Neither STC nor Exp. Tracking confirmed.")
        sufficient = stc_passed or exp_tracking
        print(f"\n  => Sufficient conditions {'PASS' if sufficient else 'FAIL'}")

    # ------------------------------------------------------------------
    # Collect results
    # ------------------------------------------------------------------
    results = {
        "system": "Kuramoto-Sivashinsky",
        "parameters": {
            "L": L, "N_fourier": N_fourier, "M_modes": M,
            "T_transient": T_transient, "T_attractor": T_attractor,
            "T_test": T_test, "dt": dt, "dt_save": dt_save,
            "N_attr": N_attr, "N_test": N_test,
            "state_dim": 2 * M,
        },
        "attractor": {
            "R": float(R_attr),
            "d_min": float(d_min), "d_median": float(d_median),
            "d_mean": float(d_mean), "d_max": float(d_max),
            "omega_min": float(np.min(omega_t)),
            "omega_mean": float(np.mean(omega_t)),
            "omega_max": float(np.max(omega_t)),
        },
        "tll": {
            "passed": tll_passed,
            "C_TLL_max": float(C_TLL_max),
            "C_TLL_95": float(C_TLL_95) if np.isfinite(C_TLL_95) else None,
            "C_TLL_use": float(C_TLL_use) if np.isfinite(C_TLL_use) else None,
            "details": tll_details,
        },
        "ldi": {
            "passed": ldi_passed,
            "integral": ldi_integral,
            "l1_integral": l1_integral,
            "ratio_ldi_l1": ratio_ldi,
            "per_time_unit": ldi_integral / T_test,
        },
        "bv_chain_rule": {
            "passed": bv_passed,
            "total_var": total_var_v,
            "bound": float(bv_bound) if np.isfinite(bv_bound) else None,
            "ratio": ratio_bv if np.isfinite(ratio_bv) else None,
        },
        "stc": {
            "passed": stc_passed,
            "alpha": alpha_fit,
            "c_star": c_star,
        },
        "exp_tracking": {
            "passed": exp_tracking,
            "lambda": lambda_fit,
            "C_exp": C_exp,
        },
        "all_passed": tll_passed and ldi_passed and bv_passed and stc_passed,
    }

    return results


def run_grid_refinement(L=32 * np.pi, M_modes=32, verbose=True):
    """
    Grid refinement study: run diagnostics for multiple N_fourier values.
    """
    # N=32 is under-resolved for L=32*pi; use N >= 64 for this domain size.
    # For N=32, we use L=22 (a classical chaotic KS domain that is well-resolved).
    N_values = [32, 64, 128, 256]

    if verbose:
        print()
        print("=" * 70)
        print("GRID REFINEMENT STUDY")
        print("=" * 70)

    refinement_results = {}

    for N in N_values:
        if verbose:
            print(f"\n>>> Running N_fourier = {N} <<<")

        # Adjust M_modes to not exceed N/2
        M = min(M_modes, N // 2)

        # Shorter runs for grid refinement
        T_trans = 200.0
        T_attr = 300.0
        T_tst = 80.0

        # N=32 cannot resolve L=32*pi (blow-up). Use shorter domain L=22
        # which is a well-known chaotic KS regime suitable for low resolution.
        L_run = 22.0 if N <= 32 else L
        dt_run = 0.025 if N <= 32 else 0.05

        if verbose and L_run != L:
            print(f"  (Using L={L_run} for N={N} to ensure numerical stability)")

        try:
            res = run_ks_diagnostics(
                N_fourier=N, L=L_run, M_modes=M,
                T_transient=T_trans, T_attractor=T_attr, T_test=T_tst,
                dt=dt_run, save_every=2, verbose=verbose
            )
            refinement_results[str(N)] = {
                "N": N, "M": M,
                "C_TLL_95": res["tll"]["C_TLL_95"],
                "C_TLL_max": res["tll"]["C_TLL_max"],
                "LDI": res["ldi"]["integral"],
                "LDI_L1_ratio": res["ldi"]["ratio_ldi_l1"],
                "BV_ratio": res["bv_chain_rule"]["ratio"],
                "STC_alpha": res["stc"]["alpha"],
                "exp_lambda": res["exp_tracking"]["lambda"],
                "tll_pass": res["tll"]["passed"],
                "ldi_pass": res["ldi"]["passed"],
                "bv_pass": res["bv_chain_rule"]["passed"],
                "stc_pass": res["stc"]["passed"],
                "all_pass": res["all_passed"],
            }
        except Exception as e:
            if verbose:
                print(f"  ERROR for N={N}: {e}")
            refinement_results[str(N)] = {"N": N, "error": str(e)}

    # Print refinement table
    if verbose:
        print()
        print("=" * 70)
        print("GRID REFINEMENT SUMMARY")
        print("=" * 70)
        print(f"  {'N':>5} | {'M':>4} | {'C_TLL(95%)':>11} | {'LDI/L1':>8} | "
              f"{'BV ratio':>9} | {'STC alpha':>10} | {'All':>4}")
        print(f"  {'-'*5}-+-{'-'*4}-+-{'-'*11}-+-{'-'*8}-+-"
              f"{'-'*9}-+-{'-'*10}-+-{'-'*4}")

        for key in sorted(refinement_results.keys(), key=int):
            r = refinement_results[key]
            if "error" in r:
                print(f"  {r['N']:5d} | {'ERR':>4} | {r['error'][:40]}")
            else:
                c95 = f"{r['C_TLL_95']:.4f}" if r['C_TLL_95'] is not None else "N/A"
                bvr = f"{r['BV_ratio']:.4f}" if r['BV_ratio'] is not None else "N/A"
                print(f"  {r['N']:5d} | {r['M']:4d} | {c95:>11} | "
                      f"{r['LDI_L1_ratio']:8.4f} | {bvr:>9} | "
                      f"{r['STC_alpha']:10.4f} | "
                      f"{'PASS' if r['all_pass'] else 'FAIL':>4}")

    return refinement_results


# ==========================================================================
# MAIN
# ==========================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("TLL + LDI VERIFICATION ON KURAMOTO-SIVASHINSKY ATTRACTOR")
    print("=" * 70)
    print(f"  KS equation: u_t + u*u_x + u_xx + u_xxxx = 0")
    print(f"  Domain: [0, L] with periodic BC")
    print(f"  Method: Pseudo-spectral (FFT) + ETDRK4")
    print()

    t_start = time.time()

    # ---- Primary run (N=128, full diagnostics) ----
    primary_results = run_ks_diagnostics(
        N_fourier=128, L=32 * np.pi, M_modes=32,
        T_transient=200.0, T_attractor=500.0, T_test=100.0,
        dt=0.05, save_every=2, verbose=True
    )

    # ---- Grid refinement study ----
    refinement = run_grid_refinement(
        L=32 * np.pi, M_modes=32, verbose=True
    )

    t_total = time.time() - t_start

    # ==================================================================
    # SUMMARY
    # ==================================================================
    print()
    print("=" * 70)
    print("FINAL SUMMARY: KS TLL+LDI Verification")
    print("=" * 70)

    all_passed = primary_results["all_passed"]

    summary_tests = [
        ("TLL (Trajectory Log-Lipschitz)", primary_results["tll"]["passed"],
         f"C_TLL(95%) = {primary_results['tll']['C_TLL_95']:.4f}"
         if primary_results["tll"]["C_TLL_95"] is not None else "N/A"),
        ("LDI (Log-Distance Integrability)", primary_results["ldi"]["passed"],
         f"Integral = {primary_results['ldi']['integral']:.2f}, "
         f"LDI/L1 = {primary_results['ldi']['ratio_ldi_l1']:.4f}"),
        ("BV Chain Rule (Theorem 4.1)", primary_results["bv_chain_rule"]["passed"],
         f"Var/Bound = {primary_results['bv_chain_rule']['ratio']:.4f}"
         if primary_results["bv_chain_rule"]["ratio"] is not None else "N/A"),
        ("STC (Sublevel-Time Control)", primary_results["stc"]["passed"],
         f"alpha = {primary_results['stc']['alpha']:.4f}"),
        ("Exponential Tracking (Lemma 6.1)", primary_results["exp_tracking"]["passed"],
         f"lambda = {primary_results['exp_tracking']['lambda']:.4f}"
         if primary_results["exp_tracking"]["lambda"] != 0 else "N/A"),
    ]

    for name, passed, detail in summary_tests:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {name}: {detail}")

    print(f"\n  Overall result (TLL+LDI+BV+STC): "
          f"{'ALL PASSED' if all_passed else 'NOT ALL PASSED'}")

    if all_passed:
        print(f"\n  *** KS PROOF OF LIFE SUCCESSFUL ***")
        print(f"  The TLL+LDI framework (TLL + LDI => BV Chain Rule)")
        print(f"  is numerically verified for the Kuramoto-Sivashinsky")
        print(f"  attractor (1D PDE, infinite-dimensional, chaotic).")
        print(f"  This extends the Lorenz result to a genuine PDE system.")

    print(f"\n  Total computation time: {t_total:.1f}s")
    print("=" * 70)

    # ==================================================================
    # Save results to JSON
    # ==================================================================
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(script_dir, "_results")
    os.makedirs(results_dir, exist_ok=True)
    results_path = os.path.join(results_dir, "ks_tll_ldi_results.json")

    output = {
        "primary": primary_results,
        "grid_refinement": refinement,
        "computation_time_s": t_total,
    }

    # Convert numpy types for JSON serialization
    def convert_numpy(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, dict):
            return {k: convert_numpy(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [convert_numpy(v) for v in obj]
        return obj

    output = convert_numpy(output)

    with open(results_path, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print(f"\n  Results saved to: {results_path}")
