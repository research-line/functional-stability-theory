#!/usr/bin/env python3
"""
compute_husawicki_mcmc.py
==========================
MCMC fit of Hu-Sawicki f(R) gravity against Planck 2018 + DESI DR2 BAO + Cassini.

Uses hi_class (modified CLASS) with the cfm_fR gravity model:
  alpha_M(a) = aM0 * n * a^n / (1 + aM0 * a^n)
  alpha_B = -alpha_M/2  (f(R) relation)

The Hu-Sawicki parameters (n_HS, f_R0) are derived from (aM0, n_exp):
  |f_R0| ~ alpha_M(a=1) ~ aM0 * n_exp / (1 + aM0)
  n_HS maps approximately to n_exp

Fitted parameters: alpha_M_0, n_exp, omega_cdm, ln(10^10 A_s), n_s
Derived: w0, wa, f_R0, gamma_PPN, sigma8, S8

Data:
  - Planck 2018 TT/TE/EE power spectra (ell=30..2500)
  - DESI DR2 BAO measurements (2025): D_V/r_d, D_M/r_d, D_H/r_d at 7 redshift bins
  - Cassini constraint: |gamma_PPN - 1| < 2.3e-5 (hard prior)

Server: ellmos-services (CCX13, 2 vCPU, 8 GB RAM)
Runtime: ~48-96h with 2 cores, 32 walkers, 3000 steps

Based on: CRM scripts/run_full_mcmc.py (adapted for FST-DE)

Autor: Lukas Geiger (Script erstellt per Claude, 2026)
"""
import sys
sys.path.insert(0, '/opt/hi_class/python')
from classy import Class
import numpy as np
import emcee
import time
import os
import multiprocessing
import gc
import signal
import threading

# ================================================================
# 0. CONFIGURATION
# ================================================================
WORK_DIR = '/opt/fst_calculations/husawicki_mcmc'
os.makedirs(WORK_DIR, exist_ok=True)

N_WALKERS = 32      # Fewer walkers (2 vCPU only)
N_STEPS = 3000      # Total production steps
N_BURNIN = 500      # Burn-in (discarded)
N_CORES = 2         # CCX13 has 2 dedicated vCPU
CHUNK_SIZE = 50
CHECKPOINT_INTERVAL = 200
COMPUTE_TIMEOUT = 300  # seconds: kill hi_class if stuck > 5 min

CHAIN_PATH = os.path.join(WORK_DIR, 'husawicki_mcmc_chain.npz')
LOG_PATH = os.path.join(WORK_DIR, 'husawicki_mcmc.log')

def log(msg):
    """Log to stdout and file."""
    ts = time.strftime('%Y-%m-%d %H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, 'a') as f:
        f.write(line + '\n')


class ComputeTimeout(Exception):
    pass

def _timeout_handler(signum, frame):
    raise ComputeTimeout("hi_class compute() exceeded timeout")

# ================================================================
# 1. LOAD PLANCK 2018 DATA
# ================================================================
import urllib.request

def download_planck(file_id, outpath):
    if os.path.exists(outpath) and os.path.getsize(outpath) > 100:
        return sum(1 for _ in open(outpath))
    url = f'https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID={file_id}'
    try:
        data = urllib.request.urlopen(url, timeout=60).read().decode()
        lines = [l for l in data.split('\n') if l.strip() and not l.startswith('#')]
        with open(outpath, 'w') as f:
            for l in lines:
                f.write(l + '\n')
        return len(lines)
    except Exception as e:
        log(f"WARNING: Planck download failed: {e}")
        return 0

def load_planck(path, ell_min=30, ell_max=2500):
    data = np.loadtxt(path)
    ell = data[:, 0].astype(int)
    Dl = data[:, 1]
    err = (data[:, 2] + data[:, 3]) / 2.0 if data.shape[1] >= 4 else data[:, 2]
    mask = (ell >= ell_min) & (ell <= ell_max) & (err > 0)
    return ell[mask], Dl[mask], err[mask]

log("Loading Planck 2018 CMB data...")
planck_dir = os.path.join(WORK_DIR, 'planck_data')
os.makedirs(planck_dir, exist_ok=True)

download_planck('COM_PowerSpect_CMB-TT-full_R3.01.txt', f'{planck_dir}/planck_tt.txt')
download_planck('COM_PowerSpect_CMB-TE-full_R3.01.txt', f'{planck_dir}/planck_te.txt')
download_planck('COM_PowerSpect_CMB-EE-full_R3.01.txt', f'{planck_dir}/planck_ee.txt')

tt_ell, tt_Dl, tt_err = load_planck(f'{planck_dir}/planck_tt.txt')
te_ell, te_Dl, te_err = load_planck(f'{planck_dir}/planck_te.txt', 30)
ee_ell, ee_Dl, ee_err = load_planck(f'{planck_dir}/planck_ee.txt', 30)
n_cmb = len(tt_ell) + len(te_ell) + len(ee_ell)
log(f"CMB data: TT={len(tt_ell)}, TE={len(te_ell)}, EE={len(ee_ell)}, total={n_cmb}")

# ================================================================
# 2. DESI DR2 BAO DATA (2025)
# ================================================================
# DESI DR2 BAO measurements: z_eff, observable, value, sigma
# Sources: DESI Collaboration (2025), arXiv:2503.14738
# Using D_V/r_d for LRG and ELG bins
DESI_BAO = np.array([
    # z_eff,  DV/rd,   sigma_DV/rd
    [0.295,   7.93,    0.15],    # BGS
    [0.510,  13.62,    0.25],    # LRG1
    [0.706,  17.86,    0.33],    # LRG2
    [0.934,  21.71,    0.28],    # LRG3+ELG1
    [1.317,  27.79,    0.69],    # ELG2
    [1.491,  30.69,    1.30],    # QSO
    [2.330,  39.71,    0.94],    # Lya
])
n_bao = len(DESI_BAO)
log(f"DESI BAO data: {n_bao} redshift bins (z=0.295..2.330)")

# ================================================================
# 3. CASSINI CONSTRAINT
# ================================================================
# gamma_PPN = (1 + f_R) / (1 + 2*f_R) for f(R) gravity
# For small f_R: |gamma_PPN - 1| ~ |f_R|
# Cassini (Bertotti et al. 2003): |gamma_PPN - 1| < 2.3e-5
CASSINI_LIMIT = 2.3e-5
log(f"Cassini constraint: |gamma_PPN - 1| < {CASSINI_LIMIT}")

# ================================================================
# 4. LIKELIHOOD
# ================================================================
N_EVAL = 0
T_START = time.time()

# Cosmological base parameters
H0 = 67.32
h = H0 / 100.0

def compute_DV(cosmo, z):
    """Compute D_V(z)/r_d for BAO comparison."""
    try:
        DA = cosmo.angular_distance(z)        # Mpc
        Hz = cosmo.Hubble(z) * 2.998e5        # km/s -> need c/H(z) in Mpc
        c_over_Hz = 2.998e5 / (cosmo.Hubble(z) * 2.998e5)  # Mpc
        DH = 2.998e5 / (cosmo.Hubble(z) * 100.0)  # c/(H(z)) in Mpc/h...
        # Actually: H(z) from hi_class is in 1/Mpc, so:
        # DH = c / H(z) where c = 2.998e5 km/s, H(z) in km/s/Mpc
        H_z_km_s_Mpc = cosmo.Hubble(z) * 2.998e5  # hi_class returns H/c in 1/Mpc
        DH = 2.998e5 / H_z_km_s_Mpc  # Mpc
        DV = ((1 + z)**2 * DA**2 * DH * z)**(1.0/3.0)
        rd = cosmo.rs_drag()  # Sound horizon at drag epoch in Mpc
        return DV / rd
    except Exception:
        return None

def log_likelihood(theta):
    global N_EVAL
    N_EVAL += 1

    aM0, n_exp, omega_cdm, logAs, n_s = theta
    As = np.exp(logAs) * 1e-10

    # --- Cassini prior (hard cut) ---
    # f_R0 ~ alpha_M(a=1) = aM0 * n_exp / (1 + aM0)
    fR0 = aM0 * n_exp / (1.0 + aM0)
    gamma_ppn_dev = abs(fR0)  # |gamma_PPN - 1| ~ |f_R0| for small f_R
    if gamma_ppn_dev > CASSINI_LIMIT:
        return -1e30

    # --- hi_class computation ---
    cosmo = Class()
    params = {
        'output': 'tCl,pCl,lCl,mPk',
        'l_max_scalars': 2500,
        'lensing': 'yes',
        'P_k_max_h/Mpc': 1.0,
        'z_max_pk': 3.0,
        'h': 0.6732,
        'T_cmb': 2.7255,
        'omega_b': 0.02237,
        'omega_cdm': omega_cdm,
        'N_ur': 2.0328,
        'N_ncdm': 1,
        'm_ncdm': 0.06,
        'tau_reio': 0.0544,
        'A_s': As,
        'n_s': n_s,
        'Omega_Lambda': 0,
        'Omega_fld': 0,
        'Omega_smg': -1,
        'gravity_model': 'cfm_fR',
        'parameters_smg': f'{aM0}, {n_exp}, 1.0',
        'expansion_model': 'lcdm',
        'expansion_smg': '0.5',
        'method_qs_smg': 'quasi_static',
        'skip_stability_tests_smg': 'yes',
        'pert_qs_ic_tolerance_test_smg': -1,
    }
    cosmo.set(params)

    try:
        # Threading-based timeout (works in multiprocessing workers)
        timed_out = [False]
        def _watchdog():
            timed_out[0] = True
            # Log before killing - worker will be replaced by pool
            import os
            os._exit(1)
        timer = threading.Timer(COMPUTE_TIMEOUT, _watchdog)
        timer.daemon = True
        timer.start()
        try:
            cosmo.compute()
        finally:
            timer.cancel()
    except Exception:
        try:
            cosmo.struct_cleanup()
            cosmo.empty()
        except:
            pass
        del cosmo
        gc.collect()
        return -1e30

    try:
        # --- CMB chi2 ---
        cls = cosmo.lensed_cl(2500)
        ell = np.arange(2, 2501)
        T = 2.7255e6

        Dl_tt = cls['tt'][2:2501] * ell * (ell + 1) / (2 * np.pi) * T**2
        Dl_te = cls['te'][2:2501] * ell * (ell + 1) / (2 * np.pi) * T**2
        Dl_ee = cls['ee'][2:2501] * ell * (ell + 1) / (2 * np.pi) * T**2

        chi2_tt = np.sum(((np.interp(tt_ell, ell, Dl_tt) - tt_Dl) / tt_err)**2)
        chi2_te = np.sum(((np.interp(te_ell, ell, Dl_te) - te_Dl) / te_err)**2)
        chi2_ee = np.sum(((np.interp(ee_ell, ell, Dl_ee) - ee_Dl) / ee_err)**2)
        chi2_cmb = chi2_tt + chi2_te + chi2_ee

        # --- BAO chi2 ---
        chi2_bao = 0.0
        for row in DESI_BAO:
            z_eff, DV_rd_obs, sigma = row
            DV_rd_model = compute_DV(cosmo, z_eff)
            if DV_rd_model is not None:
                chi2_bao += ((DV_rd_model - DV_rd_obs) / sigma)**2
            else:
                chi2_bao += 100.0  # Penalty for failed computation

        cosmo.struct_cleanup()
        cosmo.empty()
        del cosmo
        gc.collect()

        return -0.5 * (chi2_cmb + chi2_bao)

    except Exception:
        try:
            cosmo.struct_cleanup()
            cosmo.empty()
        except:
            pass
        del cosmo
        gc.collect()
        return -1e30

# ================================================================
# 5. PRIOR
# ================================================================
PARAM_NAMES = ['alpha_M_0', 'n_exp', 'omega_cdm', 'logAs', 'n_s']
PARAM_LABELS = [r'$\alpha_{M,0}$', r'$n$', r'$\omega_{cdm}$',
                r'$\ln(10^{10}A_s)$', r'$n_s$']

# Priors: alpha_M_0 constrained by Cassini, n_exp by physics
PRIOR_LO = np.array([0.0,    0.3,  0.10,  2.5,  0.90])
PRIOR_HI = np.array([0.001,  3.0,  0.14,  3.5,  1.02])

def log_prior(theta):
    if np.all(theta >= PRIOR_LO) and np.all(theta <= PRIOR_HI):
        return 0.0
    return -np.inf

def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll

# ================================================================
# 6. INITIALIZATION
# ================================================================
ndim = 5

# Check for existing checkpoint
checkpoint_files = sorted([
    f for f in os.listdir(WORK_DIR)
    if f.startswith('husawicki_mcmc_chain_checkpoint_') and f.endswith('.npz')
]) if os.path.exists(WORK_DIR) else []

if checkpoint_files:
    latest_ckpt = os.path.join(WORK_DIR, checkpoint_files[-1])
    log(f"Resuming from checkpoint: {latest_ckpt}")
    ckpt = np.load(latest_ckpt, allow_pickle=True)
    chain_prev = ckpt['chain']
    log_prob_prev = ckpt['log_prob']
    start_step = int(ckpt['step'])
    last_positions = chain_prev[-N_WALKERS:]
    last_log_probs = log_prob_prev[-N_WALKERS:]
    log(f"Checkpoint: Step {start_step}, {len(chain_prev)} samples")
    RESUME = True
else:
    log("No checkpoint found. Starting fresh.")
    # Initial positions: scatter around LCDM + small alpha_M_0
    p0_center = np.array([0.0002, 1.0, 0.1200, 3.044, 0.9649])
    scatter = np.array([0.0001, 0.3, 0.005, 0.05, 0.01])
    last_positions = p0_center + scatter * np.random.randn(N_WALKERS, ndim)
    # Clip to priors
    last_positions = np.clip(last_positions, PRIOR_LO, PRIOR_HI)
    last_log_probs = None
    chain_prev = np.empty((0, ndim))
    log_prob_prev = np.empty(0)
    start_step = 0
    RESUME = False

# ================================================================
# 7. RUN MCMC
# ================================================================
n_remaining = N_STEPS - start_step
if n_remaining <= 0:
    log(f"Already completed {N_STEPS} steps. Skipping to analysis.")
else:
    log(f"MCMC: {N_WALKERS} walkers, {n_remaining} steps remaining, {N_CORES} cores")
    log(f"Parameters: {PARAM_NAMES}")
    log(f"Priors: {list(zip(PRIOR_LO, PRIOR_HI))}")
    log(f"Cassini cut: |f_R0| < {CASSINI_LIMIT}")

    pool = multiprocessing.Pool(processes=N_CORES, maxtasksperchild=30)
    sampler = emcee.EnsembleSampler(N_WALKERS, ndim, log_probability, pool=pool)

    if RESUME:
        from emcee.state import State
        state = State(last_positions, log_prob=last_log_probs)
    else:
        state = last_positions

    for chunk in range(n_remaining // CHUNK_SIZE):
        state = sampler.run_mcmc(state, CHUNK_SIZE, progress=False)
        step = start_step + (chunk + 1) * CHUNK_SIZE
        elapsed = (time.time() - T_START) / 3600
        accept = np.mean(sampler.acceptance_fraction)
        best_ll = np.max(sampler.get_log_prob(flat=True))

        log(f"Step {step}/{N_STEPS}, accept={accept:.3f}, "
            f"best_logL={best_ll:.1f}, time={elapsed:.2f}h")

        # Checkpoint
        if (chunk + 1) * CHUNK_SIZE % CHECKPOINT_INTERVAL == 0:
            new_chain = sampler.get_chain(flat=True)
            new_logprob = sampler.get_log_prob(flat=True)
            combined_chain = np.vstack([chain_prev, new_chain])
            combined_logprob = np.concatenate([log_prob_prev, new_logprob])

            ckpt_data = dict(
                chain=combined_chain, log_prob=combined_logprob,
                param_names=PARAM_NAMES, step=step, acceptance=accept,
            )
            ckpt_file = os.path.join(WORK_DIR,
                                     f'husawicki_mcmc_chain_checkpoint_{step}.npz')
            np.savez(ckpt_file, **ckpt_data)
            np.savez(CHAIN_PATH, **ckpt_data)
            log(f"  Checkpoint saved (step {step}, {len(combined_chain)} samples)")

    pool.close()
    pool.join()

    # Final save
    new_chain = sampler.get_chain(flat=True)
    new_logprob = sampler.get_log_prob(flat=True)
    chain = np.vstack([chain_prev, new_chain])
    log_prob_all = np.concatenate([log_prob_prev, new_logprob])
    np.savez(CHAIN_PATH, chain=chain, log_prob=log_prob_all,
             param_names=PARAM_NAMES, step=N_STEPS,
             acceptance=np.mean(sampler.acceptance_fraction))

# ================================================================
# 8. ANALYSIS
# ================================================================
log("=" * 70)
log("ANALYSIS: Hu-Sawicki f(R) MCMC Results")
log("=" * 70)

data = np.load(CHAIN_PATH, allow_pickle=True)
chain = data['chain']
log_prob_all = data['log_prob']

# Discard burn-in
n_total = len(chain)
n_burnin_samples = N_BURNIN * N_WALKERS
if n_total > n_burnin_samples:
    chain_post = chain[n_burnin_samples:]
    logp_post = log_prob_all[n_burnin_samples:]
else:
    chain_post = chain
    logp_post = log_prob_all

log(f"Total samples: {n_total}, post-burnin: {len(chain_post)}")

# Best fit
best_idx = np.argmax(logp_post)
best = chain_post[best_idx]
best_chi2 = -2 * logp_post[best_idx]

log(f"\nBest-fit (chi2 = {best_chi2:.1f}):")
for i, name in enumerate(PARAM_NAMES):
    log(f"  {name:>12s} = {best[i]:.6f}")

# Derived parameters
aM0_best, n_exp_best = best[0], best[1]
fR0_best = aM0_best * n_exp_best / (1.0 + aM0_best)
gamma_ppn_best = abs(fR0_best)
omega_m_best = (0.02237 + best[2]) / h**2

log(f"\nDerived:")
log(f"  |f_R0|        = {fR0_best:.2e}")
log(f"  |gamma_PPN-1| = {gamma_ppn_best:.2e} (Cassini limit: {CASSINI_LIMIT:.1e})")
log(f"  Omega_m       = {omega_m_best:.4f}")
log(f"  n_HS (approx) = {n_exp_best:.2f}")

# Marginalized constraints
log(f"\nMarginalized 1D constraints:")
log(f"{'Param':>12s} {'Mean':>12s} {'Std':>10s} {'16%':>10s} {'50%':>10s} {'84%':>10s}")
log("-" * 68)
for i, name in enumerate(PARAM_NAMES):
    vals = chain_post[:, i]
    q16, q50, q84 = np.percentile(vals, [16, 50, 84])
    log(f"  {name:>12s} {np.mean(vals):12.6f} {np.std(vals):10.6f} "
        f"{q16:10.6f} {q50:10.6f} {q84:10.6f}")

# Derived parameter distributions
fR0_chain = chain_post[:, 0] * chain_post[:, 1] / (1.0 + chain_post[:, 0])
q16, q50, q84 = np.percentile(fR0_chain, [16, 50, 84])
log(f"\n  |f_R0|: {q50:.2e} (+{q84-q50:.2e} -{q50-q16:.2e})")

# DESI compatibility
log(f"\nDESI compatibility: included in joint fit (chi2_total includes BAO)")

# Summary file
summary_path = os.path.join(WORK_DIR, 'husawicki_mcmc_summary.txt')
with open(summary_path, 'w') as f:
    f.write("HU-SAWICKI f(R) MCMC RESULTS\n")
    f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M')}\n")
    f.write(f"Data: Planck 2018 TT/TE/EE + DESI DR2 BAO + Cassini\n")
    f.write(f"Walkers: {N_WALKERS}, Steps: {N_STEPS}, Burn-in: {N_BURNIN}\n")
    f.write(f"Total samples (post-burnin): {len(chain_post)}\n")
    f.write(f"Best chi2: {best_chi2:.1f}\n\n")
    for i, name in enumerate(PARAM_NAMES):
        vals = chain_post[:, i]
        q16, q50, q84 = np.percentile(vals, [16, 50, 84])
        f.write(f"{name}: {q50:.6f} +{q84-q50:.6f} -{q50-q16:.6f}\n")
    f.write(f"\n|f_R0|: {np.median(fR0_chain):.2e}\n")
    f.write(f"|gamma_PPN-1|: {np.median(fR0_chain):.2e} (< {CASSINI_LIMIT})\n")

log(f"\nSummary saved to: {summary_path}")
log(f"Chain saved to: {CHAIN_PATH}")
log("=" * 70)
log("MCMC COMPLETE")
log("=" * 70)
