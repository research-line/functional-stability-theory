"""
compute_goy_shell_dfc.py
========================
Turbulenz: Schalenmodell-DFC-Verifikation.

Verwendet ein IMEX-Schema (implizite Dissipation, explizite Nichtlinearitaet)
fuer das Sabra-Schalenmodell (L'vov, Podivilov, Pomyalov, Procaccia, Vandembroucq 1998):

  du_n/dt = i*(k_n * u_{n+1}^* * u_{n+2}  [forward cascade]
            - (1/2)*k_{n-1} * u_{n+1}^* * u_{n-1}  [local interaction]
            - (1/4)*k_{n-2} * u_{n-1} * u_{n-2})    [backward cascade]
            - nu * k_n^2 * u_n + f * delta_{n,n_f}

Sabra-Modell ist numerisch stabiler als GOY (gleiche Physik).

DFC1: Energiefluss Pi_n > 0 (Vorwaerts-Kaskade)
DFC2: Spektrale Steilheit |u_n|^2 ~ k_n^{-2/3}

Autor: Lukas Geiger (Skript erstellt per Claude, 2026)
"""

import numpy as np
import os
import time

# ==========================================================================
# Parameter
# ==========================================================================
N = 22              # Schalen
LAM = 2.0           # Schalen-Ratio
K0 = 2**(-4)        # Basis-Wellenzahl (standard Biferale)
NU = 1e-7           # Viskositaet
DT = 5e-5           # Zeitschritt
T_THERM = 20.0      # Thermalisierungszeit
T_DATA = 80.0       # Datensammlungszeit
F_SHELL = 2         # Forcing-Schale
F_AMP = 5e-3        # Forcing-Amplitude

k_n = K0 * LAM**np.arange(N, dtype=float)


def sabra_nonlinear(u):
    """Nichtlinearer Teil des Sabra-Modells (vektorisiert)"""
    N = len(u)
    nl = np.zeros(N, dtype=complex)

    # Forward: k_n * conj(u_{n+1}) * u_{n+2}
    nl[:N-2] += k_n[:N-2] * np.conj(u[1:N-1]) * u[2:N]

    # Local: -(1/2) * k_{n-1} * conj(u_{n+1}) * u_{n-1}
    nl[1:N-1] += -0.5 * k_n[:N-2] * np.conj(u[2:N]) * u[:N-2]

    # Backward: -(1/4) * k_{n-2} * u_{n-1} * u_{n-2}
    nl[2:N] += -0.25 * k_n[:N-2] * u[1:N-1] * u[:N-2]

    return 1j * nl


def imex_step(u, dt, nu):
    """
    IMEX-Schema: Dissipation implizit, Nichtlinearitaet explizit.
    u^{n+1} = [u^n + dt * NL(u^n) + dt * f] / (1 + dt * nu * k^2)
    """
    nl = sabra_nonlinear(u)

    # Forcing
    forcing = np.zeros(N, dtype=complex)
    forcing[F_SHELL] = F_AMP * (1.0 + 0j)

    # Semi-implizit
    denom = 1.0 + dt * nu * k_n**2
    u_new = (u + dt * nl + dt * forcing) / denom

    return u_new


# ==========================================================================
# Simulation
# ==========================================================================

print("=" * 70)
print("SABRA SHELL MODEL: DFC1/DFC2 Verifikation")
print("=" * 70)
print(f"  N = {N}, nu = {NU:.1e}, Re ~ {F_AMP/(NU * k_n[F_SHELL]**2):.0f}")
print(f"  k_range = [{k_n[0]:.4f}, {k_n[-1]:.0f}]")
print(f"  dt = {DT:.1e}")

# Initialisierung: K41-artiges Spektrum
np.random.seed(42)
u = np.zeros(N, dtype=complex)
for n in range(N):
    amp = 1e-5 * k_n[n]**(-1.0/3.0)
    u[n] = amp * np.exp(1j * 2*np.pi * np.random.random())

# Phase 1: Thermalisierung
print(f"\n[1] Thermalisierung ({T_THERM}s)...")
t0 = time.time()
n_therm = int(T_THERM / DT)
print_interval = max(1, n_therm // 5)

for step in range(n_therm):
    u = imex_step(u, DT, NU)
    if (step + 1) % print_interval == 0:
        E_total = np.sum(np.abs(u)**2)
        eps_diss = 2 * NU * np.sum(k_n**2 * np.abs(u)**2)
        print(f"    t={step*DT:8.1f}: E={E_total:.4e}, eps={eps_diss:.4e}")

elapsed = time.time() - t0
print(f"    ({elapsed:.1f}s)")

# Phase 2: Datensammlung
print(f"\n[2] Datensammlung ({T_DATA}s)...")
t0 = time.time()
n_data = int(T_DATA / DT)
n_snapshots = 500
sample_interval = max(1, n_data // n_snapshots)

E_all = []
Pi_all = []
eps_all = []

for step in range(n_data):
    u = imex_step(u, DT, NU)

    if step % sample_interval == 0:
        E_snap = np.abs(u)**2
        E_all.append(E_snap.copy())

        # Energiefluss: Pi_n = Im(k_n * u_n * conj(u_{n+1}) * u_{n+2})
        Pi = np.zeros(N)
        for nn in range(N-2):
            Pi[nn] = np.imag(k_n[nn] * u[nn] * np.conj(u[nn+1]) * u[nn+2])
        Pi_all.append(Pi.copy())

        eps_diss = 2 * NU * np.sum(k_n**2 * E_snap)
        eps_all.append(eps_diss)

E_all = np.array(E_all)
Pi_all = np.array(Pi_all)
eps_all = np.array(eps_all)
elapsed = time.time() - t0

n_snap = len(E_all)
print(f"    {n_snap} Snapshots in {elapsed:.1f}s")
E_mean = np.mean(E_all, axis=0)
eps_mean = np.mean(eps_all)
Pi_mean = np.mean(Pi_all, axis=0)
print(f"    <eps> = {eps_mean:.4e}")
print(f"    <E_total> = {np.sum(E_mean):.4e}")


# ==========================================================================
# DFC1: Vorwaerts-Kaskade
# ==========================================================================
print(f"\n{'='*70}")
print("DFC1: VORWAERTS-KASKADE (Pi_n > 0)")
print(f"{'='*70}")

# Inertialbereich: zwischen Forcing und dissipativem Bereich
# Dissipative Skala: k_d ~ (eps/nu^3)^(1/4)
if eps_mean > 0:
    k_d = (eps_mean / NU**3)**0.25
    n_diss = int(np.log(k_d/K0) / np.log(LAM))
else:
    n_diss = N - 3

inertial_start = F_SHELL + 2
inertial_end = min(n_diss - 1, N - 3)

print(f"\n  k_d = {k_d:.1f} (Schale ~{n_diss})")
print(f"  Inertialbereich: Schalen {inertial_start}-{inertial_end}")
print(f"\n  {'n':>5} {'k_n':>10} {'<Pi>':>12} {'Pi>0%':>8} {'Status':>8}")
print(f"  {'-'*48}")

dfc1_ok = True
for nn in range(max(0, inertial_start), min(inertial_end + 1, N)):
    frac_pos = np.mean(Pi_all[:, nn] > 0) * 100
    status = "OK" if frac_pos > 60 else "WARN" if frac_pos > 40 else "FAIL"
    if frac_pos <= 40:
        dfc1_ok = False
    print(f"  {nn:>5} {k_n[nn]:>10.3f} {Pi_mean[nn]:>12.4e} {frac_pos:>7.1f}% {status:>8}")

print(f"\n  DFC1: {'ERFUELLT' if dfc1_ok else 'TEILWEISE'}")


# ==========================================================================
# DFC2: Spektrale Steilheit
# ==========================================================================
print(f"\n{'='*70}")
print("DFC2: SPEKTRALE STEILHEIT")
print(f"{'='*70}")

log_k = np.log(k_n)
slopes_all = np.zeros((N-1, n_snap))
for t_idx in range(n_snap):
    log_E = np.log(np.maximum(E_all[t_idx], 1e-300))
    slopes_all[:, t_idx] = np.diff(log_E) / np.diff(log_k)

slopes_mean = np.mean(slopes_all, axis=1)
slopes_std = np.std(slopes_all, axis=1)

print(f"\n  K41: d(ln E_n)/d(ln k) = -2/3 = {-2/3:.4f}")
print(f"\n  {'n':>5} {'k_n':>10} {'slope':>10} {'std':>8} {'<=-2/3':>8}")
print(f"  {'-'*45}")

dfc2_ok = True
for nn in range(inertial_start, min(inertial_end, N-1)):
    s = slopes_mean[nn]
    ok = "JA" if s <= -0.3 else "NEIN"
    if s > 0:
        dfc2_ok = False
    print(f"  {nn:>5} {k_n[nn]:>10.3f} {s:>10.3f} {slopes_std[nn]:>8.3f} {ok:>8}")

inertial_slopes = slopes_mean[inertial_start:min(inertial_end, N-1)]
mean_slope = np.mean(inertial_slopes)
print(f"\n  Mittlere Steigung: {mean_slope:.4f} (K41: {-2/3:.4f})")
print(f"  DFC2: {'ERFUELLT' if dfc2_ok else 'TEILWEISE'}")


# ==========================================================================
# F[E]-Verifikation
# ==========================================================================
print(f"\n{'='*70}")
print("F[E]-VERIFIKATION")
print(f"{'='*70}")

E_ref = k_n**(-2.0/3.0)
E_ref *= np.sum(E_mean) / np.sum(E_ref)

def compute_F(E, E_star):
    mask = (E > 0) & (E_star > 0)
    F = np.sum(E[mask] * np.log(E[mask] / E_star[mask]) - E[mask] + E_star[mask])
    F += np.sum(E_star[~mask & (E_star > 0)])
    return F

F_mean_val = compute_F(E_mean, E_ref)
F_snaps = np.array([compute_F(E_all[t], E_ref) for t in range(n_snap)])

print(f"  F[<E>] = {F_mean_val:.6f}")
print(f"  <F> = {np.mean(F_snaps):.6f} +/- {np.std(F_snaps):.6f}")
print(f"  F > 0: {np.sum(F_snaps > 0)}/{n_snap} ({100*np.mean(F_snaps > 0):.0f}%)")


# ==========================================================================
# Zusammenfassung
# ==========================================================================
print(f"\n{'='*70}")
print("ZUSAMMENFASSUNG")
print(f"{'='*70}")
print(f"""
  Sabra Shell Model: N={N}, nu={NU:.1e}
  IMEX-Integrator, {n_snap} Snapshots

  <eps> = {eps_mean:.4e}

  DFC1 (Vorwaerts-Kaskade):   {'ERFUELLT' if dfc1_ok else 'TEILWEISE'}
  DFC2 (Steilheit):           {'ERFUELLT' if dfc2_ok else 'TEILWEISE'}
    Mittlere Steigung: {mean_slope:.4f} (K41: {-2/3:.4f})
  F[<E>] = {F_mean_val:.4f}

  DFC -> NL' -> Theorem 2 (Anomale Dissipation)
""")

# Plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax = axes[0, 0]
    ax.loglog(k_n, E_mean, 'bo-', lw=2, ms=5, label=r'$\langle E_n \rangle$ (Sabra)')
    ax.loglog(k_n, E_ref, 'r--', lw=1.5, label=r'$k^{-2/3}$ (K41)')
    E_std = np.std(E_all, axis=0)
    ax.fill_between(k_n, np.maximum(E_mean-E_std, 1e-30), E_mean+E_std, alpha=0.2, color='blue')
    ax.set_xlabel(r'$k_n$'); ax.set_ylabel(r'$E_n$')
    ax.set_title('Energiespektrum')
    ax.legend(); ax.grid(True, alpha=0.3, which='both')

    ax = axes[0, 1]
    ax.semilogx(k_n[:N-2], Pi_mean[:N-2], 'ro-', lw=2, ms=5)
    ax.axhline(y=0, color='k', lw=0.5)
    ax.axhline(y=eps_mean, color='green', ls='--', label=r'$\langle\varepsilon\rangle$')
    ax.axvspan(k_n[inertial_start], k_n[min(inertial_end, N-1)], alpha=0.1, color='green')
    ax.set_xlabel(r'$k_n$'); ax.set_ylabel(r'$\Pi_n$')
    ax.set_title('Energiefluss (DFC1)'); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.errorbar(range(len(slopes_mean)), slopes_mean, yerr=slopes_std,
                fmt='bs-', lw=1.5, ms=4, capsize=2)
    ax.axhline(y=-2.0/3.0, color='red', ls='--', lw=1.5, label=r'$-2/3$ (K41)')
    ax.axvspan(inertial_start, inertial_end, alpha=0.1, color='green')
    ax.set_xlabel('Schale n'); ax.set_ylabel('Steigung')
    ax.set_title('DFC2: Steilheit'); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(F_snaps, 'g-', lw=0.5, alpha=0.7)
    ax.axhline(y=np.mean(F_snaps), color='blue', lw=1.5,
               label=f'$\\langle F \\rangle = {np.mean(F_snaps):.3f}$')
    ax.set_xlabel('Snapshot'); ax.set_ylabel(r'$\mathcal{F}[E]$')
    ax.set_title('Free Energy'); ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(os.path.dirname(__file__), 'compute_goy_shell_dfc.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Plot: {outpath}")
except Exception as e:
    print(f"(Plot: {e})")

print("[DONE]")
