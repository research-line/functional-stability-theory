"""
DS3 Stress Test: Dissipative Selection on Lorenz Attractor
Tests whether TV(v_eps) remains bounded as eps -> 0
"""
import numpy as np
from scipy.integrate import solve_ivp
import os

# Lorenz parameters
SIGMA, RHO, BETA = 10.0, 28.0, 8.0/3.0
T_TOTAL = 50.0  # total integration time
T_TRANSIENT = 10.0  # discard transient
DT = 0.01

def lorenz(t, state):
    x, y, z = state
    return [SIGMA*(y-x), x*(RHO-z)-y, x*y-BETA*z]

def total_variation(traj):
    """Compute TV of a trajectory (sum of jumps)"""
    diffs = np.diff(traj, axis=0)
    return np.sum(np.sqrt(np.sum(diffs**2, axis=1)))

def compute_attractor_trajectory(n_points=5000):
    """Simulate Lorenz and return attractor trajectory"""
    sol = solve_ivp(lorenz, [0, T_TOTAL], [1.0, 1.0, 1.0],
                    t_eval=np.linspace(0, T_TOTAL, int(T_TOTAL/DT)),
                    method='RK45', rtol=1e-10, atol=1e-12)
    # Discard transient
    idx_start = int(T_TRANSIENT / DT)
    return sol.y[:, idx_start:].T  # shape: (N, 3)

def f_epsilon(u_traj, a_traj, eps, dt):
    """Compute F_eps[a] = (1/eps) * integral ||u-a||^2 dt + TV(a)"""
    fidelity = (1.0/eps) * np.sum(np.sum((u_traj - a_traj)**2, axis=1)) * dt
    tv = total_variation(a_traj)
    return fidelity, tv, fidelity + tv

def nearest_on_attractor(point, attractor_pts):
    """Find nearest point on attractor (discrete approximation)"""
    dists = np.sum((attractor_pts - point)**2, axis=1)
    return attractor_pts[np.argmin(dists)]

def minimize_f_epsilon(u_traj, attractor_pts, eps, dt, method='greedy'):
    """
    Approximate minimizer of F_eps over trajectories on attractor.
    Greedy method: at each time step, balance fidelity vs TV.
    """
    N = len(u_traj)
    v = np.zeros_like(u_traj)

    # Start with nearest point
    v[0] = nearest_on_attractor(u_traj[0], attractor_pts)

    for i in range(1, N):
        # Candidate 1: stay (no TV cost)
        stay = v[i-1]
        fid_stay = np.sum((u_traj[i] - stay)**2) / eps

        # Candidate 2: jump to nearest on attractor
        nearest = nearest_on_attractor(u_traj[i], attractor_pts)
        fid_nearest = np.sum((u_traj[i] - nearest)**2) / eps
        tv_jump = np.sqrt(np.sum((nearest - v[i-1])**2))

        # Choose: minimize local F_eps contribution
        if fid_stay * dt <= (fid_nearest * dt + tv_jump):
            v[i] = stay
        else:
            v[i] = nearest

    return v

# Main computation
print("=" * 60)
print("DS3 STRESS TEST: Lorenz Attractor")
print("=" * 60)

# Step 1: Generate attractor
print("\n[1] Generating Lorenz attractor...")
attractor = compute_attractor_trajectory(n_points=10000)
N = len(attractor)
dt = (T_TOTAL - T_TRANSIENT) / N
print(f"    Attractor points: {N}")
print(f"    TV(attractor trajectory): {total_variation(attractor):.2f}")

# Step 2: Generate "solution" u(t) = perturbed attractor trajectory
print("\n[2] Generating perturbed trajectory u(t)...")
# Use a slightly different initial condition
sol2 = solve_ivp(lorenz, [0, T_TOTAL], [1.1, 0.9, 1.0],
                 t_eval=np.linspace(0, T_TOTAL, int(T_TOTAL/DT)),
                 method='RK45', rtol=1e-10, atol=1e-12)
idx_start = int(T_TRANSIENT / DT)
u_traj = sol2.y[:, idx_start:].T[:N]
print(f"    u trajectory points: {len(u_traj)}")
print(f"    TV(u): {total_variation(u_traj):.2f}")

# Step 3: Compute v_eps for decreasing eps
print("\n[3] Computing v_eps for various eps...")
epsilons = [10.0, 1.0, 0.1, 0.01, 0.001, 0.0001]
results = []

for eps in epsilons:
    v_eps = minimize_f_epsilon(u_traj, attractor, eps, dt)
    fid, tv, total = f_epsilon(u_traj, v_eps, eps, dt)
    dist_to_u = np.sqrt(np.mean(np.sum((u_traj - v_eps)**2, axis=1)))
    results.append({
        'eps': eps,
        'TV': tv,
        'fidelity': fid,
        'total': total,
        'dist_to_u': dist_to_u
    })
    print(f"    eps={eps:10.4f}  |  TV(v_eps)={tv:12.2f}  |  "
          f"(1/eps)*Fid={fid:12.2f}  |  ||u-v||_avg={dist_to_u:.4f}")

# Step 4: Analyze DS3
print("\n" + "=" * 60)
print("DS3 ANALYSIS")
print("=" * 60)

tvs = [r['TV'] for r in results]
eps_vals = [r['eps'] for r in results]

# Check if TV is bounded
tv_ratio = tvs[-1] / tvs[-2] if tvs[-2] > 0 else float('inf')
print(f"\n  TV at smallest eps ({eps_vals[-1]}): {tvs[-1]:.2f}")
print(f"  TV at largest eps ({eps_vals[0]}):  {tvs[0]:.2f}")
print(f"  TV ratio (last/second-last):    {tv_ratio:.2f}")

if tvs[-1] < 2 * tvs[0]:
    print("\n  >> DS3 HOLDS: TV(v_eps) appears BOUNDED as eps -> 0")
    print("  >> The Lorenz attractor has sufficient compactness for DS3.")
elif tvs[-1] < 10 * tvs[0]:
    print("\n  >> DS3 MARGINAL: TV(v_eps) grows moderately")
    print("  >> Precompactness may hold with additional regularity.")
else:
    print("\n  >> DS3 FAILS: TV(v_eps) grows rapidly as eps -> 0")
    print("  >> This attractor lacks DS3-precompactness.")

# Step 5: Plot (if matplotlib available)
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot 1: TV vs eps
    axes[0].loglog(eps_vals, tvs, 'ro-', linewidth=2, markersize=8)
    axes[0].set_xlabel('epsilon')
    axes[0].set_ylabel('TV(v_epsilon)')
    axes[0].set_title('Total Variation vs epsilon')
    axes[0].grid(True, alpha=0.3)
    axes[0].axhline(y=total_variation(u_traj), color='blue', linestyle='--',
                     label='TV(u)')
    axes[0].legend()

    # Plot 2: Distance to u vs eps
    dists = [r['dist_to_u'] for r in results]
    axes[1].loglog(eps_vals, dists, 'bs-', linewidth=2, markersize=8)
    axes[1].set_xlabel('epsilon')
    axes[1].set_ylabel('||u - v_epsilon||_avg')
    axes[1].set_title('Fidelity vs epsilon')
    axes[1].grid(True, alpha=0.3)

    # Plot 3: Lorenz attractor with selection
    ax3 = fig.add_subplot(133, projection='3d')
    ax3.plot(attractor[::10, 0], attractor[::10, 1], attractor[::10, 2],
             'b-', alpha=0.2, linewidth=0.5, label='Attractor')
    # Show v_eps for eps=0.01
    v_mid = minimize_f_epsilon(u_traj[:500], attractor, 0.01, dt)
    ax3.plot(v_mid[:, 0], v_mid[:, 1], v_mid[:, 2],
             'r-', linewidth=1, label='v_{0.01}')
    ax3.set_title('Selection on Lorenz')
    ax3.legend(fontsize=8)

    plt.tight_layout()
    outpath = r'C:\Users\User\OneDrive\.RESEARCH\Natur&Technik\3 Folgebeweise\Navier-Stokes\compute_ds3_lorenz.png'
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\n  Plot saved to: {outpath}")
except ImportError:
    print("\n  (matplotlib not available, skipping plot)")

print("\n[DONE]")
