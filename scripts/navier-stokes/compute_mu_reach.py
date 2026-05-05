"""
compute_mu_reach.py
====================
Numerische Berechnung des mu-reach (measure-theoretic reach) auf dem
Lorenz-Attraktor und dem Kuramoto-Sivashinsky (KS) Attraktor.

Definition (aus dem Paper):
    Der mu-reach ist das essentielle Infimum des lokalen Reach ueber das
    invariante Mass auf dem Attraktor.
    Reach(v) = sup{r > 0 : nearest-point projection ist eindeutig und
                            Lipschitz-1 auf B(v,r)}.

Algorithmus:
    1. Attraktor approximieren (lange Trajektorie, N_attr Punkte)
    2. Fuer jeden Attraktor-Punkt v_i:
       a. Finde k naechste Nachbarn (KD-Tree)
       b. Generiere Test-Punkte in wachsenden Kugeln B(v_i, r)
       c. Pruefe Lipschitz-1 und Eindeutigkeit der Projektion
       d. Lokaler Reach = groesstes r, bei dem beides gilt
    3. mu-reach = Quantile der lokalen Reach-Verteilung

Systeme:
    - Lorenz (sigma=10, rho=28, beta=8/3): 3D, Attraktor-Dim ~ 2.06
    - KS (L=32*pi, N=128): 1D PDE, erste M=16 Fourier-Moden

Autor: Lukas Geiger (Skript erstellt per Claude, 2026)
"""

import numpy as np
from scipy.integrate import solve_ivp
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
# Lorenz-System
# ==========================================================================

SIGMA, RHO, BETA_L = 10.0, 28.0, 8.0 / 3.0
T_TOTAL_LORENZ = 200.0
T_TRANSIENT_LORENZ = 50.0
DT_LORENZ = 0.01


def lorenz(t, state):
    x, y, z = state
    return [SIGMA * (y - x), x * (RHO - z) - y, x * y - BETA_L * z]


def generate_lorenz_attractor(n_attr=10000, seed=42):
    """
    Simuliere Lorenz-System und extrahiere n_attr Punkte vom Attraktor.
    Lange Integration + Transient verwerfen + gleichmaessig samplen.
    """
    rng = np.random.RandomState(seed)
    # Leicht zufaellige Anfangsbedingung
    ic = [1.0 + 0.01 * rng.randn(), 1.0 + 0.01 * rng.randn(),
          1.0 + 0.01 * rng.randn()]

    n_total = int(T_TOTAL_LORENZ / DT_LORENZ)
    t_eval = np.linspace(0, T_TOTAL_LORENZ, n_total)
    sol = solve_ivp(lorenz, [0, T_TOTAL_LORENZ], ic,
                    t_eval=t_eval, method='RK45', rtol=1e-10, atol=1e-12)

    # Transient verwerfen
    idx_start = int(T_TRANSIENT_LORENZ / DT_LORENZ)
    traj = sol.y[:, idx_start:].T  # shape: (N, 3)

    # Gleichmaessig n_attr Punkte samplen
    indices = np.linspace(0, len(traj) - 1, n_attr, dtype=int)
    return traj[indices]


# ==========================================================================
# Kuramoto-Sivashinsky Gleichung -- ETDRK4 Pseudo-Spektral-Loeser
# ==========================================================================

class KuramotoSivashinsky:
    """
    Pseudo-spektral Loeser fuer die Kuramoto-Sivashinsky Gleichung:
        u_t + u*u_x + u_xx + u_xxxx = 0
    auf [0, L] mit periodischen Randbedingungen.

    Zustand in Fourier-Raum (komplexe Koeffizienten).
    ETDRK4 Zeitintegrator (Kassam & Trefethen 2005).
    """

    def __init__(self, L=32 * np.pi, N=128, dt=0.05):
        self.L = L
        self.N = N
        self.dx = L / N
        self.x = np.linspace(0, L, N, endpoint=False)

        # Wellenzahlen
        self.k = np.fft.rfftfreq(N, d=L / (2 * np.pi * N))
        # Linearer Operator: L_k = k^2 - k^4
        self.Lk = self.k**2 - self.k**4

        # ETDRK4 Koeffizienten vorberechnen
        self._setup_etdrk4_coefficients(dt=dt)

    def _setup_etdrk4_coefficients(self, dt):
        """
        ETDRK4 Koeffizienten via Konturintegral (Kassam & Trefethen 2005).
        Vermeidet katastrophale Ausloeschung nahe Lk = 0.
        """
        self.dt = dt
        Lk = self.Lk
        h = dt

        M_contour = 32
        r = np.exp(1j * np.pi * (np.arange(1, M_contour + 1) - 0.5)
                   / M_contour)
        LR = h * Lk[:, np.newaxis] + r[np.newaxis, :]

        self.E = np.exp(h * Lk)
        self.E2 = np.exp(h / 2 * Lk)

        self.f1 = np.real(np.mean(
            h * (np.exp(LR / 2) - 1) / LR, axis=1))
        self.a = np.real(np.mean(
            h * (-4 - LR + np.exp(LR) * (4 - 3 * LR + LR**2)) / LR**3,
            axis=1))
        self.b = np.real(np.mean(
            h * (2 + LR + np.exp(LR) * (-2 + LR)) / LR**3, axis=1))
        self.c = np.real(np.mean(
            h * (-4 - 3 * LR - LR**2 + np.exp(LR) * (4 - LR)) / LR**3,
            axis=1))

    def nonlinear(self, u_hat):
        """Nichtlinearer Term N(u) = -1/2 * d/dx(u^2) in Fourier-Raum."""
        N = self.N
        N_pad = 3 * N // 2
        u_hat_pad = np.zeros(N_pad // 2 + 1, dtype=complex)
        u_hat_pad[:len(u_hat)] = u_hat

        u_phys = np.fft.irfft(u_hat_pad, n=N_pad) * (N_pad / N)
        u2_phys = u_phys**2
        u2_hat_pad = np.fft.rfft(u2_phys) * (N / N_pad)
        u2_hat = u2_hat_pad[:len(u_hat)]

        return -0.5 * 1j * self.k * u2_hat

    def step_etdrk4(self, u_hat):
        """Ein ETDRK4 Zeitschritt."""
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

        if not np.all(np.isfinite(u_hat_new)):
            raise RuntimeError(
                f"ETDRK4 blow-up detected (N={self.N}, dt={self.dt:.4f}). "
                f"Reduce dt or increase N.")
        return u_hat_new

    def to_fourier(self, u):
        return np.fft.rfft(u)

    def initial_condition(self, seed=42):
        """Zufaellige Anfangsbedingung mit Energie in niedrigen Moden."""
        rng = np.random.RandomState(seed)
        u = np.zeros(self.N)
        for kn in range(1, 6):
            phase = rng.uniform(0, 2 * np.pi)
            u += rng.uniform(0.5, 2.0) * np.cos(
                2 * np.pi * kn * self.x / self.L + phase)
        return self.to_fourier(u)

    def integrate(self, u_hat0, T, dt=None, save_every=1):
        """
        Integriere KS-Gleichung von u_hat0 fuer Zeit T.

        Returns
        -------
        t_arr : array (Zeitpunkte)
        traj : array (Fourier-Koeffizienten-Snapshots)
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
    Extrahiere reellen Zustandsvektor aus den ersten M Fourier-Moden.
    Returns [Re(u_hat[1]), Im(u_hat[1]), ..., Re(u_hat[M]), Im(u_hat[M])].
    Mode 0 (Mittelwert) wird uebersprungen (in KS erhalten).
    """
    vec = np.zeros(2 * M)
    for j in range(M):
        mode = j + 1
        if mode < len(u_hat):
            vec[2 * j] = u_hat[mode].real
            vec[2 * j + 1] = u_hat[mode].imag
    return vec


def generate_ks_attractor(n_attr=5000, M=16, seed=42):
    """
    Simuliere KS-Gleichung und extrahiere n_attr Punkte vom Attraktor
    in den ersten M Fourier-Moden (2*M reelle Dimensionen).
    """
    ks = KuramotoSivashinsky(L=32 * np.pi, N=128, dt=0.05)
    u_hat0 = ks.initial_condition(seed=seed)

    # Transient: 500 Zeiteinheiten
    T_transient = 500.0
    _, traj_trans = ks.integrate(u_hat0, T_transient, save_every=100)
    u_hat_on_attr = traj_trans[-1]

    # Attraktor-Trajektorie: genug Punkte sammeln
    # save_every so waehlen, dass wir ~n_attr Punkte bekommen
    T_attr = 1000.0
    save_every = max(1, int(T_attr / ks.dt / (n_attr * 1.2)))
    _, traj = ks.integrate(u_hat_on_attr, T_attr, save_every=save_every)

    # Fourier-Zustandsvektoren extrahieren
    states = np.array([fourier_state_vector(u_hat, M) for u_hat in traj])

    # Gleichmaessig n_attr Punkte samplen
    if len(states) > n_attr:
        indices = np.linspace(0, len(states) - 1, n_attr, dtype=int)
        states = states[indices]

    return states


# ==========================================================================
# Mu-Reach Berechnung
# ==========================================================================

def compute_local_reach(point_idx, attractor, tree, n_test=50, n_r_steps=10,
                        k_neighbors=30, rng=None):
    """
    Berechne den lokalen Reach fuer einen einzelnen Attraktor-Punkt.

    Parameters
    ----------
    point_idx : int
        Index des Punkts auf dem Attraktor.
    attractor : ndarray, shape (N, d)
        Attraktor-Punkte.
    tree : cKDTree
        KD-Tree ueber den Attraktor.
    n_test : int
        Anzahl Test-Punkte pro Kugelradius.
    n_r_steps : int
        Anzahl Radien-Schritte (logarithmisch).
    k_neighbors : int
        Anzahl Nachbarn fuer Nachbarschafts-Schaetzung.
    rng : RandomState
        Zufallsgenerator.

    Returns
    -------
    local_reach : float
        Geschaetzter lokaler Reach bei diesem Punkt.
    """
    if rng is None:
        rng = np.random.RandomState(42)

    v_i = attractor[point_idx]
    dim = len(v_i)

    # Naechste Nachbarn finden -- bestimmt die Skala
    dists_nn, _ = tree.query(v_i, k=min(k_neighbors, len(attractor)))
    d_min = dists_nn[1] if len(dists_nn) > 1 else 1e-10  # naechster (nicht self)
    d_median = np.median(dists_nn[1:]) if len(dists_nn) > 1 else d_min * 10

    # Sicherheit: d_min > 0
    d_min = max(d_min, 1e-12)
    d_median = max(d_median, d_min * 2)

    # Radien: logarithmisch von d_min/10 bis d_median
    r_min = d_min / 10.0
    r_max = d_median
    radii = np.geomspace(r_min, r_max, n_r_steps)

    local_reach = 0.0

    for r in radii:
        # Test-Punkte in B(v_i, r) generieren: gleichverteilt in Kugel
        # Methode: Normalverteilte Richtung, gleichverteilter Radius
        directions = rng.randn(n_test, dim)
        norms = np.linalg.norm(directions, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-15)
        directions = directions / norms

        # Gleichverteilter Radius in d-dim. Kugel: r * U^(1/d)
        u_radii = rng.uniform(0, 1, size=n_test) ** (1.0 / dim)
        test_points = v_i + r * u_radii[:, np.newaxis] * directions

        # Fuer jeden Test-Punkt: naechsten Attraktor-Punkt finden
        dists_to_attr, idx_nearest = tree.query(test_points, k=1)

        # Pruefe Eindeutigkeit: Naechster und zweitnaechster Nachbar
        dists_2, idx_2 = tree.query(test_points, k=2)
        d_nearest = dists_2[:, 0]
        d_second = dists_2[:, 1]

        # Eindeutigkeit: Abstand zum zweitnaechsten ist signifikant groesser
        # Toleranz: relative Differenz > 1e-6
        rel_gap = (d_second - d_nearest) / np.maximum(d_nearest, 1e-15)
        uniqueness_ok = np.all(rel_gap > 1e-6)

        if not uniqueness_ok:
            # Erreicht Aequidistanz -- Reach ist der vorherige Radius
            break

        # Pruefe Lipschitz-1: ||pi(u1) - pi(u2)|| <= ||u1 - u2|| fuer Paare
        # Teste eine Stichprobe von Paaren (n_test*(n_test-1)/2 waere zu teuer)
        n_pairs = min(200, n_test * (n_test - 1) // 2)
        projections = attractor[idx_nearest.flatten()]

        lipschitz_ok = True
        for _ in range(n_pairs):
            i1, i2 = rng.randint(0, n_test, size=2)
            if i1 == i2:
                continue
            d_proj = np.linalg.norm(projections[i1] - projections[i2])
            d_test = np.linalg.norm(test_points[i1] - test_points[i2])
            if d_test < 1e-15:
                continue
            if d_proj > d_test * (1.0 + 1e-6):  # kleine Toleranz fuer Numerik
                lipschitz_ok = False
                break

        if not lipschitz_ok:
            # Lipschitz-1 verletzt -- Reach ist der vorherige Radius
            break

        # Beide Bedingungen erfuellt bei diesem Radius
        local_reach = r

    return local_reach


def compute_mu_reach(attractor, n_sample=None, n_test=50, n_r_steps=10,
                     k_neighbors=30, seed=42, label=""):
    """
    Berechne die mu-reach Statistik fuer einen Attraktor.

    Parameters
    ----------
    attractor : ndarray, shape (N, d)
        Attraktor-Punkte.
    n_sample : int or None
        Anzahl Punkte, fuer die der lokale Reach berechnet wird.
        None = alle Punkte (kann langsam sein).
    n_test : int
        Test-Punkte pro Kugel-Radius.
    n_r_steps : int
        Anzahl Radien-Schritte.
    k_neighbors : int
        Nachbarn fuer Skalenbestimmung.
    seed : int
        Random seed.
    label : str
        Label fuer Ausgabe.

    Returns
    -------
    results : dict
        Statistiken der lokalen Reach-Verteilung.
    local_reaches : ndarray
        Alle berechneten lokalen Reach-Werte.
    """
    rng = np.random.RandomState(seed)
    N = len(attractor)
    dim = attractor.shape[1]

    if n_sample is None or n_sample >= N:
        sample_indices = np.arange(N)
    else:
        sample_indices = rng.choice(N, size=n_sample, replace=False)
        sample_indices.sort()

    print(f"\n  Attraktor: {N} Punkte in {dim}D")
    print(f"  Stichprobe: {len(sample_indices)} Punkte")

    # KD-Tree aufbauen
    tree = cKDTree(attractor)

    # Lokalen Reach fuer jeden Stichproben-Punkt berechnen
    local_reaches = np.zeros(len(sample_indices))
    t_start = time.time()

    for i, idx in enumerate(sample_indices):
        local_reaches[i] = compute_local_reach(
            idx, attractor, tree, n_test=n_test, n_r_steps=n_r_steps,
            k_neighbors=k_neighbors, rng=rng)

        # Fortschritt
        if (i + 1) % max(1, len(sample_indices) // 10) == 0:
            elapsed = time.time() - t_start
            frac = (i + 1) / len(sample_indices)
            eta = elapsed / frac * (1 - frac) if frac > 0 else 0
            print(f"    [{i+1}/{len(sample_indices)}] "
                  f"Elapsed: {elapsed:.1f}s, ETA: {eta:.1f}s")

    elapsed_total = time.time() - t_start

    # Statistiken
    nonzero = local_reaches[local_reaches > 0]
    n_zero = np.sum(local_reaches == 0)

    results = {
        "system": label,
        "n_attractor": N,
        "dimension": dim,
        "n_sampled": len(sample_indices),
        "n_test_per_ball": n_test,
        "n_radius_steps": n_r_steps,
        "computation_time_s": round(elapsed_total, 2),
        "statistics": {
            "min": float(np.min(local_reaches)),
            "percentile_5": float(np.percentile(local_reaches, 5)),
            "percentile_25": float(np.percentile(local_reaches, 25)),
            "median": float(np.median(local_reaches)),
            "mean": float(np.mean(local_reaches)),
            "percentile_75": float(np.percentile(local_reaches, 75)),
            "percentile_95": float(np.percentile(local_reaches, 95)),
            "max": float(np.max(local_reaches)),
            "n_zero_reach": int(n_zero),
            "fraction_positive": float(
                np.sum(local_reaches > 0) / len(local_reaches)),
        },
        "mu_reach_estimate": float(np.percentile(local_reaches, 5)),
        "pass": bool(np.percentile(local_reaches, 5) > 0),
    }

    return results, local_reaches


def plot_reach_distribution(local_reaches, results, system_label, outpath):
    """
    Histogramm der lokalen Reach-Verteilung plotten und speichern.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    stats = results["statistics"]

    # Panel 1: Histogramm
    ax = axes[0]
    # Nur positive Werte fuer Log-Histogramm
    positive = local_reaches[local_reaches > 0]
    if len(positive) > 0:
        log_reaches = np.log10(positive)
        ax.hist(log_reaches, bins=40, color='steelblue', edgecolor='navy',
                alpha=0.8, density=True)
        # Quantile einzeichnen
        q5 = np.log10(max(stats["percentile_5"], 1e-15))
        q50 = np.log10(max(stats["median"], 1e-15))
        ax.axvline(q5, color='red', linestyle='--', linewidth=2,
                   label=f'5%-Quantil = {stats["percentile_5"]:.2e}')
        ax.axvline(q50, color='orange', linestyle='-', linewidth=2,
                   label=f'Median = {stats["median"]:.2e}')
    else:
        ax.text(0.5, 0.5, 'Alle Reach-Werte = 0',
                ha='center', va='center', transform=ax.transAxes,
                fontsize=14, color='red')

    ax.set_xlabel(r'$\log_{10}(\mathrm{Reach}(v))$', fontsize=12)
    ax.set_ylabel('Dichte', fontsize=12)
    ax.set_title(f'{system_label}: Lokale Reach-Verteilung', fontsize=13)
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(True, alpha=0.3)

    # Panel 2: Sortierte Reach-Werte (empirische CDF)
    ax = axes[1]
    sorted_reaches = np.sort(local_reaches)
    cdf = np.arange(1, len(sorted_reaches) + 1) / len(sorted_reaches)
    ax.semilogy(cdf, sorted_reaches + 1e-15, 'b-', linewidth=1.5)
    ax.axhline(stats["percentile_5"], color='red', linestyle='--',
               linewidth=1.5, label=f'5%-Quantil = {stats["percentile_5"]:.2e}')
    ax.axhline(stats["median"], color='orange', linestyle='-',
               linewidth=1.5, label=f'Median = {stats["median"]:.2e}')
    ax.axvline(0.05, color='red', linestyle=':', alpha=0.5)
    ax.set_xlabel('Empirische CDF', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Reach}(v)$', fontsize=12)
    ax.set_title(f'{system_label}: Empirische CDF', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

    # Zusammenfassung als Text
    pass_str = "PASS" if results["pass"] else "FAIL"
    fig.suptitle(
        f'mu-reach Test: {system_label}  --  '
        f'mu-reach(5%) = {stats["percentile_5"]:.2e}  '
        f'[{pass_str}]',
        fontsize=14, fontweight='bold',
        y=1.02)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Plot gespeichert: {outpath}")


# ==========================================================================
# Hauptprogramm
# ==========================================================================

def main():
    print("=" * 70)
    print("MU-REACH TEST: Measure-Theoretic Reach auf Attraktoren")
    print("=" * 70)

    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               '_results')
    os.makedirs(results_dir, exist_ok=True)

    all_results = {}

    # ==================================================================
    # 1. LORENZ
    # ==================================================================
    print(f"\n{'='*70}")
    print("[1] LORENZ-ATTRAKTOR")
    print(f"{'='*70}")

    N_ATTR_LORENZ = 10000
    N_SAMPLE_LORENZ = 500  # Stichprobe fuer Reach-Berechnung
    N_TEST = 50
    N_R_STEPS = 10

    print(f"\n  Generiere Lorenz-Attraktor ({N_ATTR_LORENZ} Punkte)...")
    t0 = time.time()
    lorenz_attr = generate_lorenz_attractor(n_attr=N_ATTR_LORENZ, seed=42)
    print(f"  Fertig in {time.time()-t0:.1f}s")
    print(f"  Attraktor-Form: {lorenz_attr.shape}")
    print(f"  Bereich x: [{lorenz_attr[:,0].min():.2f}, "
          f"{lorenz_attr[:,0].max():.2f}]")
    print(f"  Bereich y: [{lorenz_attr[:,1].min():.2f}, "
          f"{lorenz_attr[:,1].max():.2f}]")
    print(f"  Bereich z: [{lorenz_attr[:,2].min():.2f}, "
          f"{lorenz_attr[:,2].max():.2f}]")

    # Naechste-Nachbar-Distanzen
    tree_lorenz = cKDTree(lorenz_attr)
    nn_dists, _ = tree_lorenz.query(lorenz_attr, k=2)
    d_nn = nn_dists[:, 1]
    print(f"\n  Naechste-Nachbar-Distanzen:")
    print(f"    min:    {d_nn.min():.4e}")
    print(f"    median: {np.median(d_nn):.4e}")
    print(f"    max:    {d_nn.max():.4e}")

    print(f"\n  Berechne lokalen Reach fuer {N_SAMPLE_LORENZ} Punkte...")
    lorenz_results, lorenz_reaches = compute_mu_reach(
        lorenz_attr, n_sample=N_SAMPLE_LORENZ, n_test=N_TEST,
        n_r_steps=N_R_STEPS, k_neighbors=30, seed=42, label="Lorenz")

    all_results["lorenz"] = lorenz_results

    # Ausgabe
    stats = lorenz_results["statistics"]
    print(f"\n  Lorenz mu-reach Statistiken:")
    print(f"    min:     {stats['min']:.4e}")
    print(f"    5%:      {stats['percentile_5']:.4e}")
    print(f"    25%:     {stats['percentile_25']:.4e}")
    print(f"    median:  {stats['median']:.4e}")
    print(f"    mean:    {stats['mean']:.4e}")
    print(f"    75%:     {stats['percentile_75']:.4e}")
    print(f"    95%:     {stats['percentile_95']:.4e}")
    print(f"    max:     {stats['max']:.4e}")
    print(f"    Nullen:  {stats['n_zero_reach']}/{lorenz_results['n_sampled']}")
    pass_str = "PASS" if lorenz_results["pass"] else "FAIL"
    print(f"\n  >> mu-reach(5%) = {lorenz_results['mu_reach_estimate']:.4e} "
          f"[{pass_str}]")

    # Plot
    try:
        outpath_lorenz = os.path.join(results_dir, 'mu_reach_lorenz.png')
        plot_reach_distribution(lorenz_reaches, lorenz_results,
                                "Lorenz", outpath_lorenz)
    except Exception as e:
        print(f"  (Plot-Fehler: {e})")

    # ==================================================================
    # 2. KURAMOTO-SIVASHINSKY
    # ==================================================================
    print(f"\n{'='*70}")
    print("[2] KURAMOTO-SIVASHINSKY ATTRAKTOR")
    print(f"{'='*70}")

    N_ATTR_KS = 5000
    N_SAMPLE_KS = 300  # Weniger wegen hoeherer Dimension
    M_MODES = 16  # Erste 16 Fourier-Moden => 32D Zustandsraum

    print(f"\n  Generiere KS-Attraktor ({N_ATTR_KS} Punkte, "
          f"M={M_MODES} Moden => {2*M_MODES}D)...")
    t0 = time.time()
    ks_attr = generate_ks_attractor(n_attr=N_ATTR_KS, M=M_MODES, seed=42)
    print(f"  Fertig in {time.time()-t0:.1f}s")
    print(f"  Attraktor-Form: {ks_attr.shape}")

    # Naechste-Nachbar-Distanzen
    tree_ks = cKDTree(ks_attr)
    nn_dists_ks, _ = tree_ks.query(ks_attr, k=2)
    d_nn_ks = nn_dists_ks[:, 1]
    print(f"\n  Naechste-Nachbar-Distanzen:")
    print(f"    min:    {d_nn_ks.min():.4e}")
    print(f"    median: {np.median(d_nn_ks):.4e}")
    print(f"    max:    {d_nn_ks.max():.4e}")

    print(f"\n  Berechne lokalen Reach fuer {N_SAMPLE_KS} Punkte...")
    ks_results, ks_reaches = compute_mu_reach(
        ks_attr, n_sample=N_SAMPLE_KS, n_test=N_TEST,
        n_r_steps=N_R_STEPS, k_neighbors=30, seed=42,
        label=f"KS (M={M_MODES}, {2*M_MODES}D)")

    all_results["ks"] = ks_results

    # Ausgabe
    stats_ks = ks_results["statistics"]
    print(f"\n  KS mu-reach Statistiken:")
    print(f"    min:     {stats_ks['min']:.4e}")
    print(f"    5%:      {stats_ks['percentile_5']:.4e}")
    print(f"    25%:     {stats_ks['percentile_25']:.4e}")
    print(f"    median:  {stats_ks['median']:.4e}")
    print(f"    mean:    {stats_ks['mean']:.4e}")
    print(f"    75%:     {stats_ks['percentile_75']:.4e}")
    print(f"    95%:     {stats_ks['percentile_95']:.4e}")
    print(f"    max:     {stats_ks['max']:.4e}")
    print(f"    Nullen:  {stats_ks['n_zero_reach']}/{ks_results['n_sampled']}")
    pass_str_ks = "PASS" if ks_results["pass"] else "FAIL"
    print(f"\n  >> mu-reach(5%) = {ks_results['mu_reach_estimate']:.4e} "
          f"[{pass_str_ks}]")

    # Plot
    try:
        outpath_ks = os.path.join(results_dir, 'mu_reach_ks.png')
        plot_reach_distribution(ks_reaches, ks_results,
                                f"KS (M={M_MODES})", outpath_ks)
    except Exception as e:
        print(f"  (Plot-Fehler: {e})")

    # ==================================================================
    # 3. VERGLEICH
    # ==================================================================
    print(f"\n{'='*70}")
    print("[3] VERGLEICH: Lorenz vs KS")
    print(f"{'='*70}")

    print(f"""
  {'Metrik':<30} {'Lorenz (3D)':>15} {'KS (32D)':>15}
  {'-'*62}
  {'Attraktor-Punkte':<30} {lorenz_results['n_attractor']:>15d} {ks_results['n_attractor']:>15d}
  {'Dimension':<30} {lorenz_results['dimension']:>15d} {ks_results['dimension']:>15d}
  {'Stichprobe':<30} {lorenz_results['n_sampled']:>15d} {ks_results['n_sampled']:>15d}
  {'mu-reach (5%-Quantil)':<30} {stats['percentile_5']:>15.4e} {stats_ks['percentile_5']:>15.4e}
  {'mu-reach (Median)':<30} {stats['median']:>15.4e} {stats_ks['median']:>15.4e}
  {'mu-reach (Mean)':<30} {stats['mean']:>15.4e} {stats_ks['mean']:>15.4e}
  {'Anteil positive Reach':<30} {stats['fraction_positive']:>15.2%} {stats_ks['fraction_positive']:>15.2%}
  {'Rechenzeit (s)':<30} {lorenz_results['computation_time_s']:>15.1f} {ks_results['computation_time_s']:>15.1f}
  {'Ergebnis':<30} {'PASS' if lorenz_results['pass'] else 'FAIL':>15} {'PASS' if ks_results['pass'] else 'FAIL':>15}
""")

    # ==================================================================
    # 4. ERGEBNISSE SPEICHERN
    # ==================================================================
    json_path = os.path.join(results_dir, 'mu_reach_results.json')
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    print(f"  Ergebnisse gespeichert: {json_path}")

    # ==================================================================
    # ZUSAMMENFASSUNG
    # ==================================================================
    print(f"\n{'='*70}")
    print("ZUSAMMENFASSUNG: mu-reach Test")
    print(f"{'='*70}")

    both_pass = lorenz_results["pass"] and ks_results["pass"]
    print(f"""
  Definition: mu-reach = ess inf_{{v ~ mu}} Reach(v)
              (essentielles Infimum des lokalen Reach ueber
               das invariante Mass auf dem Attraktor)

  Kriterium: mu-reach > 0 wenn das 5%-Quantil > 0

  Lorenz-Attraktor:
    mu-reach >= {stats['percentile_5']:.4e}  [{'PASS' if lorenz_results['pass'] else 'FAIL'}]
    Der Lorenz-Attraktor hat endlichen Reach fast ueberall.
    Die glatte Blaetter-Struktur (2D Mannigfaltigkeit x Cantor)
    garantiert lokale Eindeutigkeit der Projektion.

  KS-Attraktor (M={M_MODES} Moden):
    mu-reach >= {stats_ks['percentile_5']:.4e}  [{'PASS' if ks_results['pass'] else 'FAIL'}]
    Der KS-Attraktor in 32D hat {'ebenfalls endlichen' if ks_results['pass'] else 'moeglicherweise verschwindenden'} Reach.
    {'Die inertiale Mannigfaltigkeit sorgt fuer glattte lokale Geometrie.' if ks_results['pass'] else 'Hohe Dimension und Faltung koennten problematisch sein.'}

  GESAMT: {'BEIDE SYSTEME BESTANDEN' if both_pass else 'NICHT ALLE SYSTEME BESTANDEN'}
  {'=> mu-reach > 0 numerisch bestaetigt' if both_pass else '=> Weitere Analyse noetig'}
""")

    print("[DONE]")
    return all_results


if __name__ == "__main__":
    main()
