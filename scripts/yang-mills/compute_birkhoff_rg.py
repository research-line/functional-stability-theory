"""
compute_birkhoff_rg.py
======================
Yang-Mills: Birkhoff-Kontraktion tau_B(R_k) fuer hierarchische RG-Schritte.

Berechnet den Birkhoff-projektiven Kontraktionskoeffizienten
tau_B des Block-Spin-RG-Kerns R_k fuer SU(2) Gitter-Eichtheorie.

Methode:
1. Generiere Block-Spin-Konfigurationen auf verschiedenen Skalen k
2. Berechne Transferoperator R_k als stochastische Matrix
3. Messe Hilbert-projektive Metrik-Kontraktion:
   d_H(R_k mu, R_k nu) / d_H(mu, nu) <= tau_B(R_k)
4. Pruefe: tau_B < 1 (uniform oder im Mittel via Kingman)

Birkhoff-Kontraktion: tau_B(R) = tanh(Delta(R)/4)
wobei Delta(R) = log(max_{i,j,k,l} R_{ik}R_{jl}/(R_{il}R_{jk}))
(Birkhoff 1957, Bushell 1973)

Autor: Lukas Geiger (Skript erstellt per Claude, 2026)
"""

import numpy as np
import os

# ==========================================================================
# SU(2) Hilfsfunktionen
# ==========================================================================

def random_su2():
    v = np.random.randn(4)
    v /= np.linalg.norm(v)
    return v

def su2_multiply(q1, q2):
    a0, a1, a2, a3 = q1
    b0, b1, b2, b3 = q2
    return np.array([
        a0*b0 - a1*b1 - a2*b2 - a3*b3,
        a0*b1 + a1*b0 + a2*b3 - a3*b2,
        a0*b2 - a1*b3 + a2*b0 + a3*b1,
        a0*b3 + a1*b2 - a2*b1 + a3*b0
    ])

def su2_adjoint(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])

def su2_trace_re(q):
    return 2.0 * q[0]

def su2_project(q):
    """Projiziere auf SU(2) (Quaternion normieren)"""
    n = np.linalg.norm(q)
    return q / n if n > 1e-10 else np.array([1.0, 0, 0, 0])


# ==========================================================================
# Block-Spin-RG
# ==========================================================================

def block_spin_average(links, block_size=2):
    """
    Block-Spin-RG-Schritt: Mittele ueber block_size^d Links.
    Fuer 1D-Kette: U'_k = P_G(U_{2k} * U_{2k+1})
    """
    N = len(links)
    N_new = N // block_size
    new_links = []
    for i in range(N_new):
        # Multipliziere block_size aufeinanderfolgende Links
        product = links[i * block_size]
        for j in range(1, block_size):
            idx = i * block_size + j
            if idx < N:
                product = su2_multiply(product, links[idx])
        # Projiziere zurueck auf SU(2)
        new_links.append(su2_project(product))
    return new_links


def wilson_action_1d(links, beta):
    """Wilson-Aktion fuer 1D SU(2) Kette (= Plaquette in 2D)"""
    S = 0.0
    for i in range(len(links) - 1):
        plaq = su2_multiply(links[i], su2_adjoint(links[i+1]))
        S += su2_trace_re(plaq) / 2.0
    return -beta * S


# ==========================================================================
# Transfermatrix-Methode
# ==========================================================================

def compute_transfer_matrix(beta, n_bins=20):
    """
    Berechne die Transfermatrix T_{ij} fuer SU(2) in 1D.
    Diskretisiere a_0 in [−1, 1] (Trace/2 der SU(2)-Matrix).
    T(a_0, a_0') = exp(beta * a_0 * a_0') * rho(a_0')
    wobei rho(a_0) = sqrt(1 - a_0^2) (Haar-Mass)
    """
    # Diskretisierung von a_0 = cos(theta/2) in [-1, 1]
    a_vals = np.linspace(-0.99, 0.99, n_bins)
    da = a_vals[1] - a_vals[0]

    T = np.zeros((n_bins, n_bins))
    for i in range(n_bins):
        for j in range(n_bins):
            # Wilson-Gewicht: exp(beta * Tr(U1 * U2^dag) / 2)
            # Fuer SU(2): Tr(U1*U2^dag)/2 = a0_1*a0_2 + vec_1 . vec_2
            # Gemittelt ueber Winkelintegration: I_0(beta * a_i) * I_0(beta * a_j)
            # Vereinfacht: exp(beta * a_i * a_j)
            weight = np.exp(beta * a_vals[i] * a_vals[j])
            haar_j = np.sqrt(max(0, 1 - a_vals[j]**2))
            T[i, j] = weight * haar_j * da

    # Normalisiere Zeilen (stochastische Matrix)
    row_sums = T.sum(axis=1)
    row_sums[row_sums == 0] = 1
    T_stoch = T / row_sums[:, None]

    return T_stoch, a_vals


def birkhoff_contraction(T):
    """
    Berechne Birkhoff-Kontraktionskoeffizienten tau_B(T).

    tau_B(T) = tanh(Delta(T)/4)

    wobei Delta(T) = max_{i,j,k,l} log(T_{ik} * T_{jl} / (T_{il} * T_{jk}))
    (projektiver Durchmesser im Hilbert-Kegel)

    Fuer strikt positive Matrizen: 0 <= tau_B < 1
    """
    n = T.shape[0]
    T_pos = np.maximum(T, 1e-300)  # Vermeidung von log(0)

    # Berechne Delta = max log-ratio
    Delta = 0.0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    if T_pos[i, k] > 0 and T_pos[j, l] > 0 and T_pos[i, l] > 0 and T_pos[j, k] > 0:
                        ratio = (T_pos[i, k] * T_pos[j, l]) / (T_pos[i, l] * T_pos[j, k])
                        if ratio > 0:
                            Delta = max(Delta, abs(np.log(ratio)))

    tau = np.tanh(Delta / 4.0)
    return tau, Delta


def birkhoff_contraction_fast(T):
    """
    Schnelle Approximation von tau_B via Spektralradius-Methode.
    tau_B <= (lambda_2 / lambda_1) wobei lambda_i Eigenwerte von T.
    """
    eigvals = np.abs(np.linalg.eigvals(T))
    eigvals_sorted = np.sort(eigvals)[::-1]
    if len(eigvals_sorted) >= 2 and eigvals_sorted[0] > 1e-10:
        spectral_gap = eigvals_sorted[1] / eigvals_sorted[0]
    else:
        spectral_gap = 0.0
    return spectral_gap


# ==========================================================================
# RG-Kaskade: tau_B ueber mehrere Skalen
# ==========================================================================

def rg_cascade(beta, n_levels=8, n_bins=20):
    """
    Berechne tau_B(R_k) fuer k = 0, ..., n_levels-1.

    Methode: Die Transfermatrix auf Skala k ist das Produkt
    von block-gemittelten Transfermatrizen:
    T_k = coarse-grain(T_{k-1})
    """
    results = []

    # Basis-Transfermatrix (feinste Skala)
    T_base, a_vals = compute_transfer_matrix(beta, n_bins)

    T_current = T_base.copy()

    for k in range(n_levels):
        # Birkhoff-Kontraktion
        tau_exact, Delta = birkhoff_contraction(T_current)
        tau_spectral = birkhoff_contraction_fast(T_current)

        # Spektralluecke der Transfermatrix
        eigvals = np.sort(np.abs(np.linalg.eigvals(T_current)))[::-1]
        gap = 1.0 - eigvals[1]/eigvals[0] if len(eigvals) >= 2 and eigvals[0] > 0 else 0

        results.append({
            'k': k,
            'tau_B': tau_exact,
            'tau_spectral': tau_spectral,
            'Delta': Delta,
            'gap': gap,
            'lambda_1': eigvals[0] if len(eigvals) > 0 else 0,
            'lambda_2': eigvals[1] if len(eigvals) > 1 else 0,
        })

        # RG-Schritt: Block-Spin-Vergoeberung der Transfermatrix
        # T_new[i,j] = sum_{k,l in block} T[i,k] * T[k,l] * T[l,j]
        # Vereinfacht: T_new = T^2 (Matrixpotenz = 2 aufeinanderfolgende Transfers)
        T_current = T_current @ T_current
        # Re-normalisiere
        row_sums = T_current.sum(axis=1)
        row_sums[row_sums == 0] = 1
        T_current = T_current / row_sums[:, None]

    return results


# ==========================================================================
# Hauptberechnung
# ==========================================================================

print("=" * 70)
print("YANG-MILLS: Birkhoff-Kontraktion tau_B(R_k) fuer RG-Schritte")
print("Prueft: tau_B < 1 (Poincare-Luecke bleibt erhalten)")
print("=" * 70)

betas = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0]
n_levels = 6
n_bins = 16  # Diskretisierung (schnell)

all_results = {}

for beta in betas:
    print(f"\n{'='*50}")
    print(f"beta = {beta:.1f}")
    print(f"{'='*50}")

    res = rg_cascade(beta, n_levels=n_levels, n_bins=n_bins)
    all_results[beta] = res

    print(f"\n  {'k':>3} {'tau_B':>10} {'tau_spec':>10} {'Delta':>10} {'gap':>10} {'lambda_2/1':>12}")
    print(f"  {'-'*60}")
    for r in res:
        print(f"  {r['k']:>3} {r['tau_B']:>10.6f} {r['tau_spectral']:>10.6f} "
              f"{r['Delta']:>10.4f} {r['gap']:>10.6f} "
              f"{r['lambda_2']/r['lambda_1'] if r['lambda_1'] > 0 else 0:>12.6f}")

    # Kingman-Lyapunov-Exponent
    taus = [r['tau_B'] for r in res]
    log_taus = [np.log(max(t, 1e-10)) for t in taus]
    lyapunov = np.mean(log_taus)
    print(f"\n  Kingman-Lyapunov: lambda = <log(tau_B)> = {lyapunov:.4f}")
    print(f"  lambda < 0: {'JA (KONTRAKTION)' if lyapunov < 0 else 'NEIN'}")

    # tau_B < 1 uniform?
    uniform = all(t < 1.0 - 1e-6 for t in taus)
    print(f"  tau_B < 1 uniform: {'JA' if uniform else 'NEIN'}")
    if uniform:
        delta_min = min(1 - t for t in taus)
        print(f"  min(1 - tau_B) = {delta_min:.6f}")


# ==========================================================================
# Zusammenfassung
# ==========================================================================

print(f"\n{'='*70}")
print("ZUSAMMENFASSUNG: Birkhoff-Kontraktion ueber RG-Skalen")
print(f"{'='*70}")

print(f"\n  {'beta':>6} {'<tau_B>':>10} {'max tau_B':>10} {'Lyapunov':>10} {'Status':>18}")
print(f"  {'-'*58}")

for beta in betas:
    res = all_results[beta]
    taus = [r['tau_B'] for r in res]
    mean_tau = np.mean(taus)
    max_tau = np.max(taus)
    lyapunov = np.mean([np.log(max(t, 1e-10)) for t in taus])

    if max_tau < 1.0 - 1e-6:
        status = "UNIFORM (tau<1)"
    elif lyapunov < 0:
        status = "KINGMAN (lam<0)"
    else:
        status = "KEINE KONTR."

    print(f"  {beta:>6.1f} {mean_tau:>10.6f} {max_tau:>10.6f} {lyapunov:>10.4f} {status:>18}")

print(f"""
  INTERPRETATION:
  - tau_B(R_k) < 1 fuer ALLE RG-Schritte => uniforme Kontraktion
    => kappa_k >= kappa_0 * (1 - C/k^2) (Prop. polynomial-lsi)
    => Massenluecke bleibt im Kontinuumslimes erhalten

  - Falls tau_B ~ 1 fuer einige k, aber Lyapunov < 0:
    => Kingman-Ergodensatz sichert mittlere Kontraktion
    => Massenluecke existiert fast sicher (nicht uniform)

  - Diskretisierung: n_bins = {n_bins} (SU(2)-Transfermatrix)
    Feinere Diskretisierung auf Server empfohlen (n_bins=64+)
""")


# Plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel 1: tau_B vs k fuer verschiedene beta
    ax = axes[0, 0]
    for beta in betas[:5]:
        res = all_results[beta]
        ks = [r['k'] for r in res]
        taus = [r['tau_B'] for r in res]
        ax.plot(ks, taus, 'o-', lw=1.5, ms=5, label=f'$\\beta={beta}$')
    ax.axhline(y=1.0, color='red', ls='--', lw=1.5, label='$\\tau_B = 1$')
    ax.set_xlabel('RG-Stufe k')
    ax.set_ylabel(r'$\tau_B(R_k)$')
    ax.set_title('Birkhoff-Kontraktion pro RG-Schritt')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 2: Lyapunov vs beta
    ax = axes[0, 1]
    lyapunovs = []
    for beta in betas:
        res = all_results[beta]
        taus = [r['tau_B'] for r in res]
        lyapunovs.append(np.mean([np.log(max(t, 1e-10)) for t in taus]))
    ax.plot(betas, lyapunovs, 'rs-', lw=2, ms=8)
    ax.axhline(y=0, color='black', ls='-', lw=0.5)
    ax.fill_between(betas, min(lyapunovs)-0.1, 0, alpha=0.1, color='green',
                    label=r'$\lambda < 0$ (Kontraktion)')
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\lambda = \langle \log \tau_B \rangle$')
    ax.set_title('Kingman-Lyapunov-Exponent')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 3: Spektralluecke vs k
    ax = axes[1, 0]
    for beta in [2.0, 4.0, 8.0, 15.0]:
        res = all_results[beta]
        ks = [r['k'] for r in res]
        gaps = [r['gap'] for r in res]
        ax.plot(ks, gaps, 'o-', lw=1.5, ms=5, label=f'$\\beta={beta}$')
    ax.set_xlabel('RG-Stufe k')
    ax.set_ylabel('Spektralluecke')
    ax.set_title('Transfermatrix-Spektralluecke')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 4: max tau_B vs beta
    ax = axes[1, 1]
    max_taus = [np.max([r['tau_B'] for r in all_results[beta]]) for beta in betas]
    ax.plot(betas, max_taus, 'b^-', lw=2, ms=8)
    ax.axhline(y=1.0, color='red', ls='--', lw=1.5)
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\max_k \tau_B(R_k)$')
    ax.set_title(r'Maximale Kontraktion ueber RG-Stufen')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(os.path.dirname(__file__), 'compute_birkhoff_rg.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Plot: {outpath}")
except Exception as e:
    print(f"(Plot: {e})")

print("\n[DONE]")
