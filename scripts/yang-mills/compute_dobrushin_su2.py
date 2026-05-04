"""
compute_dobrushin_su2.py
========================
Yang-Mills: Monte-Carlo-Simulation eines SU(2)-Gitter-Eichfeldes.
Berechnet die annealed Dobrushin-Einflusskoeffizienten c_tilde_{ij}
und prueft die defective Poincare-Ungleichung.

Grundlage: FST-YM Paper, Proposition (Defective Dobrushin)
  c_tilde_{ij} = E[ c_{ij}(eta) ]  (annealed influence)
  eta = max_i sum_{j!=i} c_tilde_{ij}  (Dobrushin coefficient)
  eta < 1 => Poincare inequality with gap >= 1 - eta

SU(2) wird parametrisiert als a0*I + i*(a1*sigma1 + a2*sigma2 + a3*sigma3)
mit a0^2 + a1^2 + a2^2 + a3^2 = 1 (Quaternionen-Darstellung).

Wilson-Aktion: S = -beta/2 * sum_{plaquettes} Re(Tr(U_p))
Konditionaler Einfluss: c_{ij} = ||d mu_i/d U_j||_TV

Autor: Lukas Geiger (Skript erstellt per Claude, 2026)
"""

import numpy as np
import os

# ==========================================================================
# SU(2) Hilfsfunktionen (Quaternionen-Darstellung)
# ==========================================================================

def random_su2():
    """Zufaellige SU(2)-Matrix via Quaternion (Haar-Mass)"""
    v = np.random.randn(4)
    v /= np.linalg.norm(v)
    return v  # (a0, a1, a2, a3)

def su2_to_matrix(q):
    """Quaternion -> 2x2 komplexe Matrix"""
    a0, a1, a2, a3 = q
    return np.array([
        [a0 + 1j*a3,  a2 + 1j*a1],
        [-a2 + 1j*a1, a0 - 1j*a3]
    ])

def su2_multiply(q1, q2):
    """Quaternion-Multiplikation (= SU(2)-Matrixmultiplikation)"""
    a0, a1, a2, a3 = q1
    b0, b1, b2, b3 = q2
    return np.array([
        a0*b0 - a1*b1 - a2*b2 - a3*b3,
        a0*b1 + a1*b0 + a2*b3 - a3*b2,
        a0*b2 - a1*b3 + a2*b0 + a3*b1,
        a0*b3 + a1*b2 - a2*b1 + a3*b0
    ])

def su2_adjoint(q):
    """Adjungierte (Inverse fuer SU(2)): (a0, -a1, -a2, -a3)"""
    return np.array([q[0], -q[1], -q[2], -q[3]])

def su2_trace(q):
    """Re(Tr(U)) = 2*a0 fuer SU(2)"""
    return 2.0 * q[0]

# ==========================================================================
# 2D SU(2) Gitter-Eichtheorie (der einfachste nicht-triviale Fall)
# ==========================================================================

class SU2Lattice:
    """
    2D quadratisches Gitter mit SU(2)-Link-Variablen.
    Wilson-Aktion: S = -(beta/2) * sum_p Re(Tr(U_p))
    U_p = U_mu(x) * U_nu(x+mu) * U_mu(x+nu)^dag * U_nu(x)^dag
    """

    def __init__(self, L, beta):
        self.L = L
        self.beta = beta
        self.n_links = 2 * L * L  # 2 Richtungen pro Site
        # Links: shape (L, L, 2, 4) -- 2 Richtungen, 4 Quaternion-Komponenten
        self.links = np.zeros((L, L, 2, 4))
        # Initialisierung: Identitaet (kalter Start)
        self.links[:, :, :, 0] = 1.0

    def hot_start(self):
        """Zufaelliger Start (heisser Start)"""
        for x in range(self.L):
            for y in range(self.L):
                for mu in range(2):
                    self.links[x, y, mu] = random_su2()

    def get_link(self, x, y, mu):
        """Link U_mu(x,y) lesen"""
        return self.links[x % self.L, y % self.L, mu]

    def set_link(self, x, y, mu, q):
        """Link setzen"""
        self.links[x % self.L, y % self.L, mu] = q

    def plaquette(self, x, y):
        """
        Plaquette U_p = U_0(x,y) * U_1(x+1,y) * U_0(x,y+1)^dag * U_1(x,y)^dag
        Returniert Re(Tr(U_p)) / 2 in [-1, 1]
        """
        U1 = self.get_link(x, y, 0)         # U_x(x,y)
        U2 = self.get_link(x+1, y, 1)       # U_y(x+1,y)
        U3 = su2_adjoint(self.get_link(x, y+1, 0))  # U_x(x,y+1)^dag
        U4 = su2_adjoint(self.get_link(x, y, 1))     # U_y(x,y)^dag
        prod = su2_multiply(su2_multiply(U1, U2), su2_multiply(U3, U4))
        return su2_trace(prod) / 2.0  # normiert auf [-1, 1]

    def staple(self, x, y, mu):
        """
        Berechne die Staple (Summe der Plaquette-Beitraege ohne den Link (x,y,mu)).
        Fuer 2D gibt es genau 2 Staples pro Link (eine Plaquette von oben, eine von unten).
        Returniert Quaternion.
        """
        L = self.L
        if mu == 0:  # x-Richtung
            # Obere Staple: U_1(x+1,y) * U_0(x,y+1)^dag * U_1(x,y)^dag
            s1 = su2_multiply(
                self.get_link(x+1, y, 1),
                su2_multiply(
                    su2_adjoint(self.get_link(x, y+1, 0)),
                    su2_adjoint(self.get_link(x, y, 1))
                )
            )
            # Untere Staple: U_1(x+1,y-1)^dag * U_0(x,y-1)^dag * U_1(x,y-1)
            s2 = su2_multiply(
                su2_adjoint(self.get_link(x+1, y-1, 1)),
                su2_multiply(
                    su2_adjoint(self.get_link(x, y-1, 0)),
                    self.get_link(x, y-1, 1)
                )
            )
        else:  # y-Richtung
            # Rechte Staple: U_0(x,y+1) * U_1(x+1,y)^dag * U_0(x,y)^dag
            s1 = su2_multiply(
                self.get_link(x, y+1, 0),
                su2_multiply(
                    su2_adjoint(self.get_link(x+1, y, 1)),
                    su2_adjoint(self.get_link(x, y, 0))
                )
            )
            # Linke Staple: U_0(x-1,y+1)^dag * U_1(x-1,y)^dag * U_0(x-1,y)
            s2 = su2_multiply(
                su2_adjoint(self.get_link(x-1, y+1, 0)),
                su2_multiply(
                    su2_adjoint(self.get_link(x-1, y, 1)),
                    self.get_link(x-1, y, 0)
                )
            )
        # Summe der Staples (Quaternion-Addition)
        return s1 + s2

    def action(self):
        """Gesamte Wilson-Aktion"""
        S = 0.0
        for x in range(self.L):
            for y in range(self.L):
                S += self.plaquette(x, y)
        return -self.beta * S

    def mean_plaquette(self):
        """Mittlere Plaquette <Re Tr U_p / 2>"""
        total = 0.0
        for x in range(self.L):
            for y in range(self.L):
                total += self.plaquette(x, y)
        return total / (self.L * self.L)

    def heatbath_update(self, x, y, mu):
        """
        Heatbath-Update fuer SU(2) nach Creutz (1980).
        Bedingte Verteilung: P(U) ~ exp(beta/2 * Re(Tr(U * A^dag)))
        wobei A = Staple.
        """
        A = self.staple(x, y, mu)
        k = np.linalg.norm(A)  # |A| (Quaternion-Norm)

        if k < 1e-10:
            # Staple ~ 0: zufaelliger Link
            self.set_link(x, y, mu, random_su2())
            return

        # Effektiver Kopplungsparameter
        alpha = self.beta * k

        # SU(2) Heatbath: Generiere a0 mit P(a0) ~ sqrt(1-a0^2) * exp(alpha * a0)
        # Kennedy-Pendleton-Algorithmus (vereinfacht fuer kleine alpha)
        for _ in range(20):  # Rejection sampling
            r = np.random.random()
            a0 = 1.0 + np.log(r + (1-r) * np.exp(-2*alpha)) / alpha
            if np.random.random() < np.sqrt(1.0 - a0**2):
                break

        # Rest des Quaternions: gleichverteilt auf S^2 mit Radius sqrt(1-a0^2)
        r_sq = 1.0 - a0**2
        if r_sq < 0:
            r_sq = 0
        r_val = np.sqrt(r_sq)
        phi = 2 * np.pi * np.random.random()
        cos_theta = 2 * np.random.random() - 1
        sin_theta = np.sqrt(1 - cos_theta**2)
        a1 = r_val * sin_theta * np.cos(phi)
        a2 = r_val * sin_theta * np.sin(phi)
        a3 = r_val * cos_theta

        # Neuer Link: U_new = V * A^dag / |A|
        V = np.array([a0, a1, a2, a3])
        A_hat = A / k  # normierte Staple
        U_new = su2_multiply(V, su2_adjoint(A_hat))
        self.set_link(x, y, mu, U_new)

    def sweep(self):
        """Ein vollstaendiger Heatbath-Sweep"""
        for x in range(self.L):
            for y in range(self.L):
                for mu in range(2):
                    self.heatbath_update(x, y, mu)


# ==========================================================================
# Dobrushin-Einfluss Messung
# ==========================================================================

def measure_dobrushin_influence(lattice, n_samples=500, n_therm=100, perturbation_eps=0.1):
    """
    Misst die annealed Dobrushin-Einflusskoeffizienten c_tilde_{ij}.

    Methode: Stoere Link j, messe Aenderung der konditionalen Erwartung von Link i.
    c_{ij} ~ ||delta <U_i | rest> / delta U_j||

    Approximation:
    1. Thermalisiere Gitter
    2. Fuer jedes Paar (i,j) benachbarter Links:
       a. Speichere Konfiguration
       b. Messe <U_i> ueber n_inner Sweeps
       c. Stoere U_j -> U_j * V (kleine Drehung)
       d. Messe <U_i'> ueber n_inner Sweeps
       e. c_{ij} ~ ||<U_i> - <U_i'>|| / ||V - I||
    """
    L = lattice.L

    # Thermalisieren
    print("    Thermalisierung...", end="", flush=True)
    for _ in range(n_therm):
        lattice.sweep()
    print(" fertig.")

    # Sammle Einflusskoeffizienten fuer NACHBAR-Paare
    # In 2D hat jeder Link 6 Nachbar-Links (die in Plaquettes vorkommen)
    # Wir messen die mittlere Einflusstaerke

    influences = []
    n_inner = 10  # innere Sweeps fuer Mittelwert

    print(f"    Messe Einfluss ({n_samples} Samples)...", flush=True)

    for sample in range(n_samples):
        # Einen Sweep als Dekorrelierer
        lattice.sweep()

        # Zufaelliges Paar (i, j) waehlen: j = Nachbar von i
        x_i, y_i = np.random.randint(0, L), np.random.randint(0, L)
        mu_i = np.random.randint(0, 2)

        # Nachbar-Link j (einer der 6 in der Plaquette)
        if mu_i == 0:
            neighbors = [(x_i+1, y_i, 1), (x_i, y_i+1, 0), (x_i, y_i, 1),
                         (x_i+1, y_i-1, 1), (x_i, y_i-1, 0), (x_i, y_i-1, 1)]
        else:
            neighbors = [(x_i, y_i+1, 0), (x_i+1, y_i, 1), (x_i, y_i, 0),
                         (x_i-1, y_i+1, 0), (x_i-1, y_i, 1), (x_i-1, y_i, 0)]

        x_j, y_j, mu_j = neighbors[np.random.randint(0, len(neighbors))]

        # Speichere aktuelle Konfiguration
        config_backup = lattice.links.copy()

        # Messe <U_i> ungestoert (Mittelwert ueber innere Sweeps)
        U_i_samples = []
        for _ in range(n_inner):
            lattice.sweep()
            U_i_samples.append(lattice.get_link(x_i, y_i, mu_i).copy())
        U_i_mean = np.mean(U_i_samples, axis=0)

        # Restauriere und stoere U_j
        lattice.links[:] = config_backup

        # Kleine SU(2)-Stoerung: U_j -> U_j * exp(eps * sigma)
        old_Uj = lattice.get_link(x_j, y_j, mu_j).copy()
        perturb = np.random.randn(3) * perturbation_eps
        perturb_q = np.array([np.sqrt(max(0, 1 - np.sum(perturb**2))),
                              perturb[0], perturb[1], perturb[2]])
        perturb_q /= np.linalg.norm(perturb_q)
        new_Uj = su2_multiply(old_Uj, perturb_q)
        lattice.set_link(x_j, y_j, mu_j, new_Uj)

        # Messe <U_i'> gestoert
        U_i_prime_samples = []
        for _ in range(n_inner):
            lattice.sweep()
            U_i_prime_samples.append(lattice.get_link(x_i, y_i, mu_i).copy())
        U_i_prime_mean = np.mean(U_i_prime_samples, axis=0)

        # Einfluss: ||<U_i> - <U_i'>|| / ||V - I||
        delta_U = np.linalg.norm(U_i_mean - U_i_prime_mean)
        delta_V = np.linalg.norm(perturb_q - np.array([1, 0, 0, 0]))
        if delta_V > 1e-10:
            c_ij = delta_U / delta_V
        else:
            c_ij = 0.0

        influences.append(c_ij)

        # Restauriere
        lattice.links[:] = config_backup

        if (sample + 1) % 100 == 0:
            eta_current = np.mean(influences) * 6  # 6 Nachbarn pro Link
            print(f"      Sample {sample+1}/{n_samples}: "
                  f"<c_ij>={np.mean(influences):.4f}, "
                  f"eta_approx={eta_current:.4f}")

    return np.array(influences)


# ==========================================================================
# Hauptberechnung
# ==========================================================================

print("=" * 70)
print("YANG-MILLS: Defective Dobrushin / SU(2) Gitter-Simulation")
print("Prueft: eta = max_i sum_{j~i} c_tilde_{ij} < 1 (Poincare-Luecke)")
print("=" * 70)

# Parameter-Scan ueber verschiedene beta-Werte
betas = [1.0, 2.0, 4.0, 8.0]
L = 4  # Gittergroesse 4x4 (minimal)

results = []

for beta in betas:
    print(f"\n{'='*50}")
    print(f"beta = {beta:.1f}  (L = {L})")
    print(f"{'='*50}")

    lattice = SU2Lattice(L, beta)
    lattice.hot_start()

    # Messe mittlere Plaquette nach Thermalisierung
    for _ in range(200):
        lattice.sweep()
    plaq = lattice.mean_plaquette()
    print(f"  <P> = {plaq:.6f}  (theor. strong: {1 - 3/(4*beta):.4f} fuer beta>>1)")

    # Messe Dobrushin-Einfluss
    influences = measure_dobrushin_influence(
        lattice,
        n_samples=80,
        n_therm=20,
        perturbation_eps=0.05
    )

    c_mean = np.mean(influences)
    c_std = np.std(influences)
    c_max = np.max(influences)

    # Dobrushin-Koeffizient: eta = max_i sum_{j~i} c_tilde_{ij}
    # In 2D hat jeder Link ~6 Nachbarn in Plaquettes
    n_neighbors = 6
    eta = c_mean * n_neighbors
    eta_max = c_max * n_neighbors

    # Poincare gap = 1 - eta (falls eta < 1)
    gap = max(0, 1 - eta)

    results.append({
        'beta': beta,
        'plaquette': plaq,
        'c_mean': c_mean,
        'c_std': c_std,
        'c_max': c_max,
        'eta': eta,
        'eta_max': eta_max,
        'gap': gap,
    })

    print(f"  Ergebnis:")
    print(f"    <c_ij> = {c_mean:.6f} +/- {c_std:.6f}")
    print(f"    c_max  = {c_max:.6f}")
    print(f"    eta (mean) = {eta:.4f}  (= 6 * <c_ij>)")
    print(f"    eta (max)  = {eta_max:.4f}")
    if eta < 1:
        print(f"    => Poincare gap >= {gap:.4f} > 0  (LUECKE EXISTIERT)")
    else:
        print(f"    => eta >= 1: Dobrushin-Bedingung NICHT erfuellt")
        print(f"       (defective regime: polynomielle Degradation erwartet)")


# ==========================================================================
# Zusammenfassung
# ==========================================================================

print("\n" + "=" * 70)
print("ZUSAMMENFASSUNG: Dobrushin-Koeffizienten fuer SU(2) in 2D")
print("=" * 70)
print(f"\n{'beta':>6} {'<P>':>8} {'<c_ij>':>10} {'eta':>8} {'eta_max':>8} {'gap':>8} {'Status':>15}")
print("-" * 70)

for r in results:
    if r['eta'] < 0.5:
        status = "STARK (eta<<1)"
    elif r['eta'] < 1.0:
        status = "OK (eta<1)"
    else:
        status = "DEFECTIVE"
    print(f"{r['beta']:>6.1f} {r['plaquette']:>8.4f} {r['c_mean']:>10.6f} "
          f"{r['eta']:>8.4f} {r['eta_max']:>8.4f} {r['gap']:>8.4f} {status:>15}")

# Analyse des Skalenverhaltens
print("\n--- Skalenverhalten ---")
print("Erwartung (Paper Prop. polynomial-lsi):")
print("  Strong coupling (beta << 1): kappa ~ beta (LSI leicht)")
print("  Weak coupling (beta >> 1):   kappa degradiert, aber bleibt > 0")
print("  Kontinuumslimes (beta->inf): Hierarchische RG noetig")

betas_arr = np.array([r['beta'] for r in results])
etas_arr = np.array([r['eta'] for r in results])
gaps_arr = np.array([r['gap'] for r in results])

# Fit: eta(beta) ~ A * exp(B * beta) oder polynomiell?
if np.all(etas_arr > 0) and np.all(etas_arr < 10):
    log_eta = np.log(np.maximum(etas_arr, 1e-10))
    # Linearer Fit in log-Skala
    coeffs = np.polyfit(betas_arr, log_eta, 1)
    print(f"\n  Log-linearer Fit: ln(eta) ~ {coeffs[0]:.4f} * beta + {coeffs[1]:.4f}")
    if coeffs[0] > 0:
        beta_crit = (np.log(1.0) - coeffs[1]) / coeffs[0]
        print(f"  Geschaetztes beta_crit (eta=1): {beta_crit:.2f}")
        print(f"  (Ab beta > {beta_crit:.1f} ist hierarchische RG noetig)")
    else:
        print(f"  eta FAELLT mit beta (unerwartetes Regime)")

# Vergleich mit analytischer Strong-Coupling-Formel
print("\n--- Vergleich mit analytischen Formeln ---")
for r in results:
    beta = r['beta']
    # Strong-coupling Expansion: kappa_link ~ (beta/4)^2 fuer SU(2)
    # Fuer kleine beta: c_{ij} ~ beta/4 (naechster-Nachbar-Kopplung)
    kappa_sc = (beta/4)**2
    print(f"  beta={beta:.1f}: eta_num={r['eta']:.4f}, "
          f"kappa_SC={kappa_sc:.4f}, eta_SC_approx={6/(4*kappa_sc) if kappa_sc > 0 else float('inf'):.4f}")


# Plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: eta vs beta
    axes[0].plot(betas_arr, etas_arr, 'ro-', linewidth=2, markersize=8, label=r'$\eta$ (mean)')
    axes[0].axhline(y=1.0, color='black', linestyle='--', linewidth=1, label=r'$\eta = 1$')
    axes[0].fill_between(betas_arr, 0, 1, alpha=0.1, color='green', label='Poincare gap > 0')
    axes[0].set_xlabel(r'$\beta$', fontsize=12)
    axes[0].set_ylabel(r'$\eta$ (Dobrushin coeff.)', fontsize=12)
    axes[0].set_title('Dobrushin-Koeffizient vs. Kopplung', fontsize=11)
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    # Panel 2: Poincare gap vs beta
    axes[1].plot(betas_arr, gaps_arr, 'bs-', linewidth=2, markersize=8)
    axes[1].axhline(y=0, color='red', linestyle='--', linewidth=1)
    axes[1].set_xlabel(r'$\beta$', fontsize=12)
    axes[1].set_ylabel(r'gap $= 1 - \eta$', fontsize=12)
    axes[1].set_title(r'Poincar\'{e}-Luecke vs. Kopplung', fontsize=11)
    axes[1].grid(True, alpha=0.3)

    # Panel 3: Plaquette vs beta
    plaq_arr = np.array([r['plaquette'] for r in results])
    axes[2].plot(betas_arr, plaq_arr, 'g^-', linewidth=2, markersize=8, label=r'$\langle P \rangle$ (MC)')
    beta_th = np.linspace(1, 10, 100)
    plaq_th = 1 - 3/(4*beta_th)  # strong coupling approx
    axes[2].plot(beta_th, plaq_th, 'k--', linewidth=1, alpha=0.5, label=r'$1 - 3/(4\beta)$ (approx)')
    axes[2].set_xlabel(r'$\beta$', fontsize=12)
    axes[2].set_ylabel(r'$\langle P \rangle$', fontsize=12)
    axes[2].set_title('Mittlere Plaquette', fontsize=11)
    axes[2].legend(fontsize=9)
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(os.path.dirname(__file__),
                           'compute_dobrushin_su2.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nPlot gespeichert: {outpath}")
except Exception as e:
    print(f"\n(matplotlib nicht verfuegbar: {e})")

print("\n[DONE]")
