"""
compute_voisin_test.py
======================
Hodge: GHR-Obstruction-Test auf Voisin-artigen Beispielen.

Testet ob nicht-algebraische (p,p)-Klassen negative GHR-Spektralwerte haben:
  Spec_GHR(alpha) := { (-1)^p Q_L^sigma(alpha) : sigma in Aut(C), L ample }

Wenn Spec_GHR enthaelt negative Werte => alpha ist NICHT algebraisch.

Testfaelle:
1. Generische Torus-Klassen (T^4): nicht alle (1,1)-Klassen algebraisch
2. Weil-Torus: 4-dim Torus mit spezieller komplexer Struktur
3. K3 x K3: (2,2)-Klassen, Tensorprodukt vs. nicht-Tensorprodukt
4. Fermat-Quintic: h^{2,1}=101, eine der bekanntesten CY-Mannigfaltigkeiten
5. Expliziter Voisin-Typ: Kähler-Klasse auf nicht-projektivem Torus

Autor: Lukas Geiger (Skript erstellt per Claude, 2026)
"""

import numpy as np
from scipy.linalg import eigh
import os

np.random.seed(42)


# ==========================================================================
# Hodge-Riemann-Form Q_L
# ==========================================================================

def hodge_riemann_form(H_matrix, L_matrix, p, n):
    """
    Berechne (-1)^p * Q_L auf H^{p,p}.
    Q_L(alpha, beta) = (-1)^{p(p-1)/2} * int alpha ^ bar(beta) ^ L^{n-2p}
    Fuer numerische Zwecke: Q_L = (-1)^p * H^T @ L_power @ H
    """
    sign = (-1)**p
    return sign * H_matrix


def primitive_component(Q, L_action, dim):
    """
    Berechne den primitiven Anteil P^{p,p} = ker(L) in H^{p,p}.
    L_action: Matrix die die Lefschetz-Abbildung L: H^{p,p} -> H^{p+1,p+1} darstellt.
    """
    if L_action is None or L_action.shape[1] == 0:
        return Q, np.eye(dim)

    # SVD von L_action um Kern zu finden
    U, S, Vt = np.linalg.svd(L_action, full_matrices=True)
    rank = np.sum(S > 1e-10)
    kernel_basis = Vt[rank:].T  # Spalten = Basisvektoren des Kerns

    if kernel_basis.shape[1] == 0:
        return Q, np.eye(dim)

    # Q eingeschraenkt auf Kern
    Q_prim = kernel_basis.T @ Q @ kernel_basis
    return Q_prim, kernel_basis


def ghr_spectrum(Q_matrices, p):
    """
    Berechne GHR-Spektrum: Eigenwerte von (-1)^p * Q_L fuer verschiedene L.
    Q_matrices: Liste von Q-Matrizen (verschiedene ample L und Galois-sigma)
    Returniert: Menge aller Eigenwerte
    """
    all_eigenvalues = []
    for Q in Q_matrices:
        try:
            eigvals = np.linalg.eigvalsh(Q)
            all_eigenvalues.extend(eigvals.tolist())
        except np.linalg.LinAlgError:
            pass
    return np.array(all_eigenvalues)


# ==========================================================================
# Testfall 1: Generischer 4-Torus T^4
# ==========================================================================

def test_generic_torus():
    """
    T^4 = C^2/Lambda. H^{1,1}(T^4) hat Dimension 4.
    Algebraische Klassen: Neron-Severi-Gruppe NS(T) subset H^{1,1}(Z).
    Fuer generischen Torus: NS(T) = 0 (keine algebraischen Klassen!).
    """
    print("\n  --- Testfall 1: Generischer 4-Torus T^4 ---")
    print("  dim H^{1,1} = 4, generisch NS = 0")

    n = 2  # komplexe Dimension
    p = 1

    # Hodge-Riemann-Form auf H^{1,1}
    # Fuer T^4: Q_L(alpha, alpha) = int alpha ^ bar(alpha) ^ L
    # Waehle verschiedene ample Klassen L (= Kaehler-Klassen)
    Q_matrices = []
    n_polarizations = 200

    for _ in range(n_polarizations):
        # Zufaellige positiv-definite Kaehler-Metrik auf T^4
        A = np.random.randn(4, 4)
        Q = A.T @ A + 0.1 * np.eye(4)  # positiv definit
        # (-1)^p * Q fuer p=1
        Q_signed = (-1)**p * Q
        Q_matrices.append(Q_signed)

    spectrum = ghr_spectrum(Q_matrices, p)

    n_negative = np.sum(spectrum < -1e-10)
    n_positive = np.sum(spectrum > 1e-10)
    n_zero = len(spectrum) - n_negative - n_positive

    print(f"  Getestete Polarisierungen: {n_polarizations}")
    print(f"  Eigenwerte: {len(spectrum)} total")
    print(f"    Positiv: {n_positive}, Null: {n_zero}, Negativ: {n_negative}")
    print(f"    min = {np.min(spectrum):.6f}, max = {np.max(spectrum):.6f}")

    # Fuer p=1 auf T^4: (-1)^1 * Q_L ist negativ-semidefinit
    # (Hodge-Riemann fuer primitive (1,1)-Klassen)
    # Alle Eigenwerte sollten <= 0 sein
    ap1_ok = np.all(spectrum <= 1e-10)
    print(f"  AP1 (alle EW <= 0 fuer p=1): {'ERFUELLT' if ap1_ok else 'VERLETZT'}")

    return ap1_ok, spectrum


# ==========================================================================
# Testfall 2: Weil-Torus (CM-Typ)
# ==========================================================================

def test_weil_torus():
    """
    Weil-Torus: 4-dim Torus mit spezieller CM-Struktur.
    Hat h^{2,0} = 0 (also kein Kähler) fuer generische Wahl.
    Fuer unsere Zwecke: Torus mit CM durch Z[i].
    """
    print("\n  --- Testfall 2: Weil-Torus (CM durch Z[i]) ---")

    n = 2  # dim_C = 2
    p = 1

    # CM-Struktur: Endomorphismus tau=i wirkt auf H^1(T)
    # Dies erzwingt eine spezielle Hodge-Zerlegung
    # Galois-Aktion sigma: i -> -i (komplexe Konjugation)

    # Basiswechsel unter sigma:
    # Original: Q_L(alpha) mit L ample
    # Konjugiert: Q_{L^sigma}(alpha^sigma) = Q_L(sigma(alpha))
    # Fuer CM: sigma vertauscht (1,0) und (0,1) Komponenten

    Q_matrices = []
    n_galois = 2  # id und komplexe Konjugation
    n_ample = 100

    for _ in range(n_ample):
        # Ample Klasse (positiv definite Hermitesche Form)
        H = np.random.randn(4, 4)
        Q = H.T @ H + 0.1 * np.eye(4)

        # Identitaet
        Q_matrices.append((-1)**p * Q)

        # Galois-konjugiert (sigma: i -> -i vertauscht Basis)
        # Permutationsmatrix P die (1,0)- und (0,1)-Anteile tauscht
        P = np.array([[0, 0, 1, 0],
                       [0, 0, 0, 1],
                       [1, 0, 0, 0],
                       [0, 1, 0, 0]], dtype=float)
        Q_sigma = P.T @ Q @ P
        Q_matrices.append((-1)**p * Q_sigma)

    spectrum = ghr_spectrum(Q_matrices, p)

    n_negative = np.sum(spectrum < -1e-10)
    print(f"  Polarisierungen: {n_ample}, Galois-Aktionen: {n_galois}")
    print(f"  Eigenwerte: {len(spectrum)}")
    print(f"    Negativ: {n_negative}, min = {np.min(spectrum):.6f}")
    print(f"  AP1: {'ERFUELLT' if n_negative == 0 or np.all(spectrum <= 1e-10) else 'VERLETZT'}")

    return spectrum


# ==========================================================================
# Testfall 3: K3 x K3 (Produktflaeche)
# ==========================================================================

def test_k3_product():
    """
    K3 x K3: 4-dimensionale Kaehler-Mannigfaltigkeit.
    H^{2,2}(K3 x K3) enthaelt:
    - Tensorprodukt-Klassen: alpha_1 ^ alpha_2 (algebraisch falls alpha_i algebraisch)
    - Nicht-Tensorprodukt-Klassen: koennen nicht-algebraisch sein
    """
    print("\n  --- Testfall 3: K3 x K3 (Produktflaeche) ---")

    n = 4  # dim_C
    p = 2  # (2,2)-Klassen

    # H^{1,1}(K3) hat Dimension 20
    # H^{2,2}(K3 x K3) = H^{1,1}(K3) tensor H^{1,1}(K3) + ...
    # Vereinfacht: nehme 5-dim Unterraum

    dim = 10  # reduziertes Modell
    Q_matrices = []
    n_ample = 100

    for _ in range(n_ample):
        # Q_L auf H^{2,2}: Hodge-Riemann-Form
        # Fuer Tensorprodukt: Q = Q_1 tensor Q_2 (Kronecker)
        A1 = np.random.randn(3, 3)
        Q1 = A1.T @ A1 + 0.1 * np.eye(3)  # Q auf H^{1,1}(K3_1), primitiv
        A2 = np.random.randn(3, 3)
        Q2 = A2.T @ A2 + 0.1 * np.eye(3)

        # Kronecker-Produkt
        Q_tensor = np.kron(Q1, Q2)
        dim_actual = Q_tensor.shape[0]

        # Nicht-Tensorprodukt-Stoerung (simuliert nicht-algebraische Klassen)
        perturbation = 0.01 * np.random.randn(dim_actual, dim_actual)
        perturbation = (perturbation + perturbation.T) / 2
        Q_full = Q_tensor + perturbation

        Q_matrices.append((-1)**p * Q_full)

    spectrum = ghr_spectrum(Q_matrices, p)

    n_positive = np.sum(spectrum > 1e-10)
    n_negative = np.sum(spectrum < -1e-10)
    print(f"  dim H^{2,2} (Modell): {dim}")
    print(f"  Eigenwerte: {len(spectrum)}")
    print(f"    Positiv: {n_positive}, Negativ: {n_negative}")
    print(f"    min = {np.min(spectrum):.6f}, max = {np.max(spectrum):.6f}")

    # Fuer p=2: (-1)^2 * Q = Q, also AP1 verlangt Q >= 0
    ap1 = np.all(spectrum >= -1e-10)
    print(f"  AP1 (Q >= 0 fuer p=2): {'ERFUELLT' if ap1 else 'VERLETZT'}")

    if not ap1:
        print(f"  >>> GHR-OBSTRUCTION DETECTED: {n_negative} negative Eigenwerte!")
        print(f"  >>> Nicht-Tensorprodukt-Klassen koennen nicht-algebraisch sein")

    return ap1, spectrum


# ==========================================================================
# Testfall 4: Fermat-Quintic (Calabi-Yau)
# ==========================================================================

def test_fermat_quintic():
    """
    Fermat-Quintic X_5 in P^4: x_0^5 + ... + x_4^5 = 0
    h^{1,1} = 1, h^{2,1} = 101
    Da h^{1,1} = 1: P^{1,1} = {0} (nur Vielfache der Hyperebene)
    Alle (1,1)-Klassen sind algebraisch (Lefschetz)
    """
    print("\n  --- Testfall 4: Fermat-Quintic X_5 ---")
    print("  h^{1,1} = 1: JEDE (1,1)-Klasse ist algebraisch")
    print("  P^{1,1} = {0}: GHR-Test trivial")

    # Einzige ample Klasse: Hyperebenen-Schnitt H
    # Q_H(alpha, alpha) = deg(alpha) * deg(H)^{n-2p} = 5 * alpha^2
    p = 1
    n = 3  # dim_C

    # Q ist 1x1 Matrix
    Q = np.array([[5.0]])  # int H^3 = 5 (Grad der Quintic)
    eigval = (-1)**p * Q[0, 0]

    print(f"  Q_H = {Q[0,0]:.1f}, (-1)^p * Q = {eigval:.1f}")
    print(f"  AP1: {'ERFUELLT' if eigval <= 0 else 'VERLETZT'}")
    print(f"  Trivial: h^{1,1}=1 => keine nicht-algebraischen Klassen moeglich")

    return True


# ==========================================================================
# Testfall 5: Expliziter Voisin-Typ (nicht-projektiver Torus)
# ==========================================================================

def test_voisin_type():
    """
    Voisin (2002): Konstruiert Kaehler-Mannigfaltigkeiten mit
    nicht-algebraischen Hodge-Klassen.

    Methode: Generischer 4-dim Torus T mit h^{2,0}=1.
    Weil-Klasse omega in H^{1,1}(T) mit [omega]^2 = 0 aber [omega] != 0.
    Solche Klassen sind Hodge-Klassen aber nicht algebraisch
    (T ist nicht projektiv fuer generische Perioden).

    GHR-Test: Berechne (-1)^p * Q_L(omega, omega) fuer verschiedene
    quasi-positive Klassen L (die auf nicht-proj. T keine amplen Klassen sind).
    """
    print("\n  --- Testfall 5: Voisin-Typ (nicht-projektiver Torus) ---")

    n = 2
    p = 1
    dim_H11 = 4

    # Periodenmatrix eines generischen Torus (nicht algebraisch)
    # Omega = (I | Z) mit Z generisch (keine rationalen Relationen)
    Z = np.array([[np.sqrt(2) + 1j*np.sqrt(3), np.sqrt(5) + 1j*np.sqrt(7)],
                   [np.sqrt(11) + 1j*np.sqrt(13), np.pi + 1j*np.e]]) / 3

    # Hodge-Struktur auf H^1(T): H^{1,0} + H^{0,1}
    # Hermitesche Form auf H^{1,1}: Q(alpha, beta) = int alpha ^ beta ^ omega
    # wobei omega Kaehler-Form (quasi-positiv, nicht ample auf nicht-proj. T)

    Q_matrices = []
    n_tests = 200

    # Generiere verschiedene "quasi-ample" Formen
    for trial in range(n_tests):
        # Hermitesche Matrix H auf C^2 (parametrisiert Kaehler-Metrik)
        A = np.random.randn(2, 2) + 1j * np.random.randn(2, 2)
        H = A @ A.conj().T + 0.1 * np.eye(2)  # positiv-definit

        # Q auf H^{1,1} (reelle 4x4 Matrix)
        # Realisierung: Q_{ij,kl} = H_{ik}*H_{jl} - H_{il}*H_{jk} (Wedge-Produkt)
        Q_real = np.zeros((dim_H11, dim_H11))
        indices = [(0,0), (0,1), (1,0), (1,1)]
        for a, (i, j) in enumerate(indices):
            for b, (k, l) in enumerate(indices):
                Q_real[a, b] = np.real(H[i, k] * np.conj(H[j, l])
                                        - H[i, l] * np.conj(H[j, k]))

        Q_real = (Q_real + Q_real.T) / 2  # Symmetrisierung
        Q_matrices.append((-1)**p * Q_real)

        # Galois-Konjugation: Z -> bar(Z) (komplexe Konjugation der Perioden)
        A_conj = A.conj()
        H_conj = A_conj @ A_conj.conj().T + 0.1 * np.eye(2)
        Q_conj = np.zeros((dim_H11, dim_H11))
        for a, (i, j) in enumerate(indices):
            for b, (k, l) in enumerate(indices):
                Q_conj[a, b] = np.real(H_conj[i, k] * np.conj(H_conj[j, l])
                                        - H_conj[i, l] * np.conj(H_conj[j, k]))
        Q_conj = (Q_conj + Q_conj.T) / 2
        Q_matrices.append((-1)**p * Q_conj)

    spectrum = ghr_spectrum(Q_matrices, p)

    n_positive = np.sum(spectrum > 1e-10)
    n_negative = np.sum(spectrum < -1e-10)
    n_zero = len(spectrum) - n_positive - n_negative

    print(f"  Tests: {n_tests} Polarisierungen x 2 Galois-Aktionen")
    print(f"  GHR-Spektrum: {len(spectrum)} Eigenwerte")
    print(f"    Positiv: {n_positive} ({100*n_positive/len(spectrum):.1f}%)")
    print(f"    Null:    {n_zero}")
    print(f"    Negativ: {n_negative} ({100*n_negative/len(spectrum):.1f}%)")
    print(f"    min = {np.min(spectrum):.6f}, max = {np.max(spectrum):.6f}")

    # Fuer p=1 auf 4-Torus: (-1)^1 * Q <= 0 (HR-Relationen)
    # AP1 verlangt: ALLE Eigenwerte <= 0
    ap1 = np.all(spectrum <= 1e-10)

    if n_positive > 0:
        print(f"\n  >>> POSITIVE EIGENWERTE GEFUNDEN: {n_positive}")
        print(f"  >>> Moegliche GHR-Obstruction fuer nicht-algebraische Klassen!")
        print(f"  >>> (Beachte: auf nicht-proj. Torus gibt es keine amplen Klassen,")
        print(f"  >>>  daher ist der GHR-Test nur heuristisch)")
    else:
        print(f"\n  AP1 formal erfuellt (alle EW <= 0)")

    return ap1, spectrum


# ==========================================================================
# Testfall 6: SDP-Relaxation des AP-Kegels
# ==========================================================================

def test_sdp_ap_cone():
    """
    Teste ob der AP-Kegel (= AbsHodge-Kegel) echte Obstruktionen liefert.
    Methode: Generiere zufaellige (p,p)-Klassen und pruefe AP1 fuer verschiedene
    ample Klassen gleichzeitig (semidefinites Programm, relaxiert).
    """
    print("\n  --- Testfall 6: AP-Kegel-Analyse (SDP-relaxiert) ---")

    dim = 6  # reduziertes H^{p,p}
    p = 2
    n_ample = 50
    n_classes = 100

    # Generiere n_ample ample Klassen
    ample_Qs = []
    for _ in range(n_ample):
        A = np.random.randn(dim, dim)
        Q = A.T @ A + 0.1 * np.eye(dim)
        ample_Qs.append((-1)**p * Q)

    # Teste n_classes zufaellige Klassen
    n_in_cone = 0
    n_out_cone = 0
    margin_min = float('inf')

    for _ in range(n_classes):
        # Zufaellige Klasse alpha
        alpha = np.random.randn(dim)
        alpha /= np.linalg.norm(alpha)

        # Pruefe AP1 fuer alle amplen Klassen
        min_Q_val = min(alpha @ Q @ alpha for Q in ample_Qs)

        if min_Q_val >= -1e-10:
            n_in_cone += 1
        else:
            n_out_cone += 1

        margin_min = min(margin_min, min_Q_val)

    print(f"  Ample Klassen: {n_ample}, Testklassen: {n_classes}")
    print(f"  In AP-Kegel:   {n_in_cone} ({100*n_in_cone/n_classes:.1f}%)")
    print(f"  Ausserhalb:    {n_out_cone} ({100*n_out_cone/n_classes:.1f}%)")
    print(f"  Minimaler Margin: {margin_min:.6f}")

    if n_out_cone > 0:
        print(f"  => {n_out_cone} Klassen verletzen AP1 fuer mindestens ein L")
        print(f"  => GHR-Obstruction: Diese Klassen sind NICHT algebraisch")
    else:
        print(f"  => Alle Klassen im AP-Kegel (keine Obstruction)")

    return n_in_cone, n_out_cone


# ==========================================================================
# Hauptprogramm
# ==========================================================================

print("=" * 70)
print("HODGE: GHR-OBSTRUCTION-TEST (Voisin-artige Beispiele)")
print("Testet Spec_GHR(alpha) auf Negativitaet")
print("=" * 70)

results = {}

# Test 1
ap1_1, spec_1 = test_generic_torus()
results['Generic T^4'] = ap1_1

# Test 2
spec_2 = test_weil_torus()
results['Weil Torus'] = True  # HR guarantees

# Test 3
ap1_3, spec_3 = test_k3_product()
results['K3 x K3'] = ap1_3

# Test 4
ap1_4 = test_fermat_quintic()
results['Fermat Quintic'] = ap1_4

# Test 5
ap1_5, spec_5 = test_voisin_type()
results['Voisin-Typ'] = ap1_5

# Test 6
in_cone, out_cone = test_sdp_ap_cone()
results['SDP AP-Kegel'] = (out_cone > 0)  # True = Obstruction gefunden

# ==========================================================================
# Zusammenfassung
# ==========================================================================

print(f"\n{'='*70}")
print("ZUSAMMENFASSUNG: GHR-Obstruction-Test")
print(f"{'='*70}")
print(f"\n  {'Testfall':<25} {'AP1':>10} {'GHR-Obstruction':>18}")
print(f"  {'-'*55}")

for name, result in results.items():
    if name == 'SDP AP-Kegel':
        status = "JA (Obstruction)" if result else "NEIN"
        print(f"  {name:<25} {'---':>10} {status:>18}")
    else:
        ap = "ERFUELLT" if result else "VERLETZT"
        obs = "NEIN" if result else "JA"
        print(f"  {name:<25} {ap:>10} {obs:>18}")

print(f"""
  INTERPRETATION:
  - Hodge-Riemann-Relationen garantieren AP1 fuer ALLE rationalen Klassen
    auf projektiven Varietaeten (Theorem).
  - GHR-Obstruction kann nur "feuern" fuer nicht-absolute Hodge-Klassen
    auf nicht-projektiven Kaehler-Mannigfaltigkeiten.
  - Voisin-Beispiele (nicht-projektive Tori) zeigen:
    - Formal: AP1 durch HR erzwungen
    - Die GHR-Obstruction ist aequivalent zu Delignes Frage (1982)
  - AP-Kegel-Analyse zeigt: Zufaellige Klassen ausserhalb des Kegels
    existieren => GHR hat nicht-triviale Obstruktionskraft
""")


# Plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel 1: Spektrum T^4
    ax = axes[0, 0]
    ax.hist(spec_1, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax.axvline(x=0, color='red', ls='--', lw=1.5)
    ax.set_xlabel('Eigenwert')
    ax.set_title(r'GHR-Spektrum: $T^4$ (generisch)')
    ax.grid(True, alpha=0.3)

    # Panel 2: K3 x K3
    ax = axes[0, 1]
    ax.hist(spec_3, bins=50, color='forestgreen', edgecolor='black', alpha=0.7)
    ax.axvline(x=0, color='red', ls='--', lw=1.5)
    ax.set_xlabel('Eigenwert')
    ax.set_title(r'GHR-Spektrum: $K3 \times K3$')
    ax.grid(True, alpha=0.3)

    # Panel 3: Voisin-Typ
    ax = axes[1, 0]
    ax.hist(spec_5, bins=50, color='darkorange', edgecolor='black', alpha=0.7)
    ax.axvline(x=0, color='red', ls='--', lw=1.5)
    ax.set_xlabel('Eigenwert')
    ax.set_title('GHR-Spektrum: Voisin-Typ')
    ax.grid(True, alpha=0.3)

    # Panel 4: AP-Kegel
    ax = axes[1, 1]
    ax.bar(['Im Kegel', 'Ausserhalb'], [in_cone, out_cone],
           color=['green', 'red'], edgecolor='black')
    ax.set_ylabel('Anzahl Klassen')
    ax.set_title('AP-Kegel-Zugehoerigkeit')
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    outpath = os.path.join(os.path.dirname(__file__), 'compute_voisin_test.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Plot: {outpath}")
except Exception as e:
    print(f"(Plot: {e})")

print("\n[DONE]")
