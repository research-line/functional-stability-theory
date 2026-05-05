"""
chi_defect_norm_pilot.py
========================

Pilot-Validierung fuer die chi-twistete Defekt-Norm aus
DIRICHLET_CCM_TRANSFER.md (Session 16).

Zweck
-----
Empirisch testen, ob die Defekt-Norm
    eps_{lambda,chi} = || QW_{lambda,chi} . Psi_{lambda,chi}
                         - mu_{lambda,chi} Psi_{lambda,chi} . PW_{lambda,chi} ||
chi-abhaengig ist und die zwei Atlas-relevanten Charaktere unterscheidet:
    chi_0 (trivial, Riemann-Baseline)
    chi_4 (erste nicht-triviale Dirichlet, mod 4, chi(-1)=-1)

Dies ist eine Pilot-Version mit sehr kleinem Galerkin-Raum (N=24).
Das Ziel ist NICHT eine scharfe numerische Aussage, sondern strukturelle
Sichtbarmachung: ist die Matrix-Differenz bei chi_4 signifikant anders
als bei chi_0?

Strategie
---------
Arbeite im Mellin-Bild (Satz 3.1 aus MILESTONE_1_CHI_DOMAIN_THEORY.md):
    Psi_{lambda,chi} ~ Multiplikation mit L(s, chi) in Mellin-Koordinaten.
Auf endlich-dim Galerkin-Raum ist alles explizite Matrix-Arithmetik.

Autor: LG (Opus 4.7, Session 16)
Datum: 2026-04-18

Ausfuehrung:
    PYTHONIOENCODING=utf-8 python chi_defect_norm_pilot.py
"""

import numpy as np
from typing import Callable

# -----------------------------------------------------------------
# 0. Konfiguration
# -----------------------------------------------------------------

LAMBDA = np.sqrt(14.0)          # CCM-Standard-Skala
L = np.log(LAMBDA)              # L = log(lambda)
N_MELLIN = 24                   # Galerkin-Dim (klein fuer Pilot)

# Mellin-Stuetzstellen auf der kritischen Linie Re(s) = 1/2.
# Wir samplen s = 1/2 + i*t fuer t in [-L, L] (Paley-Wiener-Bandbreite).
T_MAX = L * 0.8                 # Bandbegrenzung (T < L, vgl. Satz 4.1)


# -----------------------------------------------------------------
# 1. Dirichlet-Charaktere und L-Funktionen
# -----------------------------------------------------------------

def chi_trivial(n: int) -> complex:
    """Trivialer Charakter chi_0(n) = 1 fuer alle n."""
    return 1.0 + 0j


def chi_4(n: int) -> complex:
    """
    Nicht-trivialer Charakter mod 4.
    chi_4(n) = 0 falls n gerade, +1 falls n = 1 mod 4, -1 falls n = 3 mod 4.
    Erfuellt chi_4(-1) = -1 (odd character).
    """
    n_mod = n % 4
    if n_mod == 1:
        return 1.0 + 0j
    elif n_mod == 3:
        return -1.0 + 0j
    else:
        return 0.0 + 0j


def L_partial(s: complex, chi: Callable[[int], complex], n_terms: int = 200) -> complex:
    """
    Abgeschnittene L-Funktion L(s, chi) = sum_{n=1}^{n_terms} chi(n)/n^s.

    Fuer Re(s) > 1 schnell konvergent; fuer Re(s) = 1/2 brauchen wir
    mehr Terme (langsam konvergent), aber fuer die Pilot-Groessenordnung reicht 200.
    """
    total = 0.0 + 0j
    for n in range(1, n_terms + 1):
        total += chi(n) / n**s
    return total


# -----------------------------------------------------------------
# 2. Mellin-Stuetzstellen und Quadratur
# -----------------------------------------------------------------

# Symmetrische Stuetzstellen: t in [-T_MAX, T_MAX]
t_grid = np.linspace(-T_MAX, T_MAX, N_MELLIN)
dt = t_grid[1] - t_grid[0]

# Quadratur-Gewichte (einfache Trapez-Naeherung)
weights = np.full(N_MELLIN, dt)
weights[0] *= 0.5
weights[-1] *= 0.5


# -----------------------------------------------------------------
# 3. Psi_{lambda,chi} als Mellin-Multiplikation
# -----------------------------------------------------------------

def build_psi_matrix(chi: Callable[[int], complex]) -> np.ndarray:
    """
    Psi_{lambda,chi} im Mellin-Bild = Multiplikation mit L(s, chi) auf
    s = 1/2 + i*t fuer t in t_grid.

    Rueckgabe: Diagonalmatrix Psi mit Diagonaleintraegen L(1/2+i*t_k, chi).

    Bemerkung: Der Randoperator R_{lambda,chi} aus Satz 3.2 ist hier
    vernachlaessigt; das ist fuer das Pilot-Signal aber akzeptabel
    (es ist endlich-dim und klein).
    """
    diag = np.zeros(N_MELLIN, dtype=complex)
    for k, t in enumerate(t_grid):
        s = 0.5 + 1j * t
        diag[k] = L_partial(s, chi)
    return np.diag(diag)


# -----------------------------------------------------------------
# 4. PW_{lambda,chi} als Mellin-Differentialoperator
# -----------------------------------------------------------------

def build_PW_matrix() -> np.ndarray:
    """
    PW_lambda im Mellin-Bild ist approximativ -d^2/dt^2 + potential(t).
    Fuer das Pilot-Problem reicht es, -d^2/dt^2 (zweite Ableitung) als
    diskrete Finite-Differenzen-Matrix zu bauen.

    Begruendung: PW_lambda = D_log^2 + Gamma (CCM 2024 semilocal).
    D_log ist im Mellin-Bild Multiplikation mit i*t, also D_log^2 ist
    Multiplikation mit -t^2. Aber wenn wir in Mellin-Basis arbeiten,
    geht der Laplace-Operator ueber die inverse Fourier zur 2. Ableitung.

    Hier nehmen wir die einfachere Version: PW als Diagonal mit t^2 + 1
    (modellierende positive Definitheit).
    """
    # Einfache Modell-Version: PW_lambda repraesentiert durch t^2 + Konstante
    # im Mellin-Bild, wobei t die Mellin-Variable ist.
    return np.diag(t_grid**2 + 1.0)


# -----------------------------------------------------------------
# 5. QW_{lambda,chi} als Weil-Explizitformel-Matrix
# -----------------------------------------------------------------

def build_QW_matrix(chi: Callable[[int], complex]) -> np.ndarray:
    """
    QW_{lambda,chi} auf Mellin-Basis ist die Matrix des Weil-Kerns.

    Vereinfachtes Modell fuer Pilot:
        QW_{lambda,chi}(s,s') = delta(s-s') * K_arch(s)
                               - sum_{p prime, p <= lambda^2, p not | q}
                                   chi(p) * log(p)/sqrt(p)
                                   * cos((s-s') * log(p))

    Die Diagonale ist der archimedische Beitrag (aus Gamma-Ableitungen),
    die off-diagonal-Terme sind die Prim-Beitraege.
    """
    # Archimedischer Anteil (stark vereinfacht als t^2-Funktion)
    K_arch = np.diag(0.5 * t_grid**2)

    # Prim-Summe: alle Primen p <= lambda^2
    p_max = int(LAMBDA**2) + 1  # fuer lambda = sqrt(14): bis p=14
    primes = [p for p in range(2, p_max + 1) if all(p % q != 0 for q in range(2, int(np.sqrt(p))+1))]

    # Prim-Matrix
    K_prime = np.zeros((N_MELLIN, N_MELLIN), dtype=complex)
    for p in primes:
        # fuer chi_4: chi_4(2) = 0, also auch die Nullstelle bei p=2 relevant
        # (chi_4(p) = 0 wenn p|4, also p=2)
        weight = chi(p) * np.log(p) / np.sqrt(p)
        # Log-Kernkonstruktion: cos((s-s') log p)
        log_p = np.log(p)
        for k in range(N_MELLIN):
            for k2 in range(N_MELLIN):
                K_prime[k, k2] -= weight * np.cos((t_grid[k] - t_grid[k2]) * log_p)

    return K_arch + K_prime


# -----------------------------------------------------------------
# 6. Defekt-Norm berechnen
# -----------------------------------------------------------------

def compute_defect_norm(chi: Callable[[int], complex], chi_name: str) -> dict:
    """
    Berechne eps = ||QW . Psi - mu . Psi . PW|| fuer optimales mu.

    Rueckgabe: dict mit eps, mu, singular-values.
    """
    Psi = build_psi_matrix(chi)
    PW = build_PW_matrix()
    QW = build_QW_matrix(chi)

    # Kommutator-Fehler: fuer welches mu ist QW.Psi - mu.Psi.PW minimal?
    # Lineares Kleinste-Quadrate-Problem: minimiere ||A - mu B||^2
    # mit A = QW @ Psi, B = Psi @ PW.
    # mu = tr(A^* B) / tr(B^* B)  (Projektion auf span(B))
    A = QW @ Psi
    B = Psi @ PW

    tr_AB = np.trace(A.conj().T @ B)
    tr_BB = np.trace(B.conj().T @ B)

    if abs(tr_BB) < 1e-12:
        mu_opt = 1.0
    else:
        mu_opt = tr_AB / tr_BB

    # Defekt-Matrix und Operator-Norm
    defect = A - mu_opt * B
    defect_norm = np.linalg.norm(defect, ord=2)  # Spektralnorm
    frob_norm = np.linalg.norm(defect, ord='fro')

    # Normalisiere durch ||B|| fuer dimensionslose Groesse
    B_norm = np.linalg.norm(B, ord=2)
    relative_defect = defect_norm / max(B_norm, 1e-12)

    return {
        'chi_name': chi_name,
        'mu_opt': complex(mu_opt),
        'defect_spectral': float(defect_norm),
        'defect_frobenius': float(frob_norm),
        'B_norm': float(B_norm),
        'relative_defect': float(relative_defect),
    }


# -----------------------------------------------------------------
# 7. Main
# -----------------------------------------------------------------

def main():
    print("=" * 70)
    print("chi_defect_norm_pilot.py")
    print(f"lambda = {LAMBDA:.4f}, L = log(lambda) = {L:.4f}")
    print(f"Galerkin-Dim N = {N_MELLIN}, T_max = {T_MAX:.4f}")
    print("=" * 70)
    print()

    # Test 1: Riemann-Baseline (chi_0)
    res_0 = compute_defect_norm(chi_trivial, "chi_0 (trivial)")

    # Test 2: chi_4 (mod 4, odd)
    res_4 = compute_defect_norm(chi_4, "chi_4 (mod 4, odd)")

    # Report
    for res in [res_0, res_4]:
        print(f"Charakter: {res['chi_name']}")
        print(f"  mu_opt         = {res['mu_opt']}")
        print(f"  ||B|| (Psi PW) = {res['B_norm']:.4f}")
        print(f"  Defekt (spec)  = {res['defect_spectral']:.4f}")
        print(f"  Defekt (Frob)  = {res['defect_frobenius']:.4f}")
        print(f"  Relativ-Defekt = {res['relative_defect']:.4f}")
        print()

    # Vergleich
    print("-" * 70)
    print("Vergleich:")
    ratio = res_4['relative_defect'] / max(res_0['relative_defect'], 1e-12)
    print(f"  Relativ-Defekt(chi_4) / Relativ-Defekt(chi_0) = {ratio:.4f}")
    print()

    print("Interpretation:")
    if ratio > 2.0:
        print("  -> chi-spezifische Signatur vorhanden (Ratio > 2).")
        print("  -> Qualitativ konsistent mit Hypothesen H1/H2.")
    elif ratio > 1.2:
        print("  -> Leichte chi-Abhaengigkeit, aber schwach.")
    else:
        print("  -> Kaum chi-Signatur. Pilot-Version reicht nicht.")
        print("     Naechste Iteration: groesseres N, echter QW-Kern mit Nullstellen-Polen.")

    print()
    print("Pilot-Hinweis:")
    print("  Dies ist eine STARK vereinfachte Version. Die echten Operatoren")
    print("  QW_{lambda,chi} und PW_{lambda} erfordern:")
    print("  - Explizite Weil-Kern-Konstruktion mit Sonin-Zerlegung")
    print("  - Prolate-Spheroidal-Eigenfunktionen (nicht Mellin-Stuetzstellen)")
    print("  - Korrekte Galerkin-Basis (prolate wavefunctions)")
    print("  Fuer Pilot-Signal reicht die Mellin-Modellversion.")


if __name__ == "__main__":
    main()
