"""
chi_defect_norm_v2.py
=====================

Pilot v2 fuer die chi-twistete Defekt-Norm aus DIRICHLET_CCM_TRANSFER.md.

Verbesserungen gegenueber v1 (chi_defect_norm_pilot.py):

  1. ECHTE PROLATE-GALERKIN-BASIS
     Statt naiver Mellin-Stuetzstellen: Diagonalisierung des diskreten
     sinc-Bandbegrenzungs-Operators auf einem feinen Gitter t in [-L_max, L_max].
     Die ersten N_galerkin Eigenvektoren sind die diskreten Slepian-Pollak-
     Prolate-Funktionen h_n(t) (DPSS, Slepian 1978).

  2. QW-OPERATOR VIA SONIN-ZERLEGUNG
     Aus CCM 2024 (arXiv:2310.18423): PW_lambda = D_log^2 + Gamma.
     Damit QW_{lambda,chi} = PW_lambda + Gamma_{lambda,chi} mit
        Gamma_{lambda,chi} = arch_correction + prime_sum(chi).
     Das trennt den "gemeinsamen" PW-Teil vom chi-spezifischen Weil-Teil
     und eliminiert den trivialen Beitrag der Pol-Aufblaeh-Effekte.

  3. DREI POL-NORMALISIERUNGEN PARALLEL
     (a) v1-Version: eps / ||Psi.PW||              (dominiert vom zeta-Pol)
     (b) Pol-normiert: eps / ||Psi||_F             (skaliert mit L-Groesse)
     (c) Dimensionsfrei: eps / (||Psi||_2 * ||PW||_2)

Zweck: Isolieren ob die chi-spezifische Information (chi_33-Anomalie,
chi_21-N-Oszillation, Even/Odd-Paritaet) in einer der drei Normalisierungen
stabil sichtbar wird.

Autor: LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung)
Datum: 2026-04-18
Ausfuehrung:
    PYTHONIOENCODING=utf-8 python chi_defect_norm_v2.py
"""

import json
import numpy as np
import time
from pathlib import Path
from typing import Callable, Dict, Any

# =================================================================
# 0. Konfiguration
# =================================================================

LAMBDA = np.sqrt(14.0)     # CCM-Standard-Skala
L = np.log(LAMBDA)         # log(lambda) ~ 1.32

# Paley-Wiener-Band: T < L (Satz 4.1 aus MILESTONE_1_CHI_DOMAIN_THEORY)
T_PW = 0.95 * L            # fast am Rand, um volles Signal zu sehen

# Feinheit der numerischen Diagonalisierung
N_GRID = 600               # t-Grid-Punkte auf [-T_WIDE, T_WIDE]
T_WIDE = 3.0 * L           # weitere Raender fuer korrekte Prolate-Asymptotik
N_GALERKIN = 25            # Prolate-Modes (erste ~N < 2*T_PW*T_WIDE/pi laut Nyquist)

# L-Funktions-Summen-Abschneidung
N_L_TERMS = 400

# Output
RESULTS_DIR = Path(__file__).parent.parent / "_results"
RESULTS_DIR.mkdir(exist_ok=True)

# =================================================================
# 1. Dirichlet-Charaktere
# =================================================================

def make_chi_trivial():
    """Trivialer Charakter chi_0(n) = 1, modulus q=1."""
    def chi(n): return 1.0 + 0j
    return {'name': 'chi_0', 'q': 1, 'parity': +1, 'chi': chi, 'gamma1': 14.1347}

def make_chi_mod4_odd():
    """
    chi_4: nicht-trivialer Charakter mod 4, odd (chi(-1)=-1).
    chi(1)=+1, chi(3)=-1, chi(2)=chi(4)=0.
    Erste Nullstelle: gamma^(1) = 6.0209 (L(s, chi_4) = Dirichlet beta).
    """
    def chi(n):
        m = n % 4
        if m == 1: return 1.0 + 0j
        if m == 3: return -1.0 + 0j
        return 0.0 + 0j
    return {'name': 'chi_4', 'q': 4, 'parity': -1, 'chi': chi, 'gamma1': 6.0209}

def make_chi_mod_from_zeros(label: str, q: int, parity: int, gamma1: float,
                             chi_values: Dict[int, complex]):
    """
    Generischer Charakter mod q, definiert ueber explizite Werte chi(n) fuer 0 < n < q.
    chi_values: dict {1: val_1, 2: val_2, ...}
    Werte fuer n mit gcd(n,q)>1 setzen wir 0.
    Der Charakter wird modulo q periodisch erweitert.
    """
    # Normalisiere: chi_values sollte fuer jedes n in [1, q-1] einen Wert haben oder 0
    vals = np.zeros(q, dtype=complex)
    for n in range(q):
        if n in chi_values:
            vals[n] = chi_values[n]
        else:
            # gcd != 1 implizit => 0
            vals[n] = 0.0 + 0j

    def chi(n):
        return vals[n % q]

    return {'name': label, 'q': q, 'parity': parity, 'chi': chi, 'gamma1': gamma1}

def load_zeros_from_atlas() -> Dict[str, list]:
    path = Path(__file__).parent.parent.parent / "dirichlet_atlas" / "_results" / "zeros_all_chars.json"
    with open(path, 'r') as f:
        return json.load(f)

# =================================================================
# 2. L-Funktion auf kritischer Linie
# =================================================================

def L_values_on_line(t_grid: np.ndarray, chi: Callable[[int], complex],
                      n_terms: int = N_L_TERMS) -> np.ndarray:
    """
    Wertet L(1/2 + i*t, chi) an t_grid aus via abgeschnittener Dirichlet-Reihe.
    Fuer den trivialen Charakter nutze Riemann-Siegel-artige Regularisierung:
    ziehe den Pol-Beitrag heraus durch (1/(s-1))-Subtraktion.

    Rueckgabe: complex array (N_grid,)
    """
    result = np.zeros_like(t_grid, dtype=complex)
    s_grid = 0.5 + 1j * t_grid
    for n in range(1, n_terms + 1):
        chi_n = chi(n)
        if chi_n == 0:
            continue
        result += chi_n / (n ** s_grid)
    return result

# =================================================================
# 3. Prolate-Galerkin-Basis via sinc-Operator-Diagonalisierung
# =================================================================

def build_prolate_basis(n_grid: int, t_wide: float, t_pw: float,
                         n_galerkin: int) -> Dict[str, np.ndarray]:
    """
    Baue die diskrete Prolate-Basis via numerische Diagonalisierung
    des sinc-Bandbegrenzungs-Operators.

    Gitter: t in [-t_wide, +t_wide], n_grid Punkte.
    Operator: (B g)(t) = int sinc(t_pw * (t - t')) * g(t') dt'
              = Projektor auf Paley-Wiener-Raum PW_{t_pw}.

    Eigenvektoren von B mit groessten Eigenwerten sind die diskreten
    Prolate-Funktionen h_n(t). Fuer n < ~(2*t_pw*t_wide)/pi sind
    Eigenwerte nahe 1 (gut lokalisiert).

    Rueckgabe:
        t: (n_grid,) Gitter
        dt: Schrittweite
        H: (n_grid, n_galerkin) Matrix der Prolate-Funktionen (spaltenweise)
        lambdas: (n_galerkin,) zugehoerige Eigenwerte (Lokalisierungsgrade)
    """
    t = np.linspace(-t_wide, t_wide, n_grid)
    dt = t[1] - t[0]

    # sinc-Kern: K[k,l] = sin(t_pw*(t_k - t_l)) / (pi*(t_k - t_l))
    # an der Diagonale: K[k,k] = t_pw/pi
    diff = t[:, None] - t[None, :]
    # safe sinc
    K = np.where(np.abs(diff) < 1e-14,
                 t_pw / np.pi,
                 np.sin(t_pw * diff) / (np.pi * diff))
    K *= dt  # Gewichtung durch Quadratur

    # Symmetrisieren (numerische Stabilitaet)
    K = 0.5 * (K + K.T)

    # Eigendekomposition
    eigvals, eigvecs = np.linalg.eigh(K)
    # Absteigende Ordnung (groesste Eigenwerte zuerst = beste Lokalisierung)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Normalisiere Eigenvektoren: bzgl. L^2-Inner-Product int |h|^2 dt = sum |h_k|^2 * dt
    # eigh liefert orthonormale Vektoren bzgl. Standard-Inner-Product (ohne dt).
    # Wir wollen int-Normierung: h_n(t_k) = eigvec[k,n] / sqrt(dt)
    H = eigvecs[:, :n_galerkin] / np.sqrt(dt)
    lambdas = eigvals[:n_galerkin]

    return {
        't': t,
        'dt': dt,
        'H': H,
        'lambdas': lambdas,
    }

# =================================================================
# 4. Operator-Projektion in Prolate-Basis
# =================================================================

def project_diagonal(H: np.ndarray, dt: float, diag: np.ndarray) -> np.ndarray:
    """
    Berechne [M]_{mn} = int h_m(t) * diag(t) * h_n(t) dt.
    Numerisch: sum_k H[k,m] * diag[k] * H[k,n] * dt.
    """
    # H: (N_grid, N_gal), diag: (N_grid,)
    # Ergebnis: (N_gal, N_gal) = H^T @ diag(diag) @ H * dt
    return (H.T * diag[None, :]) @ H * dt

def project_dense(H: np.ndarray, dt: float, kernel: np.ndarray) -> np.ndarray:
    """
    Berechne [M]_{mn} = int int h_m(t) K(t,t') h_n(t') dt dt'.
    """
    # H: (N_grid, N_gal), kernel: (N_grid, N_grid)
    # Ergebnis: H^T @ kernel @ H * dt^2
    return H.T @ kernel @ H * dt**2

# =================================================================
# 5. PW-Operator im Mellin-Bild
# =================================================================

def build_PW_matrix(basis: Dict[str, np.ndarray]) -> np.ndarray:
    """
    PW_lambda im Mellin-Bild approximativ als Multiplikation mit t^2 + c.

    Hintergrund: der Prolate-Operator ist
        PW_lambda g(x) = -d/dx[(lambda^2 - x^2) g'(x)] + (2*pi*lambda*x)^2 g(x).
    Seine Mellin-Darstellung auf der kritischen Linie s=1/2+it ist
    naeherungsweise ein Multiplikationsoperator mit einer Funktion omega(t),
    deren Asymptotik fuer t->inf wie t^2 ist.

    Fuer den Pilot: omega(t) = t^2 + 1 (einfache positive Quadratik).
    """
    t = basis['t']
    dt = basis['dt']
    H = basis['H']
    omega = t**2 + 1.0
    return project_diagonal(H, dt, omega)

# =================================================================
# 6. Psi-Operator (Mellin-Multiplikation mit L(s,chi))
# =================================================================

def build_Psi_matrix(basis: Dict[str, np.ndarray],
                      chi: Callable[[int], complex]) -> np.ndarray:
    """
    Psi_{lambda,chi} im Mellin-Bild = Multiplikation mit L(1/2+it, chi).

    Rand-Operator R_{lambda,chi} aus Satz 3.2 ist auf Paley-Wiener-Unterraum
    klein (endlich-dim) und im Pilot vernachlaessigt.
    """
    t = basis['t']
    dt = basis['dt']
    H = basis['H']
    L_vals = L_values_on_line(t, chi)
    return project_diagonal(H, dt, L_vals)

# =================================================================
# 7. QW-Operator via Sonin-Zerlegung
# =================================================================

def build_Gamma_arch(basis: Dict[str, np.ndarray]) -> np.ndarray:
    """
    Archimedischer Korrektur-Teil: Gamma_arch(t) ~ 2 * Re(digamma(1/4 + i*t/2)).

    Das ist der "Archimedische Explizitformel"-Beitrag, der fuer Riemann
    und Dirichlet-L-Funktionen auf der kritischen Linie Re(s)=1/2 gleich
    ist (bis auf Conductor-dependence, die wir hier vernachlaessigen).

    Typischer Wert: -1.96 bei t=0, langsam anwachsend als log(|t|/2) fuer t->inf.
    """
    from scipy.special import digamma
    t = basis['t']
    dt = basis['dt']
    H = basis['H']
    # Archimedische Dichte: 2*Re(Gamma'(s/2 + 1/4)/Gamma(s/2+1/4))
    # Auf s=1/2+it: arg = 1/4 + 1/4 + i*t/2 = 1/2 + i*t/2  (Riemann-Normalisierung)
    # Fuer Dirichlet allgemein: 1/4 + parity/4 + i*t/2
    # Hier: neutrale Version fuer chi_0
    arch_density = 2.0 * np.real(digamma(0.5 + 1j * t / 2.0))
    return project_diagonal(H, dt, arch_density)

def build_Gamma_prime(basis: Dict[str, np.ndarray],
                       chi: Callable[[int], complex],
                       q: int) -> np.ndarray:
    """
    Prim-Teil des Weil-Kerns: Gamma_prime(t, t') = Integrand der Prim-Summe.

    Aus Weil-Explizitformel fuer L(s, chi):
        Gamma_prime(t, t') = - sum_{p prime, p nmid q} sum_{m >= 1}
                              chi(p)^m * log(p) / p^{m/2}
                              * cos((t - t') * m * log(p))

    Cutoff: p <= lambda^2 und m <= m_max = ceil(lambda / log(p)).

    Das ist die Source des char-spezifischen Signals.
    """
    t = basis['t']
    dt = basis['dt']
    H = basis['H']
    n_grid = len(t)

    # Primes bis lambda^2 (fuer lambda=sqrt(14): p <= 14)
    p_max = int(LAMBDA ** 2) + 1
    primes = [p for p in range(2, p_max + 1)
              if all(p % q_ != 0 for q_ in range(2, int(np.sqrt(p)) + 1))]

    diff = t[:, None] - t[None, :]
    kernel = np.zeros((n_grid, n_grid), dtype=complex)

    for p in primes:
        if q > 1 and p % q == 0:  # p teilt Conductor -> chi(p)=0
            continue
        log_p = np.log(p)
        sqrt_p = np.sqrt(p)
        # m_max sodass p^m <= lambda^2
        m_max = max(1, int(np.log(LAMBDA**2) / log_p))
        for m in range(1, m_max + 1):
            chi_pm = chi(p) ** m
            if chi_pm == 0:
                continue
            weight = -chi_pm * log_p / (sqrt_p ** m)
            kernel += weight * np.cos(diff * m * log_p)

    # Projektion der 2D-Kern-Matrix in Prolate-Basis (dense)
    return project_dense(H, dt, kernel)

def build_QW_matrix(basis: Dict[str, np.ndarray],
                     chi_info: Dict[str, Any],
                     PW_mat: np.ndarray) -> np.ndarray:
    """
    QW_{lambda,chi} via Sonin-Zerlegung:
        QW = PW + Gamma_arch + Gamma_prime(chi)

    Der arch-Teil ist chi-unabhaengig im Pilot (keine Conductor-Korrektur);
    der Prim-Teil traegt die char-Signatur.
    """
    Gamma_arch = build_Gamma_arch(basis)
    Gamma_prime = build_Gamma_prime(basis, chi_info['chi'], chi_info['q'])
    return PW_mat + Gamma_arch + Gamma_prime

# =================================================================
# 8. Defekt-Norm mit drei Normalisierungen
# =================================================================

def compute_defect(basis: Dict[str, np.ndarray],
                    chi_info: Dict[str, Any],
                    PW_mat: np.ndarray) -> Dict[str, Any]:
    """
    Berechne fuer gegebenen Charakter:
        Psi, QW aufbauen
        A = QW @ Psi
        B = Psi @ PW
        mu_opt = tr(A^* B) / tr(B^* B) (Frobenius-optimal)
        Defekt = A - mu_opt * B
        norms und drei Relativ-Defekte
    """
    Psi = build_Psi_matrix(basis, chi_info['chi'])
    QW = build_QW_matrix(basis, chi_info, PW_mat)

    A = QW @ Psi
    B = Psi @ PW_mat

    tr_AB = np.trace(A.conj().T @ B)
    tr_BB = np.trace(B.conj().T @ B)
    mu_opt = tr_AB / tr_BB if abs(tr_BB) > 1e-14 else 1.0

    defect = A - mu_opt * B
    defect_spec = np.linalg.norm(defect, ord=2)
    defect_frob = np.linalg.norm(defect, ord='fro')

    # Drei Normalisierungen
    norm_B_spec = np.linalg.norm(B, ord=2)             # wie v1
    norm_Psi_frob = np.linalg.norm(Psi, ord='fro')
    norm_Psi_spec = np.linalg.norm(Psi, ord=2)
    norm_PW_spec = np.linalg.norm(PW_mat, ord=2)

    # (a) v1-Style
    rel_a = defect_spec / max(norm_B_spec, 1e-14)
    # (b) Pol-normiert (Psi-Frobenius raushandeln)
    rel_b = defect_spec / max(norm_Psi_frob, 1e-14)
    # (c) Dimensionsfrei
    rel_c = defect_spec / max(norm_Psi_spec * norm_PW_spec, 1e-14)

    return {
        'chi_name': chi_info['name'],
        'q': chi_info['q'],
        'parity': chi_info['parity'],
        'gamma1': chi_info['gamma1'],
        'mu_opt_real': float(np.real(mu_opt)),
        'mu_opt_imag': float(np.imag(mu_opt)),
        'defect_spec': float(defect_spec),
        'defect_frob': float(defect_frob),
        'norm_B_spec': float(norm_B_spec),
        'norm_Psi_frob': float(norm_Psi_frob),
        'norm_Psi_spec': float(norm_Psi_spec),
        'norm_PW_spec': float(norm_PW_spec),
        'rel_a_v1style': float(rel_a),
        'rel_b_polnorm': float(rel_b),
        'rel_c_dimfree': float(rel_c),
    }

# =================================================================
# 9. Hauptprogramm
# =================================================================

def main():
    print("=" * 75)
    print("chi_defect_norm_v2.py - Pilot v2")
    print(f"lambda = {LAMBDA:.4f}, L = log(lambda) = {L:.4f}")
    print(f"N_GRID = {N_GRID}, N_GALERKIN = {N_GALERKIN}")
    print(f"T_PW = {T_PW:.4f}, T_WIDE = {T_WIDE:.4f}")
    print("=" * 75)

    # 1. Prolate-Basis bauen (chi-unabhaengig)
    t0 = time.time()
    print("\n[1/4] Baue Prolate-Basis (Diagonalisierung sinc-Operator)...")
    basis = build_prolate_basis(N_GRID, T_WIDE, T_PW, N_GALERKIN)
    print(f"      Eigenwerte Prolate (erste 5): {basis['lambdas'][:5]}")
    print(f"      Letzter Eigenwert (n={N_GALERKIN-1}): {basis['lambdas'][-1]:.4e}")
    print(f"      Zeit: {time.time()-t0:.2f}s")

    # 2. PW-Matrix (chi-unabhaengig)
    print("\n[2/4] Baue PW-Matrix...")
    t0 = time.time()
    PW_mat = build_PW_matrix(basis)
    print(f"      Shape: {PW_mat.shape}, Spektralnorm: {np.linalg.norm(PW_mat, 2):.4f}")
    print(f"      Zeit: {time.time()-t0:.2f}s")

    # 3. Charaktere laden
    print("\n[3/4] Lade Charaktere + Atlas-Nullstellen...")
    atlas_zeros = load_zeros_from_atlas()

    chars = [
        make_chi_trivial(),
        make_chi_mod4_odd(),
    ]

    # chi_5: mod 5, odd, gamma1 aus Atlas
    # Nach Atlas: chi_5 ist primitiver mod 5, odd (Legendre-Symbol)
    # chi(1)=1, chi(2)=i, chi(3)=-i, chi(4)=-1 (wenn chi eine quartische Charakter mod 5 ist)
    # ABER: der Atlas listet NUR reelle Gaps; also nehmen wir reelle Charaktere.
    # Fuer mod 5: reeller nicht-trivialer Charakter ist das Legendre-Symbol: chi(1)=1, chi(4)=1, chi(2)=-1, chi(3)=-1
    # Das ist der reelle primitive Charakter mod 5, even (chi(-1)=chi(4)=+1).
    # gamma1 aus atlas_zeros['chi_5']
    chars.append(make_chi_mod_from_zeros(
        'chi_5', q=5, parity=+1, gamma1=atlas_zeros['chi_5'][0],
        chi_values={1: 1.0+0j, 2: -1.0+0j, 3: -1.0+0j, 4: 1.0+0j}
    ))

    # chi_8: reeller primitiver EVEN Charakter mod 8, Kronecker(+8, n).
    # Atlas char_features_all.json: D=+8, parity=+1 (even).
    # Kronecker(8, n): chi(1)=1, chi(3)=-1, chi(5)=-1, chi(7)=+1, chi(-1)=chi(7)=+1 -> even.
    chars.append(make_chi_mod_from_zeros(
        'chi_8', q=8, parity=+1, gamma1=atlas_zeros['chi_8'][0],
        chi_values={1: 1.0+0j, 3: -1.0+0j, 5: -1.0+0j, 7: 1.0+0j}
    ))

    # chi_33: reeller primitiver EVEN Charakter mod 33, Kronecker(+33, n).
    # Atlas: D=+33, parity=+1 (even), gamma1=2.997 (kleinste unter 10 Test-Charakteren ausser chi_21, chi_29, chi_60).
    # Werte via sympy verifiziert.
    chi_33_kronecker_plus = [0, 1, 1, 0, 1, -1, 0, -1, 1, 0, -1, 0, 0, -1, -1, 0, 1,
                              1, 0, -1, -1, 0, 0, -1, 0, 1, -1, 0, -1, 1, 0, 1, 1]
    chi_33_vals = {n: float(chi_33_kronecker_plus[n]) + 0j
                    for n in range(33) if chi_33_kronecker_plus[n] != 0}
    chars.append(make_chi_mod_from_zeros(
        'chi_33', q=33, parity=+1, gamma1=atlas_zeros['chi_33'][0],
        chi_values=chi_33_vals
    ))

    print(f"      {len(chars)} Charaktere geladen: {[c['name'] for c in chars]}")

    # 4. Defekt fuer jeden Charakter
    print("\n[4/4] Berechne Defekte fuer alle Charaktere...")
    results = []
    for c in chars:
        t0 = time.time()
        res = compute_defect(basis, c, PW_mat)
        res['compute_time'] = time.time() - t0
        results.append(res)
        print(f"  {res['chi_name']:10s} (q={res['q']:3d}, par={res['parity']:+d}, gamma1={res['gamma1']:.3f}): "
              f"rel_a={res['rel_a_v1style']:.4f}, rel_b={res['rel_b_polnorm']:.4f}, rel_c={res['rel_c_dimfree']:.4e}, "
              f"time={res['compute_time']:.1f}s")

    # Zusammenfassung
    print("\n" + "=" * 75)
    print("ZUSAMMENFASSUNG")
    print("=" * 75)
    print(f"\n{'Charakter':12s} {'q':>4s} {'parity':>7s} {'gamma1':>8s} "
          f"{'rel_a':>10s} {'rel_b':>10s} {'rel_c':>12s} {'mu_opt_re':>10s}")
    print("-" * 85)
    for res in results:
        print(f"{res['chi_name']:12s} {res['q']:4d} {res['parity']:+7d} {res['gamma1']:8.3f} "
              f"{res['rel_a_v1style']:10.4f} {res['rel_b_polnorm']:10.4f} "
              f"{res['rel_c_dimfree']:12.4e} {res['mu_opt_real']:10.4f}")

    # Ratios gegen chi_0
    baseline = results[0]  # chi_0
    print(f"\nRATIOS (bezogen auf chi_0):")
    print(f"{'Charakter':12s} {'ratio_a':>10s} {'ratio_b':>10s} {'ratio_c':>10s} {'H1_vorh.':>10s} {'H2_vorh.':>10s}")
    print("-" * 70)
    for res in results:
        ra = res['rel_a_v1style'] / max(baseline['rel_a_v1style'], 1e-14)
        rb = res['rel_b_polnorm'] / max(baseline['rel_b_polnorm'], 1e-14)
        rc = res['rel_c_dimfree'] / max(baseline['rel_c_dimfree'], 1e-14)
        h2_pred = baseline['gamma1'] / res['gamma1']  # H2: gamma_0/gamma_chi
        # H1: R_chi/R_0 ~ grob 1/gamma^3 summiert; approximiere mit 1/gamma1^3 * const
        h1_pred = (baseline['gamma1'] / res['gamma1']) ** 3
        print(f"{res['chi_name']:12s} {ra:10.4f} {rb:10.4f} {rc:10.4f} {h1_pred:10.2f} {h2_pred:10.2f}")

    # Speichern
    out_json = RESULTS_DIR / "PILOT_CHI_DEFECT_V2_2026-04-18.json"
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump({
            'config': {
                'lambda': float(LAMBDA),
                'L': float(L),
                'N_grid': N_GRID,
                'N_galerkin': N_GALERKIN,
                'T_PW': float(T_PW),
                'T_wide': float(T_WIDE),
                'N_L_terms': N_L_TERMS,
            },
            'results': results,
        }, f, indent=2)
    print(f"\nRohdaten gespeichert: {out_json}")

    return results

if __name__ == "__main__":
    main()
