"""
dirichlet_ccm_mpmath.py
========================

VOLLE CCM-2025-Implementation mit mpmath (hochpraezise Arithmetik).

Implementiert die praezise arch-Matrix nach CCM 2025 Prop 4.3:

  W_R(V_n, V_m) = (alpha_L(m) - alpha_L(n))/(n - m)     wenn n != m
                = 2*gamma_L(n) - 2*beta_L(n)           wenn n == m

mit:
  alpha_L(n) = (1/pi) * int_0^L sin(2*pi*n*x/L) * rho(x) dx
  beta_L(n)  = (1/L)  * int_0^L x * cos(2*pi*n*x/L) * rho(x) dx
  gamma_L(n) = int_0^L (cos(2*pi*n*x/L) - exp(-x/2)) * rho(x) dx + c(L) + w(L)

  rho(x) = exp(x/2) / (exp(x) - exp(-x)) = e^{x/2}/(2 sinh x)

Alle Integrale numerisch via mpmath.quad mit 50 Dezimalstellen.

Riemann-Sanity: Nullstellen bei lambda=sqrt(14) mit N=30 sollten
deutlich besser als die ersten Implementation sein (erste Nullstelle sollte
nahe 14.13 sein, nicht 1.23).

Autor: LG (Opus 4.7 [1M], Session 16 Fortsetzung, Priorität α)
Datum: 2026-04-18
Ref: Connes-Consani-Moscovici 2025 (arXiv:2511.22755)
"""

import numpy as np
import mpmath as mp
import json
from pathlib import Path

mp.mp.dps = 50  # 50 Dezimalstellen Genauigkeit

# =================================================================
# 0. Konfiguration
# =================================================================

LAMBDA = mp.sqrt(mp.mpf(14))
N_DEFAULT = int(__import__('os').environ.get('CCM_N', 120))
# CCM-Konvention nach Eq (3.10) + Eq (4.1): Psi = W_{0,2} - W_R - sum W_p
# => M = -WR - Wp + W02, also SIGN_WR = -1 ist mathematisch korrekt.
SIGN_WR_DEFAULT = int(__import__('os').environ.get('CCM_SIGN_WR', -1))

RESULTS_DIR = Path(__file__).parent.parent / "_results"

RIEMANN_ZEROS = ['14.134725141734693790457251983562470270784257115699',
                  '21.022039638771554992628479593896902777334340524903',
                  '25.010857580145688763213790992562821818659549672558',
                  '30.424876125859513210311897530584091320181560023715',
                  '32.935061587739189690662368964074903488812715603517',
                  '37.586178158825671257217763480705332821405597350831',
                  '40.918719012147495187398126914633254395726165962777',
                  '43.327073280914999519496122165406805782645668371837',
                  '48.005150881167159727942472749427516041686844001144',
                  '49.773832477672302181916784678563724057723178299677']

CHI4_ZEROS_STR = ['6.020948406467813248841013122025116830842932188530',
                   '10.243766565498344019708550001906014923988016843716',
                   '12.988096097420338478712946059893727769149892049866',
                   '16.343297172788686238881030410923648213619985927091',
                   '18.291996171925028354003894030404571716663018006019']

# =================================================================
# 1. Grundfunktion rho
# =================================================================

def rho(x):
    """rho(x) = e^{x/2} / (e^x - e^{-x}) = e^{x/2}/(2 sinh x).

    Bei x->0 divergent wie 1/(2x), aber wir integrieren gegen Funktionen
    die bei x=0 Nullstellen haben (sin(2*pi*n*x/L) bzw (cos-1)).
    """
    return mp.exp(x/2) / (mp.exp(x) - mp.exp(-x))

# =================================================================
# 2. alpha_L(n), beta_L(n), gamma_L(n) via mpmath-Integration
# =================================================================

def alpha_L(n, L):
    """(1/pi) * int_0^L sin(2*pi*n*x/L) * rho(x) dx."""
    if n == 0:
        return mp.mpf(0)
    def integrand(x):
        return mp.sin(2*mp.pi*n*x/L) * rho(x)
    integral = mp.quad(integrand, [0, L])
    return integral / mp.pi

def beta_L(n, L):
    """(1/L) * int_0^L x * cos(2*pi*n*x/L) * rho(x) dx."""
    def integrand(x):
        return x * mp.cos(2*mp.pi*n*x/L) * rho(x)
    integral = mp.quad(integrand, [0, L])
    return integral / L

def gamma_L_func(n, L, c_plus_w):
    """int_0^L (cos(2*pi*n*x/L) - exp(-x/2)) * rho(x) dx + c(L) + w(L)."""
    def integrand(x):
        return (mp.cos(2*mp.pi*n*x/L) - mp.exp(-x/2)) * rho(x)
    integral = mp.quad(integrand, [0, L])
    return integral + c_plus_w

def compute_c_plus_w(L):
    """
    c(L) + w(L), zusammen per CCM Seite 14/15 (nach Vereinfachung):
    c(L) + w(L) = (1/2) log((e^{L/2}-1)/(e^{L/2}+1)) + atan(e^{L/2})
                  - pi/4 + gamma/2 + (1/2) log(8*pi)
    """
    eL2 = mp.exp(L/2)
    eL  = mp.exp(L)
    gamma_euler = mp.euler
    term1 = mp.mpf('0.5') * mp.log((eL2 - 1)/(eL2 + 1))
    term2 = mp.atan(eL2)
    term3 = -mp.pi/4
    term4 = gamma_euler/2
    term5 = mp.mpf('0.5') * mp.log(8*mp.pi)
    return term1 + term2 + term3 + term4 + term5

# =================================================================
# 3. Matrix-Beitraege
# =================================================================

def build_WR_matrix(N, L):
    """
    Archimedische Weil-Matrix via CCM Prop 4.3.
    W_R[n+N, m+N] = ((alpha_L(m) - alpha_L(n))/(n - m)) fuer n != m
                  = 2*gamma_L(n) - 2*beta_L(n) fuer n == m

    Gibt numpy-Array mit mpmath-Werten (komplex aber reell-symmetrisch).
    """
    print(f"  Berechne alpha/beta/gamma fuer n in [-{N}, {N}]...")
    c_plus_w = compute_c_plus_w(L)
    alpha_vals = {}
    beta_vals = {}
    gamma_vals = {}
    for n in range(-N, N+1):
        alpha_vals[n] = alpha_L(n, L)
        beta_vals[n] = beta_L(n, L)
        gamma_vals[n] = gamma_L_func(n, L, c_plus_w)

    size = 2*N + 1
    WR = np.zeros((size, size))
    for i in range(size):
        n = i - N
        for j in range(size):
            m = j - N
            if n == m:
                val = 2*gamma_vals[n] - 2*beta_vals[n]
            else:
                val = (alpha_vals[m] - alpha_vals[n]) / (n - m)
            WR[i, j] = float(val)  # Konvertiere zu double fuer numpy
    return WR

def build_W02_matrix(N, L):
    """CCM Eq (4.2): Rang-2-Pol-Matrix (nur Riemann)."""
    size = 2*N + 1
    W02 = np.zeros((size, size))
    L_mp = mp.mpf(L) if not isinstance(L, mp.mpf) else L
    L_float = float(L_mp)
    pre = 32.0 * L_float * np.sinh(L_float/4.0)**2
    for i in range(size):
        n = i - N
        for j in range(size):
            m = j - N
            num = L_float**2 - 16.0*np.pi**2*m*n
            den = (L_float**2 + 16.0*np.pi**2*m**2) * (L_float**2 + 16.0*np.pi**2*n**2)
            W02[i, j] = pre * num / den
    return W02

def q_kernel(m, n, y, L):
    """q(U_m, U_n)(y) fuer y in [0, L] per CCM Eq (2.9)."""
    if m == n:
        return 2.0 * (1.0 - y/L) * np.cos(2.0*np.pi*n*y/L)
    else:
        return (np.sin(2.0*np.pi*m*y/L) - np.sin(2.0*np.pi*n*y/L)) / (np.pi*(n-m))

def is_prime_power(k):
    if k < 2: return None, 0
    for p in range(2, int(np.sqrt(k)) + 1):
        if k % p == 0:
            m = 0; r = k
            while r % p == 0: r //= p; m += 1
            if r == 1: return p, m
            else: return None, 0
    return k, 1

def mangoldt(k):
    p, m = is_prime_power(k)
    return np.log(p) if p else 0.0

def build_Wprime_matrix(N, L, lam, chi_func, q_mod):
    """
    Prim-Beitrag per CCM Eq (4.3):
      sum W_p(V_n, V_m) = sum_{1<k<=exp(L)} Lambda(k) k^{-1/2} q(U_n, U_m)(log k)

    Fuer Dirichlet mit chi-Gewichten.
    """
    L_float = float(L) if isinstance(L, mp.mpf) else L
    size = 2*N + 1
    Wp = np.zeros((size, size))
    k_max = int(np.exp(L_float))  # = lambda^2
    for k in range(2, k_max + 1):
        lam_k = mangoldt(k)
        if lam_k == 0: continue
        if q_mod > 1:
            p, _ = is_prime_power(k)
            if p and q_mod % p == 0: continue
        chi_k = chi_func(k)
        if chi_k == 0: continue
        y = np.log(k)
        for i in range(size):
            n = i - N
            for j in range(size):
                m = j - N
                Wp[i, j] += chi_k * lam_k / np.sqrt(k) * q_kernel(m, n, y, L_float)
    return Wp

# =================================================================
# 4. Charaktere
# =================================================================

def chi_trivial(n): return 1.0
def chi_4(n):
    m = n % 4
    return 1.0 if m==1 else (-1.0 if m==3 else 0.0)
def chi_33(n):
    vals = [0,1,1,0,1,-1,0,-1,1,0,-1,0,0,-1,-1,0,1,1,0,-1,-1,0,0,-1,0,1,-1,0,-1,1,0,1,1]
    return float(vals[n % 33])

# =================================================================
# 5. QW-Matrix gesamt
# =================================================================

def build_QW(N, L, lam, chi_func, q_mod, parity, include_W02=False, cache_tag=None, sign_WR=-1):
    """QW_lambda^{chi,N} = sign_WR * W_R - sum W_p [+ W_02 for chi_0]."""
    cache_path = None
    if cache_tag is not None:
        cache_dir = Path(__file__).parent / "_cache"
        cache_dir.mkdir(exist_ok=True)
        cache_path = cache_dir / f"WR_N{N}_L{float(L):.6f}_{cache_tag}.npy"
    if cache_path is not None and cache_path.exists():
        print(f"  Lade W_R aus Cache: {cache_path.name}")
        WR = np.load(cache_path)
    else:
        print(f"  Baue W_R (arch, praezise Prop 4.3)...")
        WR = build_WR_matrix(N, L)
        if cache_path is not None:
            np.save(cache_path, WR)
            print(f"  Cache gespeichert: {cache_path.name}")
    # Fuer non-trivial chi: log(q/pi)-Offset + parity-Korrektur (noch nicht implementiert hier)
    # Dies ist Version 1: wir nehmen WR wie Riemann + signifikante log(q/pi)-Korrektur
    if q_mod > 1:
        # Rank-1 Conductor-Korrektur: log(q/pi) * |V_n~(0)|^2
        # V_n~(0) naeherungsweise: fuer V_n(u) = L^{-1/2} exp(2*pi*i*n*log(lam*u)/L):
        # V_n~(1/2) ist die Mellin an s=1/2, nicht trivial direkt
        # Vereinfacht: addiere log(q/pi) als globalen Diagonal-Shift
        conductor = float(mp.log(mp.mpf(q_mod)/mp.pi))
        # Rank-1: mit h_0 = delta_N / L^{1/2} = sum_j V_j / L, also in V-Basis
        # h_0 hat Eintraege 1/sqrt(L) in allen Positionen
        # Eigentlich: int h_n(t) * 1 dt = delta_{n,0} * L^{1/2} (constant ist V_0)
        # Vereinfacht hier: nur V_0-Diagonal bekommt conductor * L
        L_float = float(L)
        WR[N, N] += conductor * L_float  # Diagonalbeitrag auf V_0

    print(f"  Baue W_prime (Prim-Summe)...")
    Wp = build_Wprime_matrix(N, L, lam, chi_func, q_mod)

    M = sign_WR * WR - Wp
    if include_W02:
        M = M + build_W02_matrix(N, L)
    # Symmetrisieren
    M = 0.5 * (M + M.T)
    return M

# =================================================================
# 6. Grundzustand + Spektrum
# =================================================================

def find_ground_state(T, target_parity='even', degen_tol=1e-4, verbose=True):
    """Grundzustand via Projektion auf entarteten Unterraum.

    CCM Thm 5.10: xi ist einfacher Eigenvektor zum kleinsten EW mit delta_N(xi) != 0.

    Problem in der Numerik: der kleinste EW kann numerisch entartet sein (cluster
    von 2..20 EWs innerhalb tol). scipy.eigh liefert dann willkuerliche Linear-
    kombinationen, meist parity-gemischt.

    Loesung: im entarteten Unterraum D (spanned von EVs mit |w-w_min|<tol)
    konstruiere xi explizit als parity-symmetrisierte Projektion von delta_N
    (fuer even-chi) bzw. delta_N^- = sum sgn(j)V_j (fuer odd-chi).

    Falls die Projektion nicht-trivial ist (Norm > eps), ist das der echte xi.
    Falls sie null ist, existiert kein parity-korrekter Zustand mit delta_N(xi)!=0
    im Grundzustand-Cluster -- das waere eine MS1-Verletzung.
    """
    from scipy.linalg import eigh
    w, V = eigh(T)
    M = V.shape[0]
    N_half = (M - 1) // 2
    j_weights = np.arange(-N_half, N_half + 1, dtype=float)

    w_min = w[0]
    degen_idx = [k for k in range(len(w)) if abs(w[k] - w_min) < degen_tol]
    d = len(degen_idx)
    V_deg = V[:, degen_idx]                   # (M, d)

    if verbose:
        print(f"  find_ground_state: {target_parity}, w_min={w_min:.8f}, Cluster-Groesse d={d} (tol={degen_tol})")
        print(f"    Top-10 EWs: " + ", ".join(f"{w[k]:+.6f}" for k in range(min(10, len(w)))))

    # delta_N in V_n-Basis ist ~ (1,...,1)/sqrt(L); delta_N^- ~ sgn(j) bzw. j
    if target_parity == 'even':
        target = np.ones(M)
    else:
        # odd-projiziert: konsistent mit sum_j j*xi (Impuls-artige Kopplung)
        target = j_weights.copy()
    target = target / np.linalg.norm(target)

    # Projektion in den entarteten Unterraum
    coeffs = V_deg.T @ target                 # (d,)
    proj_norm = np.linalg.norm(coeffs)
    xi_raw = V_deg @ coeffs                   # (M,)

    # Parity-Symmetrisierung (kollabiert Mixed-Anteile)
    rev = xi_raw[::-1]
    if target_parity == 'even':
        xi_sym = 0.5 * (xi_raw + rev)
    else:
        xi_sym = 0.5 * (xi_raw - rev)
    sym_norm = np.linalg.norm(xi_sym)

    if verbose:
        print(f"    Projektion von delta_N ({target_parity}) in Cluster: ||P||={proj_norm:.4f}, nach parity-sym: {sym_norm:.4f}")

    if sym_norm < 1e-10:
        if verbose:
            print(f"    WARNUNG: parity-korrekte Projektion ~0 im Cluster -- MS1-Verletzung? Fallback auf k=0.")
        return w_min, V[:, 0], 0

    xi_final = xi_sym / sym_norm
    # Rayleigh-Quotient als Eigenwert-Verifikation (sollte ~ w_min sein)
    eps = float(xi_final @ T @ xi_final)

    if verbose:
        sum_xi = float(np.sum(xi_final))
        sum_j = float(np.sum(j_weights * xi_final))
        print(f"    xi konstruiert: eps={eps:.8f} (vs w_min={w_min:.8f}, Delta={eps-w_min:.2e})")
        print(f"    sum(xi)={sum_xi:+.4e}, sum(j*xi)={sum_j:+.4e}")

    return eps, xi_final, 0

def F_z(z, xi, N, L):
    """F(z) = sum_j xi_j / (z - 2*pi*j/L)."""
    total = 0.0
    for j in range(-N, N+1):
        pole = 2.0 * np.pi * j / L
        if abs(z - pole) < 1e-14: return np.inf
        total += xi[j+N] / (z - pole)
    return total

def find_zeros(xi, N, L, search_range):
    """Nullstellen von F(z) via brentq zwischen Polen."""
    from scipy.optimize import brentq
    poles = [2.0 * np.pi * j / L for j in range(-N, N+1)]
    z_min, z_max = search_range
    z_test = np.linspace(z_min, z_max, 30000)
    mask = np.ones_like(z_test, dtype=bool)
    for pole in poles:
        mask &= np.abs(z_test - pole) > 5e-4
    z_test = z_test[mask]
    F_vals = np.array([F_z(z, xi, N, L) for z in z_test])
    zeros = []
    for i in range(len(z_test)-1):
        if F_vals[i] * F_vals[i+1] < 0:
            pole_between = any(z_test[i] < p < z_test[i+1] for p in poles)
            if not pole_between:
                try:
                    root = brentq(lambda z: F_z(z, xi, N, L), z_test[i], z_test[i+1], xtol=1e-12)
                    zeros.append(root)
                except:
                    pass
    return sorted(zeros)

# =================================================================
# 7. Main
# =================================================================

def main():
    lam = LAMBDA
    N = N_DEFAULT
    L = 2 * mp.log(lam)  # L = log(lambda^2) = log(14)
    L_float = float(L)

    print("="*80)
    print("dirichlet_ccm_mpmath.py")
    print(f"mpmath Praezision: {mp.mp.dps} Dezimalstellen")
    print(f"lambda = sqrt(14) = {float(lam):.10f}")
    print(f"L = 2*log(lambda) = {L_float:.10f}")
    print(f"N = {N}, 2*pi/L = {2*np.pi/L_float:.6f}")
    print("="*80)

    # TEST 1: Riemann mit SIGN_WR=+1 (MS1-konform)
    sign_wr = SIGN_WR_DEFAULT
    print(f"\n{'#'*80}")
    print(f"RIEMANN (chi_0) Sanity-Check -- SIGN_WR = {sign_wr:+d}")
    print(f"{'#'*80}")
    T = build_QW(N, L, lam, chi_trivial, 1, +1, include_W02=True, cache_tag='riemann', sign_WR=sign_wr)
    run_riemann_diag(T, N, L_float, sign_wr)
    return

def run_riemann_diag(T, N, L_float, sign_wr):
    print(f"QW-Matrix-Norm: {np.linalg.norm(T):.4f}, diag-Werte: {np.diag(T)[N-2:N+3]}")

    # Diagnostics: Diagonal-Struktur und alle EWs
    diag_T = np.diag(T)
    print(f"Diagonale QW: min={diag_T.min():.6f}, max={diag_T.max():.6f}, spread={diag_T.max()-diag_T.min():.6e}")
    print(f"  Abweichung diag von mean: {np.std(diag_T):.4e}")
    # Off-Diagonal-Skala
    off_diag = T - np.diag(np.diag(T))
    print(f"  Off-diag: max|.|={np.max(np.abs(off_diag)):.6e}, Frobenius={np.linalg.norm(off_diag):.4f}")

    from scipy.linalg import eigh
    w_all, V_all = eigh(T)
    # delta_N-Overlap pro EV
    delta_N = np.ones(2*N+1) / np.sqrt(2*N+1)
    overlaps = V_all.T @ delta_N  # (M,)
    # Sortiere nach |overlap| descending
    ord_by_overlap = np.argsort(-np.abs(overlaps))
    print(f"\n  Top-10 Eigenvektoren nach |<delta_N|v>|:")
    for rank in range(10):
        k = ord_by_overlap[rank]
        v = V_all[:, k]
        is_even = np.allclose(v, v[::-1], atol=1e-6)
        print(f"    k={k:3d}  w={w_all[k]:+.8f}  |overlap|={abs(overlaps[k]):.4f}  even={is_even}")

    eps, xi, k_gs = find_ground_state(T, 'even')
    print(f"\nGrundzustand (v2, Projektion): k={k_gs}, eps_N={eps:.6f}, parity even: {np.allclose(xi, xi[::-1], atol=1e-6)}")
    # Normalisierung: delta_N(xi) = L^{-1/2} * sum(xi) = 1 => sum(xi) = sqrt(L)
    xi_sum = np.sum(xi)
    print(f"sum(xi) vor Norm = {xi_sum:.6f} (Ziel sqrt(L)={np.sqrt(L_float):.6f})")
    if abs(xi_sum) > 1e-12:
        xi = xi / xi_sum * np.sqrt(L_float)
    print(f"sum(xi) nach Norm = {np.sum(xi):.6f}")

    # ALTERNATIVE: nimm den EV mit maximalem |<delta_N|v>| (even)
    print(f"\n  [ALT] Teste Eigenvektor mit max |<delta_N|v>| (even):")
    for rank in range(len(w_all)):
        k = ord_by_overlap[rank]
        v = V_all[:, k]
        if np.allclose(v, v[::-1], atol=1e-6):
            print(f"    -> k={k}, w={w_all[k]:.6f}, overlap={overlaps[k]:.4f}")
            xi_alt = v * np.sign(overlaps[k])  # positive overlap convention
            xi_alt = xi_alt / np.sum(xi_alt) * np.sqrt(L_float)
            zeros_alt = find_zeros(xi_alt, N, L_float, (0.5, 2.0 * np.pi * N / L_float))
            print(f"    [ALT] Erste 8 Nullstellen:")
            known_list = [float(z) for z in RIEMANN_ZEROS[:8]]
            for i, z in enumerate(zeros_alt[:8]):
                err = z - known_list[i] if i < len(known_list) else None
                err_str = f"Fehler: {err:+.3e}" if err is not None else ""
                print(f"      z_{i+1} = {z:.10f}  {err_str}")
            break

    z_max = 2.0 * np.pi * N / L_float
    zeros = find_zeros(xi, N, L_float, (0.5, z_max))
    print(f"\nNullstellen von F(z) in (0.5, {z_max:.2f}):")
    known = [float(z) for z in RIEMANN_ZEROS[:5]]
    for i, z in enumerate(zeros[:7]):
        if i < len(known):
            err = z - known[i]
            print(f"  z_{i+1} = {z:.10f}   (bekannt: {known[i]:.10f}, Fehler: {err:+.3e})")
        else:
            print(f"  z_{i+1} = {z:.10f}")

    print("\n[alpha.1-Test aktiv: ueberspringe chi_4/chi_33 bis Riemann sanity ok]")
    return

    # TEST 2: chi_4
    print(f"\n{'#'*80}")
    print("chi_4 (odd, q=4)")
    print(f"{'#'*80}")
    T = build_QW(N, L, lam, chi_4, 4, -1, include_W02=False)
    eps, xi, k_gs = find_ground_state(T, 'odd')
    print(f"Grundzustand: eps_N={eps:.6f}, parity odd: {np.allclose(xi, -xi[::-1], atol=1e-6)}")
    zeros = find_zeros(xi, N, L_float, (0.5, z_max))
    known = [float(z) for z in CHI4_ZEROS_STR[:5]]
    print(f"\nNullstellen von F(z):")
    for i, z in enumerate(zeros[:7]):
        if i < len(known):
            err = z - known[i]
            print(f"  z_{i+1} = {z:.10f}   (bekannt: {known[i]:.10f}, Fehler: {err:+.3e})")
        else:
            print(f"  z_{i+1} = {z:.10f}")

    # TEST 3: chi_33
    print(f"\n{'#'*80}")
    print("chi_33 (even, q=33)")
    print(f"{'#'*80}")
    T = build_QW(N, L, lam, chi_33, 33, +1, include_W02=False)
    eps, xi, k_gs = find_ground_state(T, 'even')
    print(f"Grundzustand: eps_N={eps:.6f}, parity even: {np.allclose(xi, xi[::-1], atol=1e-6)}")
    # Atlas-Nullstellen fuer chi_33
    chi33_zeros = [2.997, 4.459, 6.362, 7.803, 10.169]
    zeros = find_zeros(xi, N, L_float, (0.5, z_max))
    print(f"\nNullstellen von F(z):")
    for i, z in enumerate(zeros[:7]):
        if i < len(chi33_zeros):
            err = z - chi33_zeros[i]
            print(f"  z_{i+1} = {z:.10f}   (bekannt: {chi33_zeros[i]:.6f}, Fehler: {err:+.3e})")
        else:
            print(f"  z_{i+1} = {z:.10f}")

    print("\nFertig.")

if __name__ == "__main__":
    main()
