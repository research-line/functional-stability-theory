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
4. Diagnostiziere gute Skalen, Bad-Defekte und Saturationszonen, ohne aus
   einer endlichen Folge einen Kingman-, Gap- oder Kontinuumsclaim abzuleiten.

Birkhoff-Kontraktion: tau_B(R) = tanh(Delta(R)/4)
wobei Delta(R) = log(max_{i,j,k,l} R_{ik}R_{jl}/(R_{il}R_{jk}))
(Birkhoff 1957, Bushell 1973)

Für den Zusammenhang zwischen Birkhoff-Koeffizient und Spektralverhältnis
positiver Matrizen siehe Han--Han, arXiv:1906.04875. Zur Trennung von
projektiver Nicht-Expansion und strikter Kontraktion siehe Ligonnière,
arXiv:2312.11147.

Autor: Lukas Geiger (Skript erstellt per Claude, 2026)
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

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
                for ell in range(n):
                    if (
                        T_pos[i, k] > 0
                        and T_pos[j, ell] > 0
                        and T_pos[i, ell] > 0
                        and T_pos[j, k] > 0
                    ):
                        ratio = (T_pos[i, k] * T_pos[j, ell]) / (T_pos[i, ell] * T_pos[j, k])
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
# Good-Scale-Diagnostik
# ==========================================================================

@dataclass(frozen=True)
class SaturationZone:
    """Maximales zusammenhängendes Intervall nicht-guter RG-Skalen."""

    start_scale: int
    end_scale: int
    length: int
    defect_sum: float
    max_tau: float
    min_gap_to_one: float


@dataclass(frozen=True)
class GoodScaleDiagnostics:
    """Endliche, claim-neutrale Diagnose einer Folge von Birkhoff-Koeffizienten."""

    epsilon: float
    threshold: float
    n_scales: int
    good_count: int
    good_density: float
    bad_count: int
    bad_defect_sum: float
    max_bad_run: int
    good_mask: tuple[bool, ...]
    bad_defects: tuple[float, ...]
    prefix_good_density: tuple[float, ...]
    cumulative_bad_defect: tuple[float, ...]
    saturation_zones: tuple[SaturationZone, ...]


def analyze_good_scales(taus: Sequence[float], epsilon: float) -> GoodScaleDiagnostics:
    """
    Diagnostiziere gute und nahezu saturierte Skalen einer endlichen RG-Folge.

    Eine Skala ist gut, wenn ``tau_k <= 1 - epsilon``. Der Bad-Defekt ist
    ``d_k = max(0, tau_k - (1 - epsilon))``. Eine Saturationszone ist ein
    maximales zusammenhängendes Intervall mit ``d_k > 0``.

    Die Rückgabe beschreibt nur die eingegebene endliche Folge. Insbesondere
    beweisen weder positive Good-Scale-Dichte noch endliche Defektsumme eine
    uniforme RG-Kontraktion oder einen Kontinuums-Massengap.
    """
    if not 0.0 < epsilon < 1.0:
        raise ValueError("epsilon muss strikt zwischen 0 und 1 liegen")

    tau_values = np.asarray(list(taus), dtype=float)
    if tau_values.ndim != 1 or tau_values.size == 0:
        raise ValueError("taus muss eine nichtleere eindimensionale Folge sein")
    if not np.all(np.isfinite(tau_values)):
        raise ValueError("taus darf keine NaN- oder unendlichen Werte enthalten")
    if np.any((tau_values < 0.0) | (tau_values > 1.0)):
        raise ValueError("Birkhoff-Koeffizienten müssen im Intervall [0, 1] liegen")

    threshold = 1.0 - epsilon
    good_mask_array = tau_values <= threshold
    bad_defects_array = np.maximum(tau_values - threshold, 0.0)
    bad_mask_array = bad_defects_array > 0.0

    prefix_counts = np.cumsum(good_mask_array, dtype=int)
    prefix_lengths = np.arange(1, tau_values.size + 1, dtype=float)
    prefix_good_density = prefix_counts / prefix_lengths
    cumulative_bad_defect = np.cumsum(bad_defects_array)

    zones = []
    zone_start = None
    for scale, is_bad in enumerate(bad_mask_array):
        if is_bad and zone_start is None:
            zone_start = scale
        if zone_start is not None and (not is_bad or scale == tau_values.size - 1):
            zone_end = scale if is_bad else scale - 1
            zone_taus = tau_values[zone_start : zone_end + 1]
            zone_defects = bad_defects_array[zone_start : zone_end + 1]
            zones.append(
                SaturationZone(
                    start_scale=zone_start,
                    end_scale=zone_end,
                    length=zone_end - zone_start + 1,
                    defect_sum=float(np.sum(zone_defects)),
                    max_tau=float(np.max(zone_taus)),
                    min_gap_to_one=float(np.min(1.0 - zone_taus)),
                )
            )
            zone_start = None

    good_count = int(np.sum(good_mask_array))
    bad_count = int(np.sum(bad_mask_array))
    max_bad_run = max((zone.length for zone in zones), default=0)

    return GoodScaleDiagnostics(
        epsilon=float(epsilon),
        threshold=float(threshold),
        n_scales=int(tau_values.size),
        good_count=good_count,
        good_density=good_count / int(tau_values.size),
        bad_count=bad_count,
        bad_defect_sum=float(np.sum(bad_defects_array)),
        max_bad_run=max_bad_run,
        good_mask=tuple(bool(value) for value in good_mask_array),
        bad_defects=tuple(float(value) for value in bad_defects_array),
        prefix_good_density=tuple(float(value) for value in prefix_good_density),
        cumulative_bad_defect=tuple(float(value) for value in cumulative_bad_defect),
        saturation_zones=tuple(zones),
    )


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
# Bericht und CLI
# ==========================================================================

DEFAULT_BETAS = (1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Endliche Birkhoff-/Good-Scale-Diagnostik für die hierarchische RG-Kaskade."
    )
    parser.add_argument("--epsilon", type=float, default=0.01, help="Good-Scale-Marge epsilon (Default: 0.01)")
    parser.add_argument("--n-levels", type=int, default=6, help="Anzahl der RG-Stufen (Default: 6)")
    parser.add_argument("--n-bins", type=int, default=16, help="Transfermatrix-Diskretisierung (Default: 16)")
    parser.add_argument("--betas", type=float, nargs="+", default=list(DEFAULT_BETAS), help="Zu prüfende beta-Werte")
    parser.add_argument("--no-plot", action="store_true", help="Kein PNG erzeugen oder überschreiben")
    args = parser.parse_args(argv)
    if args.n_levels < 1:
        parser.error("--n-levels muss mindestens 1 sein")
    if args.n_bins < 2:
        parser.error("--n-bins muss mindestens 2 sein")
    if not args.betas:
        parser.error("--betas benötigt mindestens einen Wert")
    if not 0.0 < args.epsilon < 1.0:
        parser.error("--epsilon muss strikt zwischen 0 und 1 liegen")
    return args


def format_saturation_zones(diagnostics: GoodScaleDiagnostics) -> str:
    if not diagnostics.saturation_zones:
        return "keine"
    return ", ".join(
        f"k={zone.start_scale}..{zone.end_scale} "
        f"(L={zone.length}, D={zone.defect_sum:.6g}, max_tau={zone.max_tau:.6f})"
        for zone in diagnostics.saturation_zones
    )


def render_plot(all_results, diagnostics_by_beta, betas, epsilon):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(12, 14))

    ax = axes[0, 0]
    for beta in betas[:5]:
        res = all_results[beta]
        ks = [r["k"] for r in res]
        taus = [r["tau_B"] for r in res]
        ax.plot(ks, taus, "o-", lw=1.5, ms=5, label=f"$\\beta={beta}$")
    ax.axhline(y=1.0 - epsilon, color="darkorange", ls="--", lw=1.5, label="$1-\\epsilon$")
    ax.axhline(y=1.0, color="red", ls=":", lw=1.0, label="$\\tau_B=1$")
    ax.set_xlabel("RG-Stufe k")
    ax.set_ylabel(r"$\tau_B(R_k)$")
    ax.set_title("Birkhoff-Koeffizient und Good-Scale-Schwelle")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    finite_mean_logs = []
    for beta in betas:
        taus = [r["tau_B"] for r in all_results[beta]]
        finite_mean_logs.append(float(np.mean(np.log(np.clip(taus, 1e-300, 1.0)))))
    ax.plot(betas, finite_mean_logs, "rs-", lw=2, ms=8)
    ax.axhline(y=0, color="black", ls="-", lw=0.5)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$N^{-1}\sum_k \log \tau_B(R_k)$")
    ax.set_title("Endlicher Mittelwert (kein Kingman-Claim)")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    stride = max(1, len(betas) // 4)
    for beta in betas[::stride][:4]:
        res = all_results[beta]
        ks = [r["k"] for r in res]
        gaps = [r["gap"] for r in res]
        ax.plot(ks, gaps, "o-", lw=1.5, ms=5, label=f"$\\beta={beta}$")
    ax.set_xlabel("RG-Stufe k")
    ax.set_ylabel("Spektrallücke")
    ax.set_title("Transfermatrix-Spektrallücke")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    max_taus = [np.max([r["tau_B"] for r in all_results[beta]]) for beta in betas]
    ax.plot(betas, max_taus, "b^-", lw=2, ms=8)
    ax.axhline(y=1.0 - epsilon, color="darkorange", ls="--", lw=1.5)
    ax.axhline(y=1.0, color="red", ls=":", lw=1.0)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\max_k \tau_B(R_k)$")
    ax.set_title("Maximaler Birkhoff-Koeffizient")
    ax.grid(True, alpha=0.3)

    ax = axes[2, 0]
    good_densities = [diagnostics_by_beta[beta].good_density for beta in betas]
    ax.plot(betas, good_densities, "go-", lw=2, ms=7)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\rho_\epsilon$")
    ax.set_title(f"Good-Scale-Dichte ($\\epsilon={epsilon:g}$)")
    ax.grid(True, alpha=0.3)

    ax = axes[2, 1]
    defect_sums = [diagnostics_by_beta[beta].bad_defect_sum for beta in betas]
    max_runs = [diagnostics_by_beta[beta].max_bad_run for beta in betas]
    ax.bar(betas, defect_sums, width=0.6, alpha=0.65, color="tab:purple", label="$D_{bad}$")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel("Bad-Defektsumme", color="tab:purple")
    ax.tick_params(axis="y", labelcolor="tab:purple")
    ax.grid(True, axis="y", alpha=0.3)
    ax_run = ax.twinx()
    ax_run.plot(betas, max_runs, "ko--", lw=1.5, label="max. Zonenlänge")
    ax_run.set_ylabel("Maximale Saturationszonenlänge", color="black")
    ax.set_title("Bad-Defekte und Saturationszonen")

    plt.tight_layout()
    outpath = Path(__file__).with_name("compute_birkhoff_rg.png")
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return outpath


def main(argv=None):
    args = parse_args(argv)
    betas = list(args.betas)

    print("=" * 78)
    print("YANG-MILLS: endliche Birkhoff-/Good-Scale-Diagnostik für RG-Schritte")
    print(f"Gut: tau_B <= 1-epsilon mit epsilon={args.epsilon:g}; Diagnose, kein Gap-Beweis")
    print("=" * 78)

    all_results = {}
    diagnostics_by_beta = {}

    for beta in betas:
        print(f"\n{'=' * 58}")
        print(f"beta = {beta:.6g}")
        print(f"{'=' * 58}")

        res = rg_cascade(beta, n_levels=args.n_levels, n_bins=args.n_bins)
        all_results[beta] = res

        print(f"\n  {'k':>3} {'tau_B':>10} {'tau_spec':>10} {'Delta':>10} {'gap':>10} {'lambda_2/1':>12}")
        print(f"  {'-' * 60}")
        for result in res:
            ratio = result["lambda_2"] / result["lambda_1"] if result["lambda_1"] > 0 else 0.0
            print(
                f"  {result['k']:>3} {result['tau_B']:>10.6f} {result['tau_spectral']:>10.6f} "
                f"{result['Delta']:>10.4f} {result['gap']:>10.6f} {ratio:>12.6f}"
            )

        taus = [result["tau_B"] for result in res]
        finite_mean_log = float(np.mean(np.log(np.clip(taus, 1e-300, 1.0))))
        diagnostics = analyze_good_scales(taus, args.epsilon)
        diagnostics_by_beta[beta] = diagnostics

        print(f"\n  endlicher Mittelwert <log(tau_B)> = {finite_mean_log:.6f} (kein Kingman-Claim)")
        print(
            f"  good scales: {diagnostics.good_count}/{diagnostics.n_scales} "
            f"(rho_epsilon={diagnostics.good_density:.6f}, Schwelle={diagnostics.threshold:.6f})"
        )
        print(
            f"  bad defects: Summe={diagnostics.bad_defect_sum:.6g}, "
            f"bad scales={diagnostics.bad_count}, max run={diagnostics.max_bad_run}"
        )
        print(f"  Saturationszonen: {format_saturation_zones(diagnostics)}")

    print(f"\n{'=' * 94}")
    print("ZUSAMMENFASSUNG: endliche Good-Scale-Diagnostik über RG-Skalen")
    print(f"{'=' * 94}")
    print(
        f"\n  {'beta':>8} {'<tau_B>':>10} {'max tau':>10} {'<log tau>':>11} "
        f"{'rho_eps':>9} {'D_bad':>10} {'Zonen':>6} {'maxRun':>7} {'Status':>14}"
    )
    print(f"  {'-' * 91}")

    for beta in betas:
        res = all_results[beta]
        taus = [result["tau_B"] for result in res]
        diagnostics = diagnostics_by_beta[beta]
        finite_mean_log = float(np.mean(np.log(np.clip(taus, 1e-300, 1.0))))
        if diagnostics.bad_count == 0:
            status = "EPS-UNIFORM"
        elif diagnostics.good_count == 0:
            status = "SATURIERT"
        else:
            status = "GEMISCHT"

        print(
            f"  {beta:>8.3g} {np.mean(taus):>10.6f} {np.max(taus):>10.6f} {finite_mean_log:>11.5f} "
            f"{diagnostics.good_density:>9.5f} {diagnostics.bad_defect_sum:>10.5g} "
            f"{len(diagnostics.saturation_zones):>6} {diagnostics.max_bad_run:>7} {status:>14}"
        )

    print(
        f"""
  INTERPRETATION DER ENDLICHEN DIAGNOSE:
  - rho_epsilon ist der beobachtete Anteil der Skalen mit tau_B <= 1-epsilon.
  - D_bad = sum_k max(0, tau_B(R_k) - (1-epsilon)) misst die gesamte Überschreitung
    der Good-Scale-Schwelle; Saturationszonen sind maximale zusammenhängende Bad-Runs.
  - Weder <log(tau_B)> < 0 noch rho_epsilon > 0 oder eine kleine endliche D_bad-Summe
    beweisen Kingman-Voraussetzungen, skalenuniforme Koerzivität, OS-Kompaktheit oder
    eine Yang--Mills-Massenlücke. Dafür fehlen weiterhin eigenständige Transfersätze.
  - Diskretisierung: n_bins={args.n_bins}, n_levels={args.n_levels}; die Werte sind
    finite numerische Diagnostik des vereinfachten SU(2)-Transfermatrixmodells.
"""
    )

    if not args.no_plot:
        try:
            outpath = render_plot(all_results, diagnostics_by_beta, betas, args.epsilon)
            print(f"Plot: {outpath}")
        except Exception as exc:
            print(f"(Plot nicht erzeugt: {exc})")

    print("\n[DONE]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
