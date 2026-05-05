# Session 19 — CCM Fourier-Basis: Ergebnisse und Durchbruch

**Datum:** 2026-04-18
**Autor:** LG
**Methode:** CCM 2025 (arXiv:2511.22755) Fourier-Basis V_n, full mpmath (dps=100)

## DURCHBRUCH: Cluster-Rand-Eigenvektor

Der Eigenvektor am **Rand des degenerierten Grundzustandsclusters** (nicht der Grundzustand selbst!) reproduziert Riemann-ζ-Nullstellen auf **10^-14** Präzision. Das ist 10^6 × besser als die IEW-Methode (Inverse-Energy-Weighting).

### Mechanismus

Die Weil-Quadratik QW hat in der kontinuierlichen Formulierung einen Nullraum (Kern). Die Fourier-Trunkierung auf dim=2N+1 hebt diesen zu einem fast-degenerierten Eigenwertcluster mit Gaps von 10^-53 bis 10^-14. Der Eigenvektor am **Clusterrand** (größter Gap innerhalb des Clusters) enthält die meiste Nullstelleninformation.

### Riemann ζ(s) — λ=3, N=30

| k | Gap von w_min | |F(γ₁)| | |F(γ₂)| | |F(γ₃)| |
|---|---|---|---|---|
| 0 | 0 | 4.4 | 0.053 | 0.002 |
| 1 | 1.9e-30 | 0.026 | 7.3e-4 | 3.5e-5 |
| 2 | 1.8e-24 | 1.2e-4 | 9.0e-6 | 6.4e-7 |
| 3 | 4.5e-19 | 4.0e-7 | 1.0e-7 | 1.1e-8 |
| **4** | **2.8e-14** | **4.0e-8** | **7.9e-9** | **9.6e-9** |

**Scan-Ergebnisse für k=4:**

| Nullstelle | LMFDB-Wert | Scan-Wert | Fehler |
|---|---|---|---|
| γ₁ | 14.134725141734693 | 14.134725141734748 | **5.5e-14** |
| γ₂ | 21.022039638771554 | 21.022039638771581 | **2.8e-14** |
| γ₃ | 25.010857580145688 | 25.010857580145725 | **3.6e-14** |
| γ₄ | 30.424876125859513 | 30.424876125860024 | **5.1e-13** |
| γ₅ | 32.935061587739189 | 32.935061587740776 | **1.6e-12** |

**Vergleich IEW (gleiche Konfiguration):**

| Nullstelle | IEW-Fehler | Cluster-Rand-Fehler | Verbesserung |
|---|---|---|---|
| γ₁ | 8.6e-10 | 5.5e-14 | **15.600×** |
| γ₂ | 2.6e-07 | 2.8e-14 | **9.300.000×** |
| γ₃ | 6.7e-06 | 3.6e-14 | **186.000.000×** |
| γ₄ | 4.6e-07 | 5.1e-13 | **900.000×** |
| γ₅ | 2.3e-08 | 1.6e-12 | **14.000×** |

## Dirichlet L-Funktionen

### L(s, χ₄) — λ=10, N=30

| Nullstelle | Scan-Fehler (IEW) | Scan-Fehler (Cluster k=9) |
|---|---|---|
| γ₁ = 6.021 | 7.9e-7 | 5.0e-7 |
| γ₂ = 10.244 | 4.7e-6 | 3.7e-6 |
| γ₃ = 12.988 | 8.5e-6 | 1.9e-6 |
| γ₄ = 16.343 | Pole-Kollision | 6.9e-4 (Pole) |
| γ₅ = 18.292 | 3.0e-6 | 3.0e-6 |

### L(s, χ₃) — λ=10, N=30

| Nullstelle | Scan-Fehler (IEW) | Scan-Fehler (Cluster k=11) |
|---|---|---|
| γ₁ = 8.040 | 1.2e-8 | **3.9e-9** |
| γ₂ = 13.923 | nicht gefunden | 2.4e-2 |

## Schlüsselerkenntnis: W02 erzeugt die Degeneracy

| Konfiguration | λ | W02? | Cluster-Größe | Best err(γ₁) |
|---|---|---|---|---|
| Riemann | 3.0 | Ja | 5 | **5.5e-14** |
| χ₄ | 3.0 | Nein | 1 (keine Degeneracy) | 6.9e-7 |
| χ₄ | 10.0 | Nein | 10 | 5.0e-7 |
| χ₃ | 10.0 | Nein | 12 | 3.9e-9 |

Die W02-Matrix (Pol bei s=1 der Riemann-ζ-Funktion) erzeugt bei kleinem λ bereits einen degenerierten Cluster. Ohne W02 (Dirichlet-Charaktere) entsteht der Cluster erst bei großem λ, wo aber Pole-Nullstellen-Kollisionen die Präzision limitieren.

## N-Skalierung (N=30 vs N=60)

| Nullstelle | N=30 err | N=60 err | Faktor |
|---|---|---|---|
| γ₁ | 5.5e-14 | 2.3e-14 | 2.4× |
| γ₂ | 2.8e-14 | 1.4e-14 | 2.0× |
| γ₃ | 3.6e-14 | 1.4e-14 | 2.5× |
| γ₄ | 5.1e-13 | 2.1e-13 | 2.4× |
| γ₅ | 1.6e-12 | 6.5e-13 | 2.5× |

- **Konvergenz:** Algebraisch (~N^{-0.7}), nicht exponentiell
- **Clustergröße:** Stabil bei 5 (unabhängig von N) — bestimmt durch W02-Struktur
- **IEW bei N=60 schlechter!** 9.6e-9 vs 8.6e-10 (N=30) — destruktive Interferenz
- **Float64-Limit:** Scan-Präzision bei 10⁻¹⁴ nahe am float64-Limit (~10⁻¹⁶)

## Offene Fragen

1. **Prolate-Basis:** Würde die CCM-eigene Prolate-Basis die Dirichlet-Limitierungen lösen?
2. **W02-Analogon für χ:** Gibt es eine Modifikation der Konstruktion, die für nicht-triviale Charaktere Degeneracy erzeugt?
3. **mpmath-Bisection:** Scan in voller mpmath-Präzision für Werte jenseits 10⁻¹⁶

## Dateien

- `_scripts/dirichlet_ccm_fourier_mp.py` — Hauptbibliothek (681 Zeilen)
- `_scripts/multi_eigenvector_analysis.py` — Multi-EV-Analyse
- `_scripts/cluster_boundary_test.py` — Cluster-Rand-Test alle L-Funktionen
- `_scripts/optimal_lambda_benchmark.py` — λ-Optimierung
- `_scripts/lambda_scan_plateau.py` — λ-Scan (Plateau-Diagnose)
- `_results/batch_beta/multi_eigenvec_riemann_lam3.log` — Durchbruch-Log
- `_results/batch_beta/cluster_boundary_ALL.log` — Vergleich alle L-Funktionen
