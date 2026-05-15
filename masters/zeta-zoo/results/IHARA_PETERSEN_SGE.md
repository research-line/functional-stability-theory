# Ihara-Zeta SGE YES-side Test: Petersen Graph

**Datum:** 2026-04-16
**Skript:** `_scripts/ihara_petersen_sge_test.py`
**Motivation:** SGE-Hypothese auf der YES-Seite testen. Petersen ist 3-regulär, 10 Knoten, Ramanujan-Graph, Aut(G) = S_5. Vorhersage: HP-BL(zeta_Petersen) = YES.

## Graph-Parameter

- n = 10 Knoten, m = 15 ungerichtete Kanten
- Alle Grade = 3 (3-regulär)
- Aut(G) = 120 Elemente (S_5, Ordnung 120)
- Knoten-Orbits: 1, gerichtete-Kanten-Orbits: 1

## Adjazenz-Spektrum

Eigenvalues von A: ['+3.000', '+1.000', '+1.000', '+1.000', '+1.000', '+1.000', '-2.000', '-2.000', '-2.000', '-2.000']

Erwartung (Petersen): {3, 1(5×), -2(4×)}. Ramanujan-Bound 2√(q-1) = 2√2 = 2.8284. Max |nichttriviale| = 2.0000.

## Bass-Hashimoto-Identität

Bei $u = (0.3+0.1j)$:
- $\det(I - uB) = (1.011447122913742-0.09688626006229914j)$
- $(1-u^2)^{m-n} \det(I - Au + Qu^2) = (1.0114471229137416-0.096886260062299j)$
- Verhältnis = (1.0000000000000004-9.517180929718275e-17j)
- **Stimmt überein: True**

## Ramanujan-RH für Ihara

- Nichttriviale Ihara-Nullstellen auf kritischem Kreis $|u| = 1/\sqrt{2} = 0.7071$: 18/19
- |u|-Bereich: [0.5000, 0.7071]
- **Klassisches Ramanujan-Resultat bestätigt** (RH-Analog für Ihara hält)

## HP-Operator-Kandidat: Graph-Laplacian Δ = D - A

- Symmetrisch (selbstadjungiert): ✓
- Eigenwerte: ['-0.000', '2.000', '2.000', '2.000', '2.000', '2.000', '5.000', '5.000', '5.000', '5.000']
- Über $2u^2 - \lambda_A u + 1 = 0$ (für 3-regulär) mit Ihara-Nullstellen verknüpft

## SGE-YES-Zertifikat

- Automorphismus σ = [0, 1, 2, 7, 5, 4, 6, 3, 9, 8]
- Permutationsmatrix P kommutiert mit A: **True**
- P kommutiert mit Δ: **True**
- Kanten-Lift P_edge kommutiert mit B: **True**
- **Nicht-skalarer kommutierender Operator existiert: True**

## Interpretation

**Ergebnis: Teilweise Bestätigung.** Einige Checks bestanden, andere nicht. Detail-Debug erforderlich.
