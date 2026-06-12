# Dual-DFC1-Ledger 2026-05-27

Abschlusszeit: `2026-05-27T10:31:31+02:00`

## Zweck

Dieser Kurzlauf operationalisiert die geplante Prüfung von `DFC1^vee`
für Turbulenz / Paper B. Er ist ein Ledger- und Negativkontrolltest,
kein DNS-Nachweis und kein neuer Satz über Navier-Stokes.

Getestet wird die Abel-duale Bedingung

```text
sum_j a_j Pi_j >= Pi0 * sum_j a_j - C_proj * (|phi_j1| + |phi_j2|)
a_j = phi_j - phi_{j+1},  phi_j = log(E_j/E_j*)
```

Parameter: `Pi0=1.0`, `C_proj_budget=0.15`,
Fenster: Shells `4-13`.

## Ergebnis-Tabelle

| Szenario | DFC2 | punktweise DFC1 | DFC1^vee | weighted/target | residual/allowance | bad corridor | Urteil |
|---|---:|---:|---:|---:|---:|---:|---|
| k41_reference | pass | pass | pass | 1 | n/a | 0 | vacuous reference control |
| smooth_forward | pass | pass | pass | 1.12811 | 0 | 0 | pass |
| alternating_tolerated | pass | fail | pass | 1.0204 | 0 | 0.292766 | dual pass despite pointwise backscatter |
| bad_corridor | pass | fail | fail | 0.635095 | 2.0956 | 0.746778 | fail: Abel-weighted backscatter corridor |
| same_mean_shuffle | fail | fail | fail | 0.492961 | 1.91343 | 0.748053 | fail: Abel-weighted backscatter corridor |

## Befund

- `smooth_forward` ist die Positivkontrolle: punktweise und duale DFC1 bestehen.
- `alternating_tolerated` zeigt den gewünschten Unterschied: einzelne negative Shell-Flüsse zerstören punktweise DFC1, aber nicht automatisch die Abel-duale Paarung.
- `bad_corridor` und `same_mean_shuffle` sind Negativkontrollen: gleicher oder hoher mittlerer Fluss reicht nicht, wenn Backscatter genau in den stark gewichteten Abel-Shells sitzt.
- Damit ist der nächste empirische Test nicht ein Mittelwert-Flux, sondern ein festes Fenster-/Waterline-Ledger mit `projection_residual`, `bad_corridor_capacity` und `same_mean_reference_risk`.

## Konsequenz für den Beweisstand

`DFC1^vee` bleibt offen und wird durch diesen Lauf nicht bewiesen. Der Fortschritt ist operativ: Die offene DR/Eyink -> LP/Wavelet-Projektionsbrücke ist jetzt in ein prüfbares Ledger zerlegt. Ein echter Nachweis braucht weiterhin Daten oder Analyse für den Projektionsfehler und eine quantitative DFC3-Schätzung.

## Maschinenlesbarer Begleitdatensatz

- `DUAL_DFC_LEDGER_2026-05-27.json`
