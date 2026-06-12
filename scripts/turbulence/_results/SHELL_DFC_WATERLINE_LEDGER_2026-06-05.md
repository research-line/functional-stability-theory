# Shell-DFC-Waterline-Ledger 2026-06-05

Abschlusszeit: `2026-06-05T23:59:59+02:00`

## Zweck

Dieser Lauf erweitert den bisherigen Toy-`DFC1^vee`-Ledger um einen
kurzen dynamischen Sabra-Shell-Smoke und Matched Controls. Er ist kein
DNS-Nachweis, kein Satz über Navier--Stokes und kein Upload-Gate.

## Externer Paperstand

- Benavides--Bustamante 2026 (`arXiv:2507.03397v2`) stützt Shell-Modelle
  als Flux-Testbett, betont aber Phasendynamik als eigene Flussdeterminante.
- Tuteri--Chibbaro--Alexakis 2026 (`arXiv:2603.11892`) zeigt, dass
  Kaskadenrichtung und Shell-Geometrie bei 2D/dual-cascade-Modellen nicht
  blind aus klassischen Shell-Modellen übertragen werden dürfen.
- Zinchenko--Schumacher 2026 (`arXiv:2508.03401`) erweitert Duchon--Robert
  auf kompressible Flüsse; das ist Kontext, aber kein direkter Eingang für
  die inkompressible Paper-B-Bridge.
- JHTDB Forced Isotropic Turbulence bleibt das nächste sinnvolle DNS-Ziel:
  `1024^3`, `R_lambda=418`, `epsilon=0.103` und gespeicherte Velocity-/Pressure-Frames.

## Sabra-Smoke-Parameter

| Parameter | Wert |
|---|---:|
| Shells | 18 |
| Snapshots | 4000 |
| `dt` | 0.0002 |
| Thermalisierung | 3 |
| Datenzeit | 8 |
| `nu` | 1e-06 |
| Forcing-Schale | 2 |
| Forcing-Amplitude | 0.02 |
| mittlere Energie | 0.021757 |
| mittlere Dissipation | 2.75542e-09 |
| Laufzeit s | 7.503 |

## Ergebnis-Tabelle

| Szenario | Fenster | DFC2 | DFC1^vee | weighted/target | residual/allowance | bad corridor | Urteil |
|---|---:|---:|---:|---:|---:|---:|---|
| sabra_transition_window_4_12 | 4-12 | fail | blocked | n/a | n/a | n/a | blocked: nonmonotone Free-Energy weights |
| sabra_tail_window_actual_order_8_14 | 8-14 | pass | fail | 0.21779 | 5.21473 | 0.898632 | fail: Abel-weighted waterline deficit |
| sabra_tail_window_sorted_high_weight_control | 8-14 | pass | pass | 4.80991 | 0 | 0.0127959 | pass control: same flux multiset only after artificial placement |
| sabra_tail_window_sorted_low_weight_control | 8-14 | pass | fail | 0.115 | 5.9 | 0.918631 | fail: Abel-weighted waterline deficit |

## Befund

- Das feste Übergangsfenster `4-12` wird blockiert, weil die
  Free-Energy-Gewichte nicht monoton sind. Das verhindert einen
  nachträglichen Kaskadenfront-Claim.
- Im Tail-Waterline-Fenster `8-14` sind die Gewichte monoton, aber die
  tatsächliche Flux-Platzierung besteht `DFC1^vee` nicht
  (`weighted/target=0.217790`, `residual/allowance=5.21473`).
- Die sortierte High-Weight-Kontrolle besteht mit demselben Flux-Multiset
  nur nach künstlicher Platzierung. Die sortierte Low-Weight-Kontrolle
  scheitert noch stärker. Damit ist nicht der Mittelwert entscheidend,
  sondern die feste Abel-gewichtete Platzierung des Flux.

## Konsequenz für den Beweisstand

`DFC1^vee` bleibt offen. Der neue Lauf stärkt die Guardrails: Shell- oder
DNS-Daten zählen erst als Evidenz, wenn Fenster, Waterline, Phase-/Flux-
Provenienz, Matched Controls und Projektionsrest vorab festgelegt sind.
Ein v1.7-Upload sollte daraus keinen stärkeren Beweisclaim ableiten.

## Maschinenlesbare Begleitdaten

- `SHELL_DFC_WATERLINE_LEDGER_2026-06-05.json`
- `SHELL_DFC_WATERLINE_LEDGER_2026-06-05.csv`
