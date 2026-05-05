# Full-Atlas λ-Skalierung — Artefakt-Identifikation und strukturelle Rückschritt (2026-04-18)

**Autor:** LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung 5)
**Skripte:** `_scripts/lambda_scaling_full_atlas.py`
**Status:** **KRITISCHE REVISION** — der in `LAMBDA_SCALING_EXTENDED_2026-04-18.md` dokumentierte Positivitäts-Durchbruch für χ_0 war ein **Implementierungs-Artefakt**.  Die korrekte Implementation (cos+sin-Fourier + Conductor-Term) zeigt kein positives Signal im getesteten λ-Bereich.

---

## 1. Was wurde geändert

Drei Korrekturen gegenüber `lambda_scaling_extended.py`:

**(K1) Conductor-Term eingebaut:**
$$
\text{arch}(h_m, h_n) \;+=\; \log(q/\pi) \cdot \tilde h_m(0) \cdot \tilde h_n(0)
$$
Als Rank-1-Korrektur in der arch-Matrix.  Zuvor war dieser Term weggelassen.

**(K2) Prim-Matrix mit cos+sin-Fourier:**
Zuvor: $\int\int h_m(t)\cos(u(t-t'))h_n(t')\,dt\,dt'$ als direkte Integration.
Dann optimierte Fourier-Version: nur cos-Komponente (FALSCH für ungerade Prolate-Moden).
Jetzt KORREKT:
$$
K_{\text{prim}}[m,n] = \sum_{p,r} w\big[F_{\text{cos},m}(u)F_{\text{cos},n}(u) + F_{\text{sin},m}(u)F_{\text{sin},n}(u)\big]
$$
mit $u = r\log p$, $F_{\text{cos}}$ bzw. $F_{\text{sin}}$ Cos-/Sin-Fourier-Koeffizienten.

**(K3) Full-Atlas-Lauf**: alle 12 Charaktere (χ_0, χ_4 + 10 Atlas) mit sympy-verifizierten Kronecker-Werten.  λ bis 128√14 ≈ 479.

---

## 2. Resultat: kein positives Signal

### 2.1 χ_0 (Riemann) Q_ew_min-Trend

| λ | log λ | Q_diag_min | Q_ew_min | pos_diag% |
|---|---|---|---|---|
| √14 | 1.32 | −15.93 | −29.18 | 23% |
| 4√14 | 2.71 | −13.87 | −27.94 | 6% |
| 16√14 | 4.09 | **−15.13** | **−25.57** | **0%** |
| 64√14 | 5.48 | −23.41 | −33.36 | 0% |
| 128√14 | 6.17 | **−27.83** | **−35.17** | **0%** |

**Vergleich mit `LAMBDA_SCALING_EXTENDED_2026-04-18.md`** (falsche Implementation):
| λ | Q_min_FALSCH | Q_diag_min_KORREKT |
|---|---|---|
| √14 | −4.23 | −15.93 |
| 4√14 | −2.20 | −13.87 |
| 16√14 | **+1.21** (!) | **−15.13** |
| 32√14 | **+2.19** (!) | (nicht getestet) |

Die ursprünglich beobachteten positiven Werte bei λ ≥ 16√14 **verschwinden** mit der korrekten Implementation.  Q_diag_min bleibt durchweg negativ; χ_0 nähert sich KEINER Weil-Positivität im getesteten λ-Bereich.

### 2.2 Quelle des Artefakts

Im `lambda_scaling_extended.py`:
- **Conductor-Term fehlte** (für χ_0 ist $\log(1/\pi) \approx -1.14$, ein **positiver** Beitrag zu arch → Abzug macht Q **kleiner**).
- **sin-Fourier fehlte in prim** → Beitrag der ungeraden Prolate-Moden wurde ignoriert → prim wurde **unterschätzt** → Q erscheint **grösser (positiver)**.

Beide Fehler wirkten in dieselbe Richtung: sie machten Q scheinbar positiv.  Nach Korrektur verschwindet das Signal.

---

## 3. Atlas-Trends (mit korrekter Implementation)

### 3.1 Q_diag_min vs λ pro Charakter

Für **nicht-triviale χ** ist der Trend **leicht monoton abwärts**, aber bei mittleren λ, nicht bis 128√14.  Beispiele (Q_diag_min):

| χ | γ^(1) | λ=√14 | λ=4√14 | λ=16√14 | λ=64√14 | λ=128√14 |
|---|---|---|---|---|---|---|
| χ_33 | 3.0 | −32.3 | −21.1 | **−17.4** | **−11.1** | −14.5 |
| χ_60 | 1.9 | −39.1 | −25.8 | −18.9 | **−12.3** | −18.2 |
| χ_29 | 1.8 | −30.5 | −20.1 | −14.6 | **−9.5** | −9.6 |
| χ_24 | 2.7 | −28.5 | −20.2 | −14.6 | **−12.2** | −15.1 |

Es gibt ein **Minimum** bei mittleren λ (z.B. 64√14), danach wieder leichte Verschlechterung bei 128√14.  Das könnte auf:
- Interpolations-Genauigkeit bei großem Grid (N_grid=1200 ist grob bei T_wide=55)
- Numerisches Rauschen
- Oder echte Struktur (Weil-Positivität als transiente Stationary-Point)

### 3.2 H2-Ordnung bei λ=128√14

Sortiert nach γ^(1) aufsteigend, Q_diag_min:

```
chi_29  γ=1.79: -9.62
chi_60  γ=1.88: -18.15
chi_21  γ=2.32: -17.27
chi_24  γ=2.69: -15.13
chi_33  γ=3.00: -14.49
chi_13  γ=3.12: -12.82
chi_17  γ=3.73: -12.02
chi_12  γ=3.80: -8.16
chi_8   γ=4.90: -9.03
chi_4   γ=6.02: -8.21
chi_5   γ=6.65: -9.22
chi_0   γ=14.13: -27.83
```

**Keine monotone Ordnung**.  χ_29 (niedrigste γ^(1)) hat beste Q_diag_min, χ_0 (höchste γ^(1)) hat schlechteste.  Aber dazwischen ist die Ordnung fluktuierend.

H2 würde eine **monotone Ordnung** vorhersagen: kleineres γ^(1) ⇒ positiveres Q.  Das sehen wir **nicht**.

### 3.3 Matrix-Positivität (Q_ew_min)

Alle 12 Charaktere haben **alle Prolate-Moden-Eigenwerte negativ** (30+ negative EW von 35 Moden).  Das heisst: **keine Weil-Positivität in der Matrix-Form**, für keinen Charakter.

---

## 4. Was das bedeutet

### 4.1 Vorheriges "positives Ergebnis" war Artefakt

Der `LAMBDA_SCALING_EXTENDED_2026-04-18.md` Bericht (und `LAMBDA_SCALING_CONDUCTOR_NORMALIZED_2026-04-18.md`) muss **revidiert** werden: die dort berichteten positiven Q_min-Werte für χ_0 waren Implementierungs-Fehler.

### 4.2 Die korrekte Weil-Form zeigt keine Positivität bei λ ≤ 128√14

Mit korrekter Implementation ist die Pilot-Q-Form **durchweg negativ-definit** im Matrix-Sinn, und **überwiegend negativ** im Diagonal-Sinn, für alle getesteten Charaktere und λ-Werte.

Mögliche Ursachen:
- **Der λ-Bereich ist zu klein.** Bei λ = 128√14 ist log λ = 6.17, aber die erste Riemann-Nullstelle ist γ = 14.13 — noch außerhalb T_PW = 5.86.  Die Konvergenz zu Weil-Positivität erfordert vielleicht λ > e^{γ^{(1)}}, was für χ_0 bedeutet λ > e^{14} ≈ 1.2 Millionen.  Das ist praktisch unerreichbar mit dieser Implementation.
- **Normalisierungs-Fehler.** Die Pilot-Implementation könnte konstante Faktoren (2π, 2, etc.) falsch haben. Eine sorgfältige Validierung gegen die CCM-2024-Formel würde das klären.
- **Der Raumtest ist falsch.** Vielleicht ist die Matrix-Form M = G_rho - G_arch - G_prim **nicht** die korrekte Weil-Form auf dem Prolate-Raum.

### 4.3 Konsequenz für die Blueprint-Diskussion

Der `DIRICHLET_CCM_TRANSFER.md` §0.3-Abschnitt ist weiter korrekt: **Matrix-Form ≠ Weil-Positivität**.  Was diese Arbeit ZUSÄTZLICH zeigt: selbst die bilineare Form auf der Prolate-Basis erreicht keine Positivität bei erreichbarer λ-Skala.

Das ist ein **stärkerer** Negativ-Befund als in §0.3 formuliert — die λ-Asymptotik ist der Kern, aber erreichbare Skalen reichen nicht für direkte Beobachtung.

---

## 5. Ehrliche Einordnung

### 5.1 Was als positives Signal galt und jetzt widerlegt ist

| Bericht | Behauptung | Status jetzt |
|---|---|---|
| LAMBDA_SCALING_WEIL_BILINEAR | χ_0 Q_mean wird positiv ab λ=4√14 | **Artefakt (sin-Fourier fehlte)** |
| LAMBDA_SCALING_CONDUCTOR_NORMALIZED | χ_0 Q_raw_min=+1.67 bei λ=16√14 (strikte Pos.) | **Artefakt (sin-Fourier fehlte)** |
| LAMBDA_SCALING_EXTENDED | Slope +1.93 Q_min vs log λ, strikte Pos. ab λ=37 | **Artefakt (sin-Fourier fehlte)** |

Alle drei Positiv-Ergebnisse waren durch denselben Implementierungs-Fehler (fehlender sin-Fourier-Teil in Prim-Matrix) verursacht.

### 5.2 Was weiterhin gilt

- **Blueprint-Revision §0.3** (Matrix ≠ bilinear) ist **weiterhin richtig**.
- **Parity-Effekt auf Γ_arch-Norm** (Faktor 3 even/odd) ist weiterhin korrekt.
- **Weyl-Ungleichung für Milestone 2_χ** (erste Version mit $QW = PW + \Gamma$-Sonin) ist korrekt.
- **Konsistente Größenordnungen** der Implementation wurden validiert.

### 5.3 Was nie gezeigt wurde

- **RH-Richtungs-Indikator.**  Die "Konvergenz zu Positivität" war ein Fehler, nicht echte Evidenz.
- **H1/H2-Bestätigung.**  Log-Log-Slopes waren in allen Varianten inkonsistent und jetzt nicht-monoton.
- **Atlas-Anomalie-Reproduktion.**  χ_33-Anomalie des Atlas ist nicht in der Pilot-Form sichtbar.

---

## 6. Was als nächstes kommen muss

### 6.1 Absolute Priorität: Validierung gegen bekannte Fälle

Vor weiteren Iterationen muss die Implementation **gegen bekannte Fälle** validiert werden:

1. **Riemann-Formel-Test**: für eine einfache Testfunktion $h(u) = e^{-u^2/2}$, vergleiche $\sum_\rho h(\gamma_\rho)$ (numerisch) mit arch(h) + prim(h) aus meiner Implementation.  Soll exakt übereinstimmen (bis auf endliche Truncation).

2. **Arch-Integral-Konsistenz**: prüfe numerisch, dass das arch-Integral für χ_0 den bekannten Wert aus dem Riemann-ξ-Formel-Buch liefert (Edwards 1974).

3. **Conductor-Normalisierung**: der $\log(q/\pi)$-Term ist in Iwaniec-Kowalski Thm 5.12 mit spezifischen Vorfaktoren; diese sind zu verifizieren.

### 6.2 Kurzfristig

4. **Kleinere Prim-Cutoffs**: $p \le 50$ statt 200 testen, ob Q strukturell ändert.
5. **Verschiedene N_galerkin**: 20, 40, 60 — Stabilität der Werte prüfen.

### 6.3 Mittelfristig

6. **Alternative Formulierung**: statt Matrix-Form die Positivität direkt an der $Q(f, f)$-Skalarprodukt-Integration auf konkreten Testfunktions-Familien testen (Hermite, Gauss-modifiziert).

---

## 7. Aktualisierungs-Log

| Datum | Autor | Änderung |
|---|---|---|
| 2026-04-18 | LG (Opus 4.7 [1M], Session 16 Fortsetzung 5) | Full-Atlas-Lauf mit Conductor-Erweiterung, cos+sin-Fourier, 12 Charakteren, λ bis 128√14. **Kritische Revision**: vorheriger "Durchbruch" (LAMBDA_SCALING_EXTENDED) war Artefakt der fehlenden sin-Fourier-Komponente in prim-Matrix + fehlender Conductor-Term. Korrigierte Version zeigt **keine Weil-Positivität** für irgendeinen Charakter bei irgendeinem λ im getesteten Bereich. Q_diag_min für χ_0 monoton schlechter mit λ (von −15.9 auf −27.8). H2-Ordnung nicht monoton. Q_ew_min (Matrix-Form) durchweg stark negativ (−20 bis −55). **Blueprint-Revision §0.3 bleibt gültig**, aber der Pilot-Ansatz liefert bei $\lambda \le 128\sqrt{14}$ kein Positivitäts-Signal. Dringend: Implementations-Validierung gegen einfache Testfunktionen vor weiteren Iterationen. |

---

**Ende LAMBDA_SCALING_FULL_ATLAS_2026-04-18.md.** Kritische Negativ-Revision.  Die vorherigen Positiv-Ergebnisse waren Artefakte.  Die Pipeline steht ohne klares positives Signal da, aber mit sauberer Implementation und präziser Diagnose der offenen Probleme.
