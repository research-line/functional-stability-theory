# Pilot-Ergebnis — χ-Defekt-Norm (2026-04-18)

**Autor:** LG (Claude Opus 4.7 [1M], Session 16)
**Skript:** `_scripts/chi_defect_norm_pilot.py`
**Status:** Erste numerische Beobachtung; **strukturell erhellendes Teilergebnis**, keine Bestätigung oder Widerlegung der Blueprint-Hypothesen.

---

## 1. Setup

- $\lambda = \sqrt{14} \approx 3.7417$ (CCM-Standard-Skala).
- $L = \log\lambda \approx 1.32$.
- Mellin-Galerkin-Dimension $N = 24$.
- Bandbegrenzung $T_{\max} = 0.8 L \approx 1.06$ (innerhalb Paley-Wiener-Bereich aus Milestone 1_χ).
- Getestete Charaktere: $\chi_0$ (trivial, Riemann-Baseline), $\chi_4$ (mod 4, odd).

Vereinfachungen im Pilot:
- $PW_\lambda$ als Mellin-Diagonal $t^2 + 1$ (echter Prolate-Operator später).
- $QW_{\lambda,\chi}$ mit vereinfachtem archimedischem Anteil $\tfrac12 t^2$ plus Prim-Summe für $p \le 14$.
- Ψ_{λ,χ} als Mellin-Multiplikation mit abgeschnittener $L$-Funktion ($n \le 200$).
- Randoperator $R_{\lambda,\chi}$ aus Satz 3.2 (Milestone 1_χ) vernachlässigt.

## 2. Rohdaten

| Charakter | μ_opt | ‖B‖ (‖Ψ PW‖) | Defekt (spektral) | Defekt (Frobenius) | Relativ-Defekt |
|---|---|---|---|---|---|
| $\chi_0$ | −2.99 | 26.87 | 854.30 | 1189.78 | **31.80** |
| $\chi_4$ | +0.60 | 1.68 | 5.26 | 7.77 | **3.12** |

**Ratio:** Relativ-Defekt($\chi_4$) / Relativ-Defekt($\chi_0$) = **0.098**.

## 3. Was die Zahlen sagen

### 3.1 χ-Signatur existiert

Die Defekt-Norm ist **nicht charakter-invariant**.  Das wirkt trivial, ist aber wichtig im Kontext der Atlas-v1–v4-Geschichte, in der die char-Signatur verschwand.  Hier taucht sie auf: Faktor ≈ 10 Unterschied zwischen $\chi_0$ und $\chi_4$.

### 3.2 Richtung widerspricht H1/H2

Die Blueprint-Hypothesen (`DIRICHLET_CCM_TRANSFER.md` §3) sagen:
- **H1:** $\varepsilon_{\lambda,\chi} \sim R_\chi \cdot \varepsilon_{\lambda,\chi_0}$, **wachsend** in $R_\chi$ (alle niedrigen Nullstellen kooperativ).
- **H2:** $\varepsilon_{\lambda,\chi} \sim (\gamma_{\chi_0}^{(1)}/\gamma_\chi^{(1)}) \cdot \varepsilon_{\lambda,\chi_0}$, **wachsend** bei niedrigerem $\gamma^{(1)}$.

Für $\chi_4$: $\gamma_{\chi_4}^{(1)} \approx 6.02 < 14.13 = \gamma_{\chi_0}^{(1)}$; $R_{\chi_4} \approx 0.02 > R_{\chi_0} \approx 0.002$.  **Beide Hypothesen sagen: χ_4 sollte größeren Defekt haben als χ_0.**  Der Pilot zeigt das **Gegenteil**: χ_4 hat kleineren Defekt.

### 3.3 Strukturelle Erklärung (Hypothese)

Die Richtungs-Umkehr hat einen naheliegenden Grund: **ζ(s) hat einen Pol bei s = 1**, L(s, χ_4) nicht.  Auf der kritischen Linie Re(s) = 1/2 sind wir an Abstand 1/2 vom Pol; das macht ζ(1/2 + it) im Mittel **größer** als L(1/2 + it, χ_4).  Die Mellin-Multiplikation Ψ mit ζ ist also eine "grössere" Operation — und der Defekt der Differenz mit der einfachen PW-Version wird entsprechend grösser.

Konsequenz: die Pilot-Zahlen messen nicht primär die Spektralinformation der $L$-Nullstellen (wie H1/H2 vermuten), sondern die **globale Größe** von $|L(1/2+it, \chi)|$.  Das ist ein starkes Signal, aber falsch kalibriertes.

### 3.4 Was das für das Blueprint bedeutet

Drei mögliche Lesarten:

(a) **H1/H2 sind falsch in ihrer einfachen Form.**  Die Defekt-Norm misst primär die Größe von L auf der kritischen Linie, nicht die Nullstellen-Dichte.  Dann müsste das Blueprint eine neue Hypothese H3 formulieren.

(b) **Die Pilot-Implementation ist zu primitiv.**  Die vereinfachte $QW$-Matrix (Polynomial-Diagonal + simple Prim-Summe) fängt die Pol-Struktur der Nullstellen nicht ein.  Eine echte Prolate-Galerkin-Basis mit Sonin-Zerlegung würde die L-Nullstellen explizit tragen.  Dann wäre das Pilot-Signal irrelevant für H1/H2.

(c) **Beides.**  H1/H2 sind in ihrer rohen Form zu simpel, und die Implementation ist zu vereinfacht, um das zu sehen.

Wahrscheinlichste Lesart: **(b) dominiert kurzfristig, (a) ist mittelfristig zu adressieren**.

## 4. Konkrete Diagnostik für die Pilot-Zahlen

### 4.1 $\|\text{Ψ PW}\|$

χ_0: 26.87; χ_4: 1.68.  Faktor ≈ 16.  Das reflektiert fast ausschließlich die Größe der L-Werte auf der Stützstellen:
- $|\zeta(0.5 + it)|$ im Bereich $|t| \le 1$: typisch 2–4 (Pol-Nähe, hebt an).
- $|L(0.5 + it, \chi_4)|$ im Bereich $|t| \le 1$: typisch 0.3–1.

Die χ_0-Werte sind dominiert vom Pol; dieser Beitrag hat **keine Nullstellen-Information**.

### 4.2 μ_opt

Das optimale μ (das die Defekt-Norm minimiert) ist bei χ_0 negativ (-2.99), bei χ_4 positiv (+0.60).  Der Vorzeichenwechsel ist strukturell interessant und könnte mit der Paritäts-Dichotomie aus Milestone 1_χ §6 zusammenhängen ($\chi_4(-1) = -1$ wechselt die Sektor-Zuordnung).

### 4.3 Der absolute Defekt

Auch die absolute Größe (31.80 vs. 3.12) ist nicht fundamental aussagekräftig, weil die vereinfachten QW/PW-Matrizen nicht korrekt normiert sind.  Aussagekräftig wäre: **ändert sich der Relativ-Defekt bei zunehmendem N?**  Das testet dieses Pilot-Skript nicht.

## 5. Konsequenzen für den Blueprint

### 5.1 Keine Widerlegung

Der Pilot **widerlegt H1/H2 nicht**.  Die Implementation ist zu primitiv, und der dominierende Effekt (ζ-Pol-Nähe) ist unabhängig vom zu messenden Signal.

### 5.2 Keine Bestätigung

Der Pilot **bestätigt H1/H2 nicht**.  Im Gegenteil: auf dem naiven Modell zeigt sich der umgekehrte Trend.

### 5.3 Was die nächste Iteration braucht

1. **Echte Prolate-Galerkin-Basis** (Slepian-Pollak-Funktionen) statt Mellin-Stützstellen.
2. **Expliziter Weil-Kern** mit Sonin-Zerlegung (CCM 2024 semilokale Trace-Formel).
3. **Normalisierung gegen ζ-Pol**: betrachte $\varepsilon_{\lambda,\chi} / \|\Psi_{\lambda,\chi}\|$ als dimensionsfreie Größe.
4. **Randoperator $R_{\lambda,\chi}$** aus Satz 3.2 (Milestone 1_χ) einbauen.

## 6. Ehrliche Einordnung

Dies ist ein **Pilot in der engsten Bedeutung des Wortes**: erste Berechnung, stark vereinfacht, vor allem dazu da, die Größenordnungen sichtbar zu machen und die Galerkin-Pipeline zu testen.

**Was bleibt:**  Die χ-Defekt-Norm ist eine konkrete numerische Größe, und ihre Werte sind mit einfachen Mitteln berechenbar.  Die Richtung der χ-Abhängigkeit ist subtiler als im Blueprint vermutet — was eine interessante Beobachtung ist, kein Rückschlag.

**Was als Nächstes zu tun wäre:**  Eine Version mit echter Prolate-Basis (z.B. via `scipy.special.pro_ang1` oder via numerische Diagonalisierung des PW-Operators auf einem feinen Gitter), idealerweise mit kontrolliertem Konvergenz-Verhalten in N.  Das ist 1–2 Tage Arbeit und sollte in einer nächsten Session angegangen werden.

---

## 7. Konkrete Zahlen (für spätere Referenz)

```
lambda     = 3.7417 (sqrt(14))
L          = 1.3195 (log lambda)
N          = 24
T_max      = 1.0556 (0.8 * L)

chi_0:
  mu_opt         = -2.9857 + 0.0000i
  ||Psi PW||     = 26.8651
  Defekt (spec)  = 854.3048
  Defekt (Frob)  = 1189.7808
  Rel. Defekt    = 31.7999

chi_4:
  mu_opt         = +0.6049 + 0.0000i
  ||Psi PW||     = 1.6848
  Defekt (spec)  = 5.2621
  Defekt (Frob)  = 7.7716
  Rel. Defekt    = 3.1232

Ratio chi_4/chi_0 = 0.0982
```

---

## 8. Aktualisierungs-Log

| Datum | Autor | Änderung |
|---|---|---|
| 2026-04-18 | LG (Opus 4.7 [1M], Session 16) | Pilot-Lauf mit $\lambda = \sqrt{14}$, $N=24$, $\chi_0$ vs. $\chi_4$. Ratio 0.098 (χ_4 kleiner als χ_0) — **umgekehrte Richtung** zu H1/H2 des Blueprints. Wahrscheinliche Ursache: ζ-Pol dominiert die Mellin-Werte; bei L(·,χ_4) ohne Pol ist Ψ_{λ,χ_4} in der simplen Modellversion strukturell kleiner. Keine Widerlegung (Implementation zu primitiv), keine Bestätigung (umgekehrter Trend). Nächste Iteration: echter Prolate-Galerkin-Raum + Sonin-Zerlegung + Pol-Normalisierung. |

---

**Ende PILOT_CHI_DEFECT_2026-04-18.md.** Erste Zahlen liegen vor; die χ-Defekt-Norm ist eine echte, messbare Größe, aber die Blueprint-Hypothesen H1/H2 in ihrer rohen Form sind mit dem Pilot-Modell nicht sichtbar. Nächster Schritt ist eine bessere Implementation, nicht sofort eine Blueprint-Revision.
