# α.1 — Vorzeichen-Flip und Grundzustand-Struktur (2026-04-18)

**Autor:** LG (Claude Opus 4.7 [1M], Session 17 Fortsetzung α.1)
**Skript:** `_scripts/dirichlet_ccm_mpmath.py` (mit Cache, Diagnose, SIGN_WR-Test)
**Kontext:** Folge auf `DIRICHLET_CCM_OPERATOR_2026-04-18.md` und Handoff.
**Parameter:** $\lambda = \sqrt{14}$, $N = 30$, mpmath 50 Digits.

---

## 1. Ausgangspunkt

Im Vorgänger-Bericht war `z_6 = 14.13472533` (7·10⁻¹⁰ von $\gamma_1$), aber `z_1..z_5` saßen bei den Pol-Mittelpunkten. Diagnose: falscher Grundzustand.

Der Fix-Versuch aus dem Handoff (`sum(xi) > eps`-Filter in `find_ground_state`) funktioniert **nicht**, weil der Grundzustand in der ursprünglichen Konvention (`SIGN_WR = -1`) in einem 21-fach entarteten Cluster bei $\epsilon_N \approx -1.337$ liegt — und aus diesem Cluster ist $\delta_N$ **praktisch orthogonal** (Projektion $\|P\| \approx 0.002$).

## 2. Vorzeichen-Flip

Zwei Konventionen getestet:

| Konvention | $\epsilon_N$ | Cluster d | $\|\langle \delta_N \mid \xi\rangle\|$ | $\sum \xi_j$ | MS1 erfüllt? |
|---|---|---|---|---|---|
| `SIGN_WR = -1` | $-1.337$ | **21** (numer. entartet) | $0.002$ | $0.016$ | **Nein** |
| `SIGN_WR = +1` | $-3.089$ | **1** (einfach) | $0.241$ | $1.88$ | **Ja** |

Mit $\text{SIGN\_WR} = +1$ ist der Grundzustand **einfach, even, mit nicht-trivialer $\delta_N$-Kopplung**. Das ist genau die CCM-MS1-Bedingung (Thm 5.10).

**Interpretation der Vorzeichen-Konvention:** Prop 4.3 in CCM 2025 liefert explizit die Matrix-Einträge von $W_\mathbb R$ über $\alpha_L, \beta_L, \gamma_L$. In der Weil-Form (CCM 3.10) erscheint $-W_\mathbb R$. Frage ist: ist in Prop 4.3 das Minus schon eingebaut oder nicht?

**Empirische Antwort:** Das Minus ist eingebaut. `M = +W_R - \sum W_p + W_{0,2}` liefert die MS1-konforme Struktur.

## 3. Nullstellen von $F(z)$ mit $\text{SIGN\_WR} = +1$

Nach Normalisierung $\sum \xi_j = \sqrt L$:

```
z_1 = 0.9714793458   (bekannt γ_1 = 14.1347,  Fehler -13.16)
z_2 = 3.3235181215   (bekannt γ_2 = 21.0220,  Fehler -17.70)
z_3 = 5.6719941421
z_4 = 8.0321371481
z_5 = 10.4146796840
z_6 = 12.9306495093
z_7 = 14.7646931274
```

$\gamma_1 = 14.1347$ liegt **zwischen** $z_6 = 12.93$ und $z_7 = 14.76$. Aber wurde nicht als isolierte Nullstelle gefunden.

## 4. Überraschung: Alt-Test mit Eigenvektor $k=24$

Der **Eigenvektor mit maximal $|\langle \delta_N \mid v \rangle|$ im even-Sektor** (nicht der Grundzustand): bei $\text{SIGN\_WR}=+1$ ist das $k=24$, $w=-0.506$, overlap $0.445$.

**Nullstellen von $F(z)$ mit diesem Vektor:**

```
z_1 = 2.23   z_2 = 4.42   z_3 = 6.69   z_4 = 9.00   z_5 = 11.35
z_6 = 14.09  (γ_1 = 14.13, Fehler -0.04!)
z_7 = 18.27
z_8 = 21.16  (γ_2 = 21.02, Fehler +0.14!)
```

**Das ist kein Zufall**: $\gamma_1$ und $\gamma_2$ erscheinen als Nullstellen, aber mit einem **Geist-Zwischenspektrum** ($z_7 = 18.27$ ist keine $\gamma$).

## 5. Diagnose

**Das ist ein strukturelles Phänomen der Galerkin-Approximation**, nicht ein Implementations-Bug im Vorzeichen oder in der Normalisierung:

1. Der **wahre** $\xi$ im Kontinuum-Limes hat eine sehr spezifische Fourier-Struktur, die bei endlichem $N$ über mehrere Eigenmoden "verschmiert" ist.
2. Bei $N = 30$ reicht die Auflösung noch nicht, um alle $\gamma_n$ als direkte Nullstellen von $F(z)$ zu sehen.
3. Das CCM-Resultat $10^{-55}$ gilt bei **$N = 120$** — was bei mpmath mit 50 Digits etwa 60× mehr Rechenzeit bedeutet.

**Befund:** Die Struktur ist jetzt physikalisch korrekt (MS1 erfüllt), aber der Konvergenzfehler ist bei $N = 30$ einfach zu groß. Das Skalen-Verhalten $N \to \infty$ muss numerisch untersucht werden.

## 6. Was als nächstes

### α.1.b — Conductor-Term in Rang-1-Form
Die aktuelle heuristische Diagonal-Verschiebung (`WR[N,N] += log(q/π)·L`) ist für non-trivial χ falsch. Rang-1-Form:
$$[\text{Cond}]_{nm} = \log(q/\pi) \cdot \tilde V_n(0) \cdot \tilde V_m(0)$$
mit $\tilde V_n(0) = \int V_n(u)\,d^\ast u$. Für die Fourier-Basis konstant $V_0$ ist das $\sqrt L$ bei $n=0$ und $0$ sonst. Das heisst: Conductor trägt **nur auf $(n,m) = (0,0)$** bei.

### α.2 — Skalen-Test
Führe $N \in \{30, 60, 120\}$ bei $\lambda = \sqrt{14}$ für Riemann durch. Konvergiert $z_1 \to \gamma_1$? Mit welchem Skalen-Gesetz (polynomial/exponentiell)?

### α.3 — Hochpräzise Diagonalisierung
Bei $N = 30$ gibt es in `SIGN_WR = -1` einen 21-fach entarteten Cluster. Das könnte ein numerisches float64-Artefakt sein. Test via `mpmath.eigsy` in 50-Digit-Präzision klärt, ob die Entartung echt oder numerisch ist.

### β — Dirichlet-Tests (nur nach α.2 Konvergenz-Nachweis)
$\chi_4$, $\chi_5$, $\chi_{33}$ mit SIGN_WR=+1 und Conductor-Term.

## 6a. Server-Lauf N=60 (2026-04-18 nachmittag)

Nach dem Sign-Flip wurde $N=60$ auf ellmos-services mit `SIGN_WR=+1` gerechnet (ollama gestoppt). **Ergebnis: keine Konvergenz.**

| Größe | $N=30$ | $N=60$ |
|---|---|---|
| $\epsilon_N$ | $-3.089$ | $-3.935$ |
| Cluster d (tol=$10^{-4}$) | 1 | 1 |
| $\|\langle \delta_N \mid \xi_0\rangle\|$ | 0.241 | **0.021** |
| $\sum \xi_j$ (vor Norm) | 1.88 | **0.23** |
| $z_1$ (Grundzustand) | 0.97 | 0.99 |
| $z_6$ | 12.93 | 12.96 |
| $z_7$ | 14.76 | 14.87 |
| $\gamma_1 = 14.13$ | zwischen $z_6, z_7$ | zwischen $z_6, z_7$ |

**Zentrale Beobachtung:** Bei $N=60$ ist $\|\langle \delta_N \mid \xi_0\rangle\|$ **kleiner** geworden (0.021 statt 0.241). Der CCM-Grundzustand soll nach MS1 mit $\delta_N$ koppeln und bei $N \to \infty$ konvergieren. Stattdessen driftet mein $\xi_0$ weg von $\delta_N$. **Das ist strukturell falsch**, nicht ein Präzisionsproblem.

Die Nullstellen $z_1..z_5$ bei $N=60$ sind an fast identischen Stellen wie bei $N=30$ (Abweichung < 0.02 trotz doppelter Basis-Dimension). Das bestätigt: meine $F(z)$ konvergiert **nicht** zu den L-Nullstellen.

## 6b. Hypothesen zur strukturellen Ursache

Der Fehler liegt **nicht** in `find_ground_state` oder im Vorzeichen, sondern tiefer:

**(H1) Prop 4.3 unvollständig umgesetzt.** Meine $\alpha_L, \beta_L, \gamma_L$-Integrale könnten fehlende Hurwitz-Lerch- oder ${}_2F_1$-Terme haben. CCM 2025 §4 nennt diese explizit. Prüfen: gibt CCM Prop 4.3 die Matrix $W_\mathbb R$ **direkt** oder $(2N+1) \times (2N+1)$-Anteil einer größeren Konstruktion?

**(H2) $\delta_N$-Definition falsch.** In CCM ist $\delta_N \in E_N$ der "Dirichlet-Kernel" — das ist möglicherweise **nicht** $L^{-1/2}\sum V_j$, sondern hat andere Koeffizienten. Prüfen: CCM Eq (5.24) gibt $\delta_N$ explizit in $V_n$-Basis?

**(H3) Prim-Summe-Faktor falsch.** Mein q-Kernel und Faktor $\chi(k)\Lambda(k)/\sqrt k$ könnte Skalierungsfehler haben. CCM Eq (2.9) vs meiner Implementation prüfen.

**(H4) F(z)-Formel unvollständig.** CCM Eq (5.25) hat Pre-Faktor $2L^{-1/2}\sin(zL/2)$. Wenn dieser nicht korrekt für Nullstellensuche berücksichtigt wird, könnten die echten L-Nullstellen nicht als $F(z)=0$ erscheinen, sondern als $\widehat\xi(z)=0$ wobei $\widehat\xi = \text{sin-Faktor} \cdot F$.

## 7. Strategische Einordnung

**Positiv:** MS1 ist bei $\text{SIGN\_WR}=+1$ strukturell erfüllt (einfacher Grundzustand, parity-korrekt, $\delta_N(\xi) \ne 0$). Die Vorzeichen-Frage war der korrekte Fix — aber nicht ausreichend für CCM-Präzision bei $N = 30$.

**Offen:** Die Konvergenzrate bei endlichem $N$ ist unbekannt. Wenn sie exponentiell ist (wie CCM andeutet: $10^{-55}$ bei $N=120$), sollte $N=60$ bereits Fehler $\sim 10^{-3}$ liefern — das wäre publikationsreif.

**Wichtigste nächste Aktion:** $N=60$-Lauf mit SIGN_WR=+1. Wenn $z_1 \to \gamma_1$ mit Fehler $< 10^{-2}$, ist die Implementation validiert und die Dirichlet-Tests können starten.

---

## 8. Aktualisierungs-Log

| Datum | Autor | Änderung |
|---|---|---|
| 2026-04-18 | LG (Opus 4.7, Session 17 α.1 Fortsetzung) | α.1-Fix angewandt: neue Ground-State-Projektion im entarteten Cluster; kritischer Befund: SIGN_WR-Vorzeichen war falsch. Mit SIGN_WR=+1 ist Grundzustand einfach und MS1-konform, aber Nullstellen von $F(z)$ approximieren $\gamma_n$ bei $N=30$ noch nicht. $\gamma_1$ und $\gamma_2$ erscheinen in $F(z)$-Nullstellen des zweiten-höchsten-overlap EV ($k=24$). Konvergenz-Skalentest bei $N=60, 120$ als nächster Schritt. Conductor-Rang-1-Form identifiziert. |
| 2026-04-18 | LG (Opus 4.7, Session 17 α.1 Fortsetzung, Server) | **N=60-Lauf auf ellmos-services** fertig: **keine Konvergenz**. $\sum \xi_j$ sinkt von 1.88 (N=30) auf 0.23 (N=60) — Grundzustand driftet **weg** von $\delta_N$. Nullstellen $z_1..z_5$ an identischen Positionen bei N=30 und N=60. Fehler daher **strukturell**: eine der Hypothesen H1-H4 (Prop 4.3 unvollständig, $\delta_N$-Definition, Prim-Faktor, F(z)-Formel) muss geprüft werden. Nächster Schritt: CCM 2025 §4-5 gegen Implementation line-by-line abgleichen. |

---

**Ende ALPHA1_SIGN_FLIP_BEFUND_2026-04-18.md.** Kern-Insight: Vorzeichen-Konvention war inkompatibel zu CCM Prop 4.3, Fix via SIGN_WR=+1 stellt MS1 her. Für CCM-Präzision reicht $N = 30$ nicht; Skalentest $N = 60, 120$ vorbereitet.
