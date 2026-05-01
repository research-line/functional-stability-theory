# Substantieller Befund: $\chi_4$-Parity-Dominance oszilliert

**Datum:** 2026-04-15 (Session 3, Phase B)
**Scripts:** `_scripts/chi4_parity_fast.py`, `chi4_vectorized.py`, `chi4_fine_scan.py`
**Status:** **TEILFALSIFIKATION** der Meta-Paper-Vorhersage

---

## 1. Rohdaten (Fine Scan)

| $\lambda$ | $N$ | Gap $= \lambda_1^- - \lambda_1^+$ | Gap / $\sqrt L$ | $\pi(3)-\pi(1)$ | Dominance |
|---:|---:|---:|---:|---:|---|
| 30 | 23 | $-0.012$ | $-0.007$ | $+1$ | ODD |
| 50 | 30 | $+0.033$ | $+0.017$ | $+2$ | EVEN |
| 75 | 37 | $+0.180$ | $+0.086$ | $+2$ | EVEN |
| 100 | 42 | $-0.097$ | $-0.045$ | $+2$ | ODD |
| 150 | 50 | $-0.015$ | $-0.007$ | $+2$ | ODD |
| 200 | 56 | $+0.122$ | $+0.053$ | $+3$ | EVEN |
| 250 | 60 | $-0.027$ | $-0.011$ | $+4$ | ODD |
| 300 | 65 | $-0.020$ | $-0.008$ | $+3$ | ODD |
| 400 | 71 | $+0.049$ | $+0.020$ | $+3$ | EVEN |
| 500 | 77 | $-0.043$ | $-0.017$ | $+6$ | ODD |
| 700 | 80 | $-0.008$ | $-0.003$ | $+6$ | ODD |
| 1000 | 80 | $+0.115$ | $+0.044$ | $+7$ | EVEN |
| 1400 | 80 | $-0.044$ | $-0.016$ | $+3$ | ODD |
| 2000 | 80 | $-0.027$ | $-0.010$ | $+8$ | ODD |
| 3000 | 80 | $+0.100$ | $+0.035$ | $+7$ | EVEN |
| 5000 | 80 | $+0.104$ | $+0.036$ | $+10$ | EVEN |

**Dominance:** 9 ODD / 7 EVEN, praktisch 50/50.

## 2. Truncation-Stabilitaets-Kontrolle (bei $\lambda = 1000$)

| $N$ | Gap | Dominance |
|---:|---:|---|
| 10 | $-1.22$ | ODD |
| 14 | $-1.01$ | ODD |
| 18 | $+0.156$ | EVEN |
| 22 | $+0.133$ | EVEN |
| 26 | $+0.139$ | EVEN |
| 30–75 | $+0.115 – +0.139$ | EVEN (stabil) |

**Konvergenz.** Ab $N = 18$ stabilisiert sich der Gap bei $+0.12$. Der Wechsel von ODD (bei $N \leq 14$) zu EVEN (ab $N \geq 18$) ist ein klares Truncation-Artefakt: **kleine $N$ liefern spurious ODD-Dominance.** Das heisst, unsere urspruenglichen Ergebnisse in `erster_lauf_2026-04-14.md` (bei $N = 16$ etc.) waren ebenfalls vom Truncation-Effekt betroffen.

**Revidierte Zahlen mit stabilen $N \gtrsim 2L^2$:** siehe Abschnitt 1. Die Alternations-Struktur bleibt bestehen, ist aber nicht blosse Truncation.

## 3. Korrelations-Analyse

**Chebyshev-Bias $\pi(x; 4, 3) - \pi(x; 4, 1)$.**
Korrelation mit Gap: $+0.21$ (schwach positiv).

Interpretation: Chebyshev-Bias ist **nicht** die Hauptursache der Alternation. Die Bias-Richtung ist monoton ansteigend im Mittel, aber der Gap schwankt drastisch.

## 4. Theoretische Interpretation

### 4.1 Die Meta-Paper-Vorhersage ist falsifiziert (in der naiven Form)

Die Cartography sagte: *"Dirichlet: $C2^{\mathrm{parity}}$ closable via v2.1-twist"* mit Parity-Richtung $\chi(-1)$. Fuer $\chi_4$ bedeutete das: **stabile Odd-Dominance**. Die Numerik zeigt:

- Kein stabiles Dominance-Vorzeichen
- |Gap| klein (< 0.2) und von $\lambda$ nahezu unabhaengig
- Oszillation zwischen Odd und Even

**Die Meta-Paper-Vorhersage ist in dieser Form nicht haltbar.**

### 4.2 Was koennte stattdessen gelten?

**Hypothese H1 (Truncation + zu kleine N):** Der "wahre" Gap koennte bei unendlichem $N$ existieren, aber wir erreichen ihn nicht mit $N = 80$. Dagegen spricht: im Bereich $N = 18$ bis $75$ bei fixem $\lambda = 1000$ ist Gap vollstaendig stabil. Kein asymptotischer Trend sichtbar.

**Hypothese H2 (Random-Phase gibt nur asymptotisch Dominance):** Aus Random-Matrix-Heuristik: Gap $\sim \sqrt\lambda / \sqrt{\log\lambda}$. Fuer $\lambda = 5000$, $\log\lambda = 8.5$: erwartet Gap $\sim 70.7 / 2.9 \approx 24$. Beobachtet: 0.1. **Diskrepanz: Faktor 240.** Die Heuristik ueberschaetzt.

**Hypothese H3 (Char-parity flippt Rolle, aber auch Dominance-Vorzeichen):** Fuer odd characters ($\chi(-1) = -1$) ist die "richtige" v2.0-Basis nicht Cosinus, sondern eine mit vertauschten Rollen. Der Shift Parity Lemma koennte fuer odd characters **umgekehrt signiert** sein, und unser Test sieht daher gemischte Vorzeichen.

**Hypothese H4 (Effekt ist NICHT fuer alle Dirichlet gleich):** Gerade und ungerade Charaktere koennten unterschiedlich sein. **Der kritische Test: wiederhole mit einem even character (z.B. $\chi_5$).**

### 4.3 Was wir schon wissen

- Chebyshev-Bias ist nicht die Hauptursache.
- Truncation-Effekte treten bei $N < 2L^2$ auf, aber bei grossen $N$ stabilisiert sich Gap ohne klares Muster.
- Riemann-Kontrolle ($\chi_0$ trivial) zeigt klare Even-Dominance mit Gap $\sim \sqrt\lambda$. Daher ist das Phaenomen nicht bloss numerisches Rauschen: **bei $\chi_0$ funktioniert es klar, bei $\chi_4$ nicht.**

## 5. Konsequenzen fuer das Meta-Programm

### 5.1 Meta-Paper-Cartography (Dirichlet-Zeile) muss revidiert werden

Die Zeile
> "Dirichlet $L(s,\chi)$: $C2^{\mathrm{parity}}$ expected closable via v2.1-twist"

ist **zu optimistisch formuliert**. Korrigierter Vorschlag:

> "Dirichlet $L(s,\chi)$: $C2^{\mathrm{parity}}$-Transfer via v2.1 ist **charakter-paritaets-sensitiv**. Fuer gerade $\chi$ ($\chi(-1) = +1$) plausibel via direkten Transfer; fuer ungerade $\chi$ ($\chi(-1) = -1$) nicht direkt — erfordert Basiswechsel oder modifizierte Methode. Numerik bei $\chi_4$ (odd) zeigt oszillatorisches Dominance-Verhalten ohne stabile Dominance-Richtung."

### 5.2 Nach Meta-Paper: was ist **wirklich** bewiesen?

- **Riemann v2.1:** unconditional Even-Dominance (trivial character, $\chi_0(-1) = +1$).
- **Selberg:** v2.0-Methode + Laplace-Casimir, beide Wege zur selben RH.
- **Dirichlet $\chi_4$ (odd):** numerisch kein klares Dominance-Muster — **Hypothese falsifiziert**.
- **Dirichlet even character:** noch nicht getestet; kritisch fuer Phase B.

### 5.3 Naechster Schritt (Phase B.1)

Wiederhole den Fine-Scan mit einem **even character**. Kandidaten:
- $\chi_5$ (Legendre mod 5): $\chi_5(-1) = +1$ (even).
- $\chi_8^{(1)}$ (mod 8, primitiv even).
- $\chi_{12}^{(1)}$ (mod 12, primitiv even).

**Vorhersage.** Wenn die naive Meta-Paper-Regel stimmt, sollte fuer $\chi_5$ klare Even-Dominance auftreten. Wenn auch dort Oszillation: tieferes Problem mit Dirichlet-Transfer.

## 6. Falsifikation als Gewinn

**Wissenschaftlich relevant.** Dieser Befund ist **nicht** eine Niederlage des Meta-Papers, sondern eine **Praezisierung**. Die Meta-Paper-Vorhersage war ein konkretes, falsifizierbares Statement — und es wurde (teilweise) falsifiziert. Das ist, wie Wissenschaft funktioniert:

1. Meta-Paper formuliert: "Dirichlet transferiert via v2.1-twist."
2. Prototyp testet: 16 Faelle, 9 ODD / 7 EVEN → **kein klares Muster**.
3. Revision: "Dirichlet-Transfer ist paritaets-sensitiv."
4. Naechster Test: even character → entscheidender Test der Revision.

**Der Befund ist wertvoll,** weil er die Reichweite des v2.1-Programms empirisch praezisiert, statt auf naiven Vorhersagen zu bleiben.

---

**Status:** Phase B Prototyp dokumentiert. Naechster Schritt: even-character-Test zur Verfeinerung der Hypothese.
