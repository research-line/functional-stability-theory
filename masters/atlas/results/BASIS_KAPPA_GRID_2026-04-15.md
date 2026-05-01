# Basis × Kappa Grid-Test: unerwartete Entkopplung

**Datum:** 2026-04-15 (Session 4)
**Script:** `_scripts/basis_kappa_full_grid.py`
**Status:** **UEBERRASCHENDER BEFUND** — widerspricht der ursprünglichen TWISTED_BASIS_ODD-Hypothese.

---

## 1. Versuchsaufbau

Voller Grid: 2 Charaktere × 2 Basen × 2 Kappas = 8 Kombinationen, 8 $\lambda$-Werte.

- **$\chi_4$** (odd, Legendre mod 4): $\chi(-1) = -1$.
- **$\chi_5$** (even, Legendre mod 5): $\chi(-1) = +1$.
- **Standard-Basis:** $\cos(\pi n t / L)$, $\sin(\pi n t / L)$ — Neumann-Rand im cos-Sektor.
- **Twisted-Basis:** $\cos(\pi(n+\tfrac{1}{2}) t / L)$, $\sin(\pi(n+\tfrac{1}{2}) t / L)$ — Dirichlet-Rand.
- **$\kappa = 0$:** even Gamma-Faktor, Diagonale mit $\psi_{\mathrm{dig}}(1/4 + i\tau/2)$.
- **$\kappa = 1$:** odd Gamma-Faktor, $\psi_{\mathrm{dig}}(3/4 + i\tau/2)$.

---

## 2. Ergebnisse (konsolidiert)

| $\chi$ | par($\chi$) | Basis | $\kappa$ | EVEN | ODD | mean gap | Dominanz |
|---|---|---|---|---:|---:|---:|---|
| $\chi_4$ | ODD | standard | 0 | 5 | 3 | $+0.24$ | EVEN |
| $\chi_4$ | ODD | standard | 1 | 4 | 4 | $+0.02$ | MIXED |
| $\chi_4$ | ODD | twisted  | 0 | 2 | 6 | $-0.23$ | **ODD** |
| $\chi_4$ | ODD | twisted  | 1 | 3 | 5 | $-0.04$ | ODD (weak) |
| $\chi_5$ | EVEN | standard | 0 | 7 | 1 | $+1.20$ | **EVEN** |
| $\chi_5$ | EVEN | standard | 1 | 6 | 2 | $+0.13$ | EVEN |
| $\chi_5$ | EVEN | twisted  | 0 | 2 | 6 | $-0.22$ | ODD |
| $\chi_5$ | EVEN | twisted  | 1 | 6 | 2 | $-0.01$ | EVEN (weak) |

---

## 3. Entscheidende Erkenntnis

### 3.1 Die Basis dominiert das Vorzeichen, nicht der Charakter

Die Zeilen mit **gleicher Basis und gleichem $\kappa$** zeigen fast **identisches Vorzeichen**, unabhängig vom Charakter:

| Basis | $\kappa$ | $\chi_4$-Gap | $\chi_5$-Gap | Korrelation |
|---|---|---:|---:|---|
| standard | 0 | $+0.24$ | $+1.20$ | **beide EVEN** |
| standard | 1 | $+0.02$ | $+0.13$ | **beide (schwach) EVEN** |
| twisted  | 0 | $-0.23$ | $-0.22$ | **beide ODD, fast identisch!** |
| twisted  | 1 | $-0.04$ | $-0.01$ | **beide ~0** |

**Interpretation:** Das Dominance-Vorzeichen ist **weitgehend charakter-unabhängig**. Es wird durch die Geometrie der Basis (Rand-Bedingungen) festgelegt, nicht durch die Arithmetik des Charakters.

### 3.2 Standard-Basis → EVEN-Dominanz

Für **beide** Charaktere gilt: Standard-Basis + $\kappa = 0$ gibt positiven Gap. Das ist nicht überraschend — es ist die **bereits bekannte** Riemann-artige Struktur.

Aber: auch mit $\kappa = 1$ (odd-Gamma) bleibt Standard-Basis EVEN-dominant, wenn auch schwächer. Der Gamma-Faktor allein **flippt das Vorzeichen nicht**.

### 3.3 Twisted + $\kappa = 0$ → ODD-Dominanz (!)

Das ist der eigentlich überraschende Befund: **die twisted basis mit even-Gamma zeigt starke ODD-Dominance** (beide Charaktere, mean gap $\approx -0.22$).

Der "Paritäts-Flip" entsteht also **nicht** durch den Charakter, sondern durch den Basis-Wechsel selbst. Die Rand-Bedingung (Dirichlet statt Neumann) inverts den Sektor-Gap.

### 3.4 Twisted + $\kappa = 1$ → Balanced

Die Kombination "twisted + odd-Gamma" gibt Gap $\approx 0$. Das deutet auf **Kompensations-Effekt**: Basis-Flip + Gamma-Flip heben sich (fast) auf.

---

## 4. Revidierte Meta-Meta-Struktur

Die einfache UBA-Vermutung aus `GAMMA_BASIS_DUALITAET.md` ist so **nicht haltbar**. Stattdessen:

### 4.1 Revidierte Regel (Hypothese v0.6)

> **Das Dominance-Vorzeichen ist ein Basis-Invariant, weitgehend unabhängig vom Charakter.**
>
> - Standard-Basis (Neumann-Rand): EVEN-dominiert.
> - Twisted-Basis (Dirichlet-Rand): ODD-dominiert (für $\kappa = 0$), balanciert für $\kappa = 1$.
>
> Die Charakter-Struktur $\chi$ verschiebt nur die **Stärke** der Dominance (über die Primsum-Frontier), nicht das Vorzeichen.

### 4.2 Konsequenz für v2.1-Universalität

Die Idee "v2.1 sieht den Gamma-Faktor" ist **differenzierter** als zuerst gedacht:

- v2.1 sieht den **Gamma-Faktor** via Diagonalelementen ($\kappa$).
- v2.1 sieht die **Basis** via Rand-Bedingungen.
- Diese beiden Effekte sind **nicht gekoppelt** wie vermutet.

Das impliziert: **die richtige Cos/Sin-Basis für Riemann ist** glücklich gewählt, weil sie mit $\kappa = 0$ (even-Gamma) kombiniert "EVEN-Dominanz" gibt. Für andere Gamma-Faktoren ist die Wahl der Basis **ein separater Parameter**, der getrennt untersucht werden muss.

### 4.3 Offene Fragen nach diesem Befund

1. **Gibt es eine Basis, die odd characters klar EVEN-dominiert?** Die twisted basis reicht nicht. Möglicherweise braucht man eine ganz andere Basis (z.B. Bessel- oder Hermite-basiert).

2. **Warum entkoppelt sich der Charakter vom Vorzeichen?** Das ist mathematisch unerwartet. Die Primsum-Gewichte sollten eigentlich das Vorzeichen beeinflussen. Eventuell ist das ein **Truncation-Artefakt** und bei $N \to \infty$ würde der Charakter wieder durchbrechen.

3. **Ist die Dirichlet-Rand-Bedingung der twisted basis der Grund?** Die Standard-Cos-Basis hat Neumann-Rand (Ableitung = 0); die twisted basis hat Dirichlet (Funktionswert = 0). Die konstante Funktion $\psi_0^{\mathrm{std}} = 1/\sqrt{2L}$ ist in der Standard-Basis, aber in der twisted basis fehlt sie → der EVEN-Sektor ist "magerer", was die Dominance ins Gegenteil kippen könnte.

---

## 5. Implikationen für Meta-Paper & Handoff

### 5.1 Was muss revidiert werden

- `TWISTED_BASIS_ODD.md`: Die Hypothese "twisted basis gibt odd chars EVEN-Dom" ist **falsch**. Die Datei muss als überholt gekennzeichnet werden.
- `GAMMA_BASIS_DUALITAET.md`: Die UBA-Vermutung muss reformuliert werden. Der Zusammenhang zwischen Gamma-Faktor und Basis ist **nicht** ein 1:1-Paar.
- Meta-Paper: keine unmittelbare Änderung nötig, da die Cartography-Zeile für odd characters ohnehin "$C2^{\mathrm{parity}}$ not directly transferable" sagt.

### 5.2 Was ist neu

Eine **Basis-Klassifikation** als separater Baustein. Die richtige Frage ist nicht "welche Basis zum Gamma", sondern "welche Basis erzwingt Rand-Bedingungen, die eine wohl-definierte Sektor-Asymmetrie erzeugen".

### 5.3 Empfehlung

Dieses Ergebnis ist ein **negatives**, aber wertvolles Resultat. Es zeigt:
- Die naive "Twisted-Basis"-Idee funktioniert nicht.
- Das v2.1-Programm für odd Dirichlet-Charaktere ist **strukturell offener** als zuerst gedacht.
- Weitere theoretische Arbeit muss die **Kopplung zwischen Charakter-Arithmetik und Basis-Geometrie** untersuchen.

---

## 6. Status

- **TWISTED_BASIS_ODD.md:** Hypothese falsifiziert. Datei bleibt als historisches Dokument mit Vermerk.
- **GAMMA_BASIS_DUALITAET.md:** UBA-Vermutung zu revidieren (v0.5).
- **Nächste Schritte:** 
  1. Warten auf Server-Asymptotik-Scan (läuft).
  2. Theoretisch untersuchen: warum dominiert die Basis so stark?

---

**Session 4 Teil 2 Fazit:** Das twisted-basis-Experiment falsifiziert eine naive Vorhersage, liefert aber eine **neue strukturelle Einsicht**: die Dominance ist ein Basis-Invariant, nicht ein Charakter-Invariant. Das ist ein nicht-triviales Detail, das die Meta-Meta-Struktur verfeinert.
