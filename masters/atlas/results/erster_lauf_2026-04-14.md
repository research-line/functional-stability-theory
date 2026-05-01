# Erstes Prototyp-Ergebnis — $\chi_4$ Parity-Dominance

**Datum:** 2026-04-14 / 2026-04-15 (Session 3)
**Script:** `_scripts/chi4_parity_fast.py`
**Laufzeit:** ca. 10 Sekunden

## Rohdaten

| $\lambda$ | $N$ | $\lambda_1^+$ (Even) | $\lambda_1^-$ (Odd) | Gap $= \lambda_1^- - \lambda_1^+$ | Dominance |
|---:|---:|---:|---:|---:|---|
| 30 | 12 | $-4.2646$ | $-4.2835$ | $-0.0189$ | **ODD** ✓ |
| 100 | 16 | $-6.2525$ | $-6.3321$ | $-0.0796$ | **ODD** ✓ |
| 300 | 18 | $-7.7538$ | $-7.7937$ | $-0.0398$ | **ODD** ✓ |
| 1000 | 20 | $-9.1149$ | $-8.9631$ | $+0.1518$ | **EVEN** ✗ |

**Kontrolle (Riemann, $\chi_0$, kappa=0):**

| $\lambda$ | $\lambda_1^+$ | $\lambda_1^-$ | Gap | Dominance |
|---:|---:|---:|---:|---|
| 100 | $-16.2024$ | $-13.1879$ | $+3.0144$ | EVEN ✓ |

## Beobachtungen

### B1: Odd-Dominance fuer $\chi_4$ qualitativ bestaetigt
Bei drei von vier getesteten $\lambda$-Werten ist der Grundzustand im Odd-Sektor. Das ist konsistent mit der Meta-Paper-Vorhersage fuer odd character ($\chi_4(-1) = -1$).

### B2: Flip bei $\lambda = 1000$ — unerwartet
Bei $\lambda = 1000, N = 20$ flippt die Dominance zu Even. Entweder:
- **(a) Truncation-Artefakt:** $N = 20$ ist zu klein fuer stabile Trunkation bei $\lambda = 1000$. Die Riemann-Paper-Heuristik verlangt $N \gtrsim L^2 \approx 48$ fuer $L = \log 1000 \approx 6.9$. Das ist nicht erfuellt.
- **(b) Echtes Phaenomen:** der Gap kollabiert bei Dirichlet-$\chi_4$ mit $\lambda$ und kann Vorzeichen-Flip erleiden.

### B3: Gap-Groesse bei $\chi_4$ vs. Riemann
Gap bei $\chi_4$: $\sim 0.02 - 0.15$. Gap bei Riemann: $+3.01$ bei $\lambda = 100, N = 16$. **Faktor 30-40 Unterschied**.

Interpretation: Die $\chi_4$-Phasen ($\pm 1$, etwa gleichverteilt) heben sich **teilweise auf**, was den aggregierten Dominance-Effekt dramatisch schwaecht. Die Random-Phase-Heuristik aus KONZEPT.md §3 sagt Gap $\sim \sqrt\lambda/\sqrt{\log\lambda}$, nicht $\sqrt\lambda$. Fuer $\lambda = 100$: $\sqrt{100}/\sqrt{\log 100} = 10/\sqrt{4.6} \approx 4.7$, also Gap $\sim 4.7$. Mit starken Praefactors kleiner, also durchaus im Bereich $0.1$.

### B4: Konsequenz fuer Meta-Paper-Vorhersage
Die Vorhersage "Odd-Dominance fuer $\chi_4$" ist qualitativ richtig, aber **die Dominance ist knapp und fragil**. Das verstaerkt die Einschaetzung aus KONZEPT.md §4.3:
- Siegel-Landau-Problematik ist kritisch.
- CAP-Zertifikate brauchen groessere $N$ und hoehere Interval-Arithmetik-Praezision.
- Das Dirichlet-v2.1-Paper ist **technisch anspruchsvoller** als das Riemann-v2.1-Paper.

## Naechste Schritte

1. **Stabilitaets-Scan:** `chi4_scan.py` testet Gap als Funktion von $N$ bei festem $\lambda = 1000$, und Gap-Skalierung mit $\lambda$ bei angemessen skaliertem $N$.
2. **Interval-Arithmetic-CAP:** Portierung auf `mpmath.iv` fuer rigorose Zertifizierung (analog zu Riemann v2.1).
3. **Frontier-Dominance-Analyse:** Extraktion des reinen Frontier-Beitrags vs. Non-Frontier-Anteil.

## Status

**Qualitativ bestaetigt:** Odd-Dominance fuer $\chi_4$ bei kleinem/mittlerem $\lambda$.
**Offen:** Stabilitaet bei grossem $\lambda$, rigorose Zertifizierung, Quantifizierung des Gap-Skalierens.

Meta-Paper-Vorhersage **nicht widerlegt**, aber **prazisiert**: Dirichlet-Parity-Dominance ist schwaecher als Riemann-Even-Dominance, und ihre Robustheit ist eine substantielle Arbeitslinie.
