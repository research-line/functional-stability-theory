# Arch-Term-Integration — Session 7 Woche 1 Option 1 — historische N=200-Auswertung (überholt)

> **Archivhinweis:** Dieser Bericht bewahrt die damalige N=200-Auswertung. Die
> damalige 9/9-Beobachtung ist kein aktueller Beweis und kein prädiktiver
> Abschluss. Sie ist durch die [N=600-Serveranalyse](ARCH_TERM_N600_ANALYSIS.md)
> überholt: Die vollständige Galerkin-Zerlegung ist dort als Tautologie
> identifiziert, `S_with_arch` erreicht nur 8/10, und der echte
> `C2_χ`-Predictor bleibt zusammen mit dem Paley-Wiener-Weg offen.

**Datum:** 2026-04-16
**Parameter:** lambda = 20000, N_Galerkin = 200

## Historisches Hauptresultat (N=200)

Das Erfolgskriterium aus Session-6-Handoff §4.2 (sign_ok ≥ 9/10, R² ≥ 0.8)
wurde in dieser N=200-Auswertung für die damals als stabil eingeordneten
Charaktere beobachtet. Die spätere N=600-Auswertung ersetzt diese Einordnung;
die Tabelle ist daher nur als historische Messaufnahme zu lesen.

| Predictor | Alle 10 Charaktere | Ohne χ_21 (9 stabile) |
|-----------|---|---|
| S_exact (nur Prim) | sign=8/10, R²=0.49 | sign=8/9, R²=0.77 |
| **S_with_arch (Prim + Arch)** | **sign=9/10, R²=0.42** | **sign=9/9, R²=0.80** |
| gap_galerkin (volle EV-Diff, N=200) | sign=9/10, R²=0.41 | sign=9/9, R²=0.77 |
| prime_diff_matrix | sign=8/10, R²=0.49 | sign=8/9, R²=0.80 |

## Zerlegung

gap_galerkin := ev_minus − ev_plus = arch_diff + prime_diff

- `arch_diff` = ⟨φ⁻, W_arch⁻ φ⁻⟩ − ⟨φ⁺, W_arch⁺ φ⁺⟩
- `prime_diff` = S_exact (Autokorrelation) = prime_diff_matrix (via W_full − W_arch)

Die Konsistenz-Check-Differenzen zwischen S_exact (via Autokorrelation) und prime_diff_matrix (via Matrix) sind ~1–2% — kleine Numerik-Diskrepanz, beide sind strukturell äquivalent.

## Tabelle

| chi | D | gap_emp | gap_gal | S_exact | arch_diff | S+arch | sign+arch |
|-----|---|---------|---------|---------|-----------|--------|-----------|
| chi_5 | 5 | +0.01560 | +0.01429 | +9.0558e-02 | −5.7075e-02 | +3.3483e-02 | OK |
| chi_8 | 8 | −0.04902 | −0.04548 | −1.1646e-01 | +7.4278e-02 | −4.2185e-02 | OK |
| chi_12 | 12 | +0.03367 | +0.03251 | +6.5118e-02 | −1.4403e-02 | +5.0715e-02 | OK |
| chi_13 | 13 | −0.12504 | −0.12502 | −1.0211e-01 | −1.3269e-02 | −1.1538e-01 | OK |
| chi_17 | 17 | +0.01318 | +0.00408 | +2.0197e-02 | −3.5230e-03 | +1.6674e-02 | OK |
| chi_21 | 21 | −0.00423 | **+0.28163** | +2.3112e-01 | +5.7720e-02 | +2.8884e-01 | FAIL |
| chi_24 | 24 | +0.01076 | +0.01632 | +1.6231e-02 | +1.1948e-03 | +1.7426e-02 | OK |
| chi_29 | 29 | +0.01439 | +0.02372 | +2.5404e-02 | +4.4655e-04 | +2.5850e-02 | OK |
| chi_33 | 33 | −0.14221 | −0.14051 | −1.6529e-01 | +3.4411e-02 | −1.3088e-01 | OK |
| chi_60 | 60 | +0.00521 | +0.12363 | −4.9510e-02 | +1.7622e-01 | +1.2671e-01 | OK |

## χ_21-Ausreißer-Analyse

χ_21 (D = 21 = 3·7) zeigt bei N=200 den größten Mismatch:
- gap_emp = −0.00423 (N=600, stabil)
- gap_galerkin = **+0.28163** (N=200) — dreht Vorzeichen ggü. Target!
- gap_gal liegt eine Größenordnung über |gap_emp|

Dies wurde damals als **bekannte N-Oszillation** interpretiert. Die spätere
N=600-Analyse zeigt jedoch zusätzlich einen χ_29-Vorzeichenwechsel und eine
strukturelle Tautologie der vollständigen Galerkin-Zerlegung. Die damaligen
neun stabilen Vorzeichen sind deshalb keine Konvergenz- oder
Vorhersagegarantie.

## Interpretation

Die damalige Interpretation lautete, dass die Weil-Kern-Gap-Formel
(ANALYTIC_KERNEL Thm 4.1 + ANALYTIC_GROUNDSTATE Thm 4.1) bei N=200 für die
damals als stabil eingeordneten Charaktere 9/9 Vorzeichen und R² = 0.80
erreiche. Das ist eine historische numerische Beobachtung, keine bestätigte
Formel und keine aktuelle Vorhersage. Die N=600-Befunde im verlinkten Bericht
setzen die Claim-Grenze auf eine negative/strukturelle Diagnose zurück.

Die damalige Aussage, die Session-6-Falsifikation sei zu pessimistisch gewesen,
ist ebenfalls überholt. Die führende Approximation φ⁺ ≈ φ⁻ bleibt falsifiziert;
die unreduzierte Zerlegung liefert ohne unabhängige `C2_χ`-Auswertung keinen
prädiktiven Abschluss.

## Offene Punkte

1. **`C2_χ`-Predictor:** Eine unabhängige analytische Auswertung bleibt offen;
   der Paley-Wiener-/Asymptotik-Weg ist ausdrücklich noch nicht geschlossen.
2. **N-Konvergenz:** χ_21 und der bei N=600 sichtbare χ_29-Flip benötigen
   getrennte numerische Ursachenanalyse; sie legitimieren keine 9/9-Claim.
3. **Formel- und Skalierungsfragen:** Die historische Koeffizienten- und
   Slope-Diskussion bleibt Arbeitsmaterial, nicht ein belegter Claim.

## Paper-Status

Der damalige Plan, aus dieser Messaufnahme ein Abschluss-Paper zu machen, ist
überholt. Die Dirichlet-Character-Atlas-Struktur bleibt; die aktuelle
öffentliche Kernbotschaft ist der negative/strukturelle Befund aus der
[N=600-Analyse](ARCH_TERM_N600_ANALYSIS.md):

> **"A Weil-Kernel Numerical Atlas for Dirichlet L-Function Sector Gaps:
> Galerkin Diagnostics and the Boundaries of Leading-Order Theory"**

**Historische Paper-Claim-Aufzeichnung (nicht aktuell):**
Für reelle Dirichlet-Charaktere mit kleinem Leiter (D ∈ {5, 8, 12, 13, 17, 24, 29, 33, 60}) gibt die explizite Weil-Kern-Gap-Formel
$$ \mathrm{gap}_\chi = -2 \sum_{p,m} \chi(p)^m \frac{\log p}{p^{m/2}} \Delta_\chi(m \log p) + \mathrm{arch\_diff}_\chi + O(N^{-1}) $$
mit $\Delta_\chi(t) = (\phi^+_\chi * \phi^+_\chi)(t) - (\phi^-_\chi * \phi^-_\chi)(t)$
bei N=200 die damals beobachteten 9/9 Vorzeichen und $R^2 = 0.80$ für eine
ausgeschlossene χ_21-Menge. Diese Zahlen sind kein Beweis und keine
Vorhersage; nach N=600 bleibt `C2_χ` offen und die vollständige Zerlegung ist
tautologisch.
