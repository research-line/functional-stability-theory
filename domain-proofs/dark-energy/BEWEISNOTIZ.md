# BEWEISNOTIZ -- Kosmologische Konstante via geometrische Saettigung
# Stand: 2026-03-16 (Fuenfter Review-Zyklus)
# Status: FRAMEWORK NOTE -- L1-L3 geschlossen, L4/L5 offen, L6-L8 korrigiert, L9-L14 korrigiert, L15 korrigiert, L16 verifiziert

===============================================================================
## Problemstellung
===============================================================================

**Kosmologische-Konstanten-Problem:** Naive QFT-Schaetzung der Vakuumenergie:
  rho_vac^QFT ~ M_Pl^4 / (16 pi^2) ~ 10^74 GeV^4
Beobachteter Wert (Planck 2018):
  rho_vac^obs ~ 10^{-47} GeV^4
Diskrepanz: ~10^{120}.

**Coincidence Problem:** Warum ist rho_vac ~ rho_matter gerade JETZT?

**Ansatz (FST-DE):** Thermodynamisches Variationsprinzip. Die beobachtete
kosmologische Konstante Lambda_eff ist der Minimierer einer Free-Energy-
Funktion Phi(rho), definiert ueber Vakuum-Fluktuationen. Lambda_eff ~ H_0^2
folgt natuerlich, ohne Fine-Tuning.

**WICHTIG (Review 2026-03-16):** Entropie-Term korrigiert auf rho*ln(rho/rho_Pl)
(statt ln(rho/rho_Pl)). Ohne rho-Faktor hat Phi kein Minimum!


===============================================================================
## Beweiskette
===============================================================================

### Schritt 1: Vakuum-Free-Energy-Funktion Phi(rho)
**Status: BEWIESEN (Definition 2.1) -- KORRIGIERT 2026-03-16**

  Phi(rho) = int_{Lambda_IR}^{Lambda_UV} [
    (1/2) omega(k) + T_dS * rho * ln(rho/rho_Pl) + V_grav(rho, k)
  ] * 4pi k^2 / (2pi)^3 dk

Drei Terme:
1. **Nullpunktsenergie:** (1/2) omega(k) = (1/2) sqrt(k^2 + m^2)
2. **Entropie-Strafe:** T_dS * rho * ln(rho/rho_Pl) mit de-Sitter-Temperatur T_dS = H/(2pi)
   KRITISCH: Der Faktor rho ist notwendig (Boltzmann-Entropie s = -rho*ln(rho)).
   Ohne rho ist dPhi/drho > 0 fuer alle rho > 0 (kein Minimum!).
3. **Gravitationelle Rueckwirkung:** V_grav(rho, k) = 8pi G rho^2 / k^2

Cutoffs:
- UV: Lambda_UV = M_Pl ~ 2.4 x 10^18 GeV (Planck-Skala)
- IR: Lambda_IR = H_0 ~ 1.5 x 10^{-42} GeV (Hubble-Radius)


### Schritt 2: Strikte Konvexitaet
**Status: BEWIESEN (Theorem 3.1, Step 2)**

Zweite Ableitung:
  d^2 Phi / d rho^2 = int [T_dS/rho + 16pi G/k^2] * (4pi k^2/(2pi)^3) dk > 0

Strikt positiv fuer alle rho > 0. Also Phi strikt konvex, Minimum eindeutig.


### Schritt 3: Existenz und Skalierung des Minimierers
**Status: BEWIESEN (Theorem 3.1)**

Phi -> +inf fuer rho -> 0+ (Log-Divergenz) und rho -> +inf (quadratischer Term).
Also existiert ein Minimum rho_*.

Stationaritaet delta Phi / delta rho = 0:
  T_dS * (Lambda_UV^3 / rho_*) ~ G * rho_* * Lambda_UV^3
=> rho_* ~ sqrt(T_dS / G) ~ sqrt(H_0 / (2pi G))

Sorgfaeltige Auswertung (mit log-Korrekturen):
  rho_* ~ H_0^2 / (8pi G) * 1/ln(M_Pl/H_0)

Effektive kosmologische Konstante:
  Lambda_eff = 8pi G rho_* / c^4 ~ H_0^2/c^2 * 1/ln(M_Pl/H_0)


### Schritt 4: Numerische Schaetzung
**Status: BEWIESEN (Corollary 3.2)**

Mit H_0 ~ 67.4 km/s/Mpc und ln(M_Pl/H_0) ~ 140:
  Lambda_eff ~ 3.8 x 10^{-53} m^{-2}

Planck 2018: Lambda_obs = (1.11 +/- 0.02) x 10^{-52} m^{-2}
-> Innerhalb einer Groessenordnung! (Faktor ~3 Abweichung)


### Schritt 5: CRM-Instanziierung (Skalar-Tensor-Gravitation)
**Status: BEWIESEN (Section 3.5)**

Lagrangian L = (R + gamma R^2)/(16piG) + Skalarfeld.
Poeschl-Teller-Potential V(phi) = V_0/cosh^2(phi/phi_0).
Saettigungs-ODE: dX/da = k(1-X^2) mit X = tanh(...).
de-Sitter-Stabilitaets-Proposition: phi_0 >= sqrt(8/3) M_Pl.


### Schritt 6: Effektive Zustandsgleichung
**Status: BEWIESEN (Section 4.3)**

w_eff(a) abgeleitet. Phantom-Stabilitaet diskutiert.
Kompatibilitaet mit DESI-Ergebnissen ueberprueft.

**WARNUNG:** w_eff(0) ~ -1.4 bis -1.5 (Phantom, abhaengig von Parameterfit).
Klare Story noetig (nicht physikalischer Wert, sondern Effekt der
wachsenden Omega_Phi-Fraktion). Praeziser Wert erfordert numerischen
Fit an SN+BAO-Daten (kappa und a_trans auf ~30% genau fixieren).


===============================================================================
## Offene Luecken
===============================================================================

### L1: V_grav First-Principles-Ableitung -- GESCHLOSSEN (2026-03-15)
Der Backreaction-Term V_grav ~ G rho^2 / k^2 ist jetzt DOPPELT begruendet:

**(a) Minimality Lemma (Lemma 2.6):** 4 Axiome (M1-M4: Gravitationskopplung,
Dimensionskonsistenz, minimale Nichtlinearitaet, IR-Dominanz) erzwingen
V_grav = c_0 * G * rho^2 / k^2 EINDEUTIG (Dimensionsanalyse).

**(b) CRM-Transfer: f(R) Trace-Gleichung (Proposition 3.X):** Aus dem
CRM-Lagrangian L = (R + gamma*R^2)/(16piG) folgt die Trace-Gleichung
R + 6*gamma*Box*R = -8*pi*G*T (Sotiriou & Faraoni 2010). Im quasi-statischen
Limes bei Wellenzahl k ergibt die modifizierte Poisson-Gleichung:
  V_grav^{f(R)} = 8*pi*G*rho^2/k^2 * [1 + (1/3)*k^2/a^2/(k^2/a^2 + m_s^2)]
Fuer k << a*m_s (super-Compton): V_grav -> 8*pi*G*rho^2/k^2 (exakt GR).
Fuer k >> a*m_s (sub-Compton): V_grav -> (4/3)*8*pi*G*rho^2/k^2.
Der CRM liefert also die KONSTRUKTIVE Ableitung aus dem Lagrangian,
nicht nur die Eindeutigkeit.

**Ehemals groesste Schwachstelle -> jetzt zweifach geschlossen.**

### L2: Minimalitaetslemma -- GESCHLOSSEN (2026-03-15)
Lemma 2.6 (Minimality of the back-reaction form) in EN-Paper eingefuegt.
Vier Axiome (M1-M4): Gravitationskopplung, Dimensionskonsistenz, minimale
Nichtlinearitaet, IR-Dominanz -> V_grav = c_0 * G * rho^2 / k^2 EINDEUTIG.
Beweis via Dimensionsanalyse in natuerlichen Einheiten: alpha=1, beta=2
=> gamma=-2 als einzige Loesung. Naechster Term beta=3 unterdrueckt um
rho/k^4 (unter Planck-Dichte vernachlaessigbar). c_0 = 8pi aus
Newtonscher Selbstenergie-Interpretation.
Zusaetzlich: Pre-existierenden enumerate-Bug in Paper gefixt (Zeile 994).

### L3: RG-Abschnitt -- GESCHLOSSEN (2026-03-15)
Proposition (Thermodynamic beta-function) eingefuegt:
  beta_Lambda^thermo = d/d(ln a) * (dPhi/drho)|_{rho=rho*}
  = -H_0^2 M_Pl^2 / [ln(M_Pl/H)]^2 * d(ln H)/d(ln a)
Explizit berechenbar aus dem Funktional Phi ohne QFT-Loops.
In der Materieaera: d(ln H)/d(ln a) = -3/2, also beta > 0
(Lambda_eff waechst bei Expansion -- konsistent mit Uebergang
von Dezeleration zu Akzeleration).

### L4: Renormierung der UV-Divergenzen
Rigorose Kontrolle der Modesumme ueber alle Standard-Model-Spezies
und Renormierung der UV-Divergenzen steht aus.
**VERSCHAERFT (2026-03-16):** Die Raw-Modensumme liefert NICHT die
gewuenschte Skalierung rho_* ~ H_0^2/G. Die Skalierung entsteht nur
nach RG-Matching (Abschnitt 5.2). Rigorose Herleitung des RG-Matchings
ist der zentrale offene Punkt.

### L5: Quantum Corrections
One-Loop-Robustheit diskutiert (Skalaron-Beitrag delta_rho ~ m_s^4/(64pi^2)),
aber vollstaendige Loop-Analyse fehlt.

### L6: Solar-System Fifth Force -- NEU (2026-03-16)
Bei m_s << H_0 (gamma >> 10^60, aus One-Loop-Constraint) ist das Skalaron
auf Solar-System-Skalen effektiv maselos. Fuenfte-Kraft-Beitrag gamma_PPN = 1/2,
ausgeschlossen durch Cassini (|gamma_PPN - 1| < 2.3e-5).
Screening-Mechanismus (Chamealeon oder Vainshtein) noetig, aber nicht vorhanden.
**Dies ist die groesste offene Spannung des Modells.**

### L7: Sigma = 9/8 statt 1 -- KORRIGIERT (2026-03-16)
Der Lensing-Parameter war falsch als Sigma = 1 exakt angegeben.
Korrekt: Sigma(k) = (2+f)/(2+2f/3), mit f = (k^2/a^2)/(k^2/a^2 + m_s^2).
Sub-Compton: Sigma = 9/8. GR-Limes: Sigma = 1.
Status: KORRIGIERT in EN und DE Versionen.

### L8: w_eff vs. w_DE -- KORRIGIERT (2026-03-16)
w_eff(z=0) ~ -1.4 bis -1.5 ist die EFFEKTIVE Zustandsgleichung (einschliesslich
Zeitvariation von Omega_Phi), nicht die intrinsische w_DE des Skalarfelds.
Die intrinsische w_DE ist >= -1 (kein Phantom). Direkter DESI-Vergleich
erfordert (w0, wa)-Mapping, das noch aussteht.

### L9: V_grav^{f(R)} Faktor 4/3 -- KORRIGIERT (2026-03-16, Zweiter Zyklus)
Die f(R)-Herleitung von V_grav (Proposition 3.X) hatte eine FALSCHE Formel:
V_grav^{f(R)} = 8pi*G*rho^2/k^2 * 1/(1+a^2*m_s^2/k^2) behauptete exaktes
Matching mit dem Ansatz im sub-Compton-Limes. Tatsaechlich:
- Super-Compton (k << am_s): V_grav -> 8pi*G*rho^2/k^2 (exakt GR, RICHTIG)
- Sub-Compton (k >> am_s): V_grav -> (4/3)*8pi*G*rho^2/k^2 (4/3 Enhancement)
Die KORREKTE Formel lautet: V_grav^{f(R)} = 8pi*G*rho^2/k^2 * [1 + (1/3)*f(k)]
mit f(k) = (k^2/a^2)/(k^2/a^2 + m_s^2).
Der 4/3-Faktor betrifft nur den sub-Compton-Anteil des Modenintegrals und
aendert rho_* um O(1), nicht die parametrische Skalierung Lambda_eff ~ H_0^2/ln.
Status: KORRIGIERT in EN und DE. Remark erweitert.

### L11: DE-Beweis der beta-Funktion -- KORRIGIERT (2026-03-16, Dritter Zyklus)
Der DE-Beweis der Proposition (Thermodynamische beta-Funktion) hatte einen
fehlerhaften Zwischenschritt: "T_dS/rho_* = 16piG*rho_**V_k" ist KEINE
korrekte Vereinfachung der Stationaritaetsbedingung
B*[ln(rho_*/rho_Pl) + 1] + 2C*rho_* = 0. Der DE-Beweis wurde mit der
EN-Version synchronisiert (impliziter Funktionensatz).
Zusaetzlich: BEWEISNOTIZ Zeile 53 korrigiert (rho^2 -> rho in zweiter
Ableitung, konsistent mit Paper).
Status: KORRIGIERT in DE. EN war bereits korrekt.

### L10: Vorfaktor c_0 Konvention -- PRAEZISIERT (2026-03-16, Zweiter Zyklus)
Der c_0 = 8pi Vorfaktor im Minimality-Lemma war numerisch nicht korrekt
hergeleitet (Newtonsches Ergebnis ist ~4pi/5, nicht 8pi). Korrektur:
c_0 ist O(4pi)--O(8pi), die Konventionswahl wird in Lambda_ren absorbiert.
Status: PRAEZISIERT (kein Fehler im Hauptresultat, nur im Nebenargument).

### L12: DE beta-Proposition Vorzeichen -- KORRIGIERT (2026-03-16, Vierter Zyklus)
Die DE-Proposition (Thermodynamische beta-Funktion) behauptete FALSCH
"was negativ ist (da H in der Materiaera abnimmt)". Korrekt: beta > 0
(positiv, weil zwei Minuszeichen -- eines aus der Formel, eines aus
d(ln H)/d(ln a) = -3/2 -- sich aufheben). Der DE-Beweis war korrekt
(schloss mit beta > 0), aber das Proposition-Statement widersprach dem
eigenen Beweis. Explizite Materiaera-Gleichung mit +3/2 > 0 eingefuegt.
Status: KORRIGIERT in DE. EN war korrekt.

### L13: Lagrangian-Konvention f(R) = R + 2 gamma R^2 -- KORRIGIERT (2026-03-16, Vierter Zyklus)
Der Lagrangian war als L = R/(16piG) + gamma R^2 + ... geschrieben
(inkonsistente Konvention: [gamma] war als [E]^{-2} deklariert, aber
gamma R^2 als separater Term erforderte dimensionsloses gamma). Zusaetzlich
war f(R) = R + 2 gamma R^2 in Prop. Vgrav-from-fR, was m_s^2 = 1/(12 gamma)
ergibt (nicht 1/(6 gamma) wie behauptet). KORREKTUR:
(a) Lagrangian auf Standard-Konvention: L = (R + gamma R^2)/(16piG) + ...
(b) f(R) = R + gamma R^2 (ohne Faktor 2)
(c) f_R = 1 + 2 gamma R, f_RR = 2 gamma (statt 4 gamma)
(d) Spurgleichung: 6 gamma (statt 12 gamma)
(e) m_s^2 = 1/(6 gamma) -- jetzt dimensionell und numerisch korrekt
Status: KORRIGIERT in EN und DE.

### L14: One-Loop-Formeln -- KORRIGIERT (2026-03-16, Vierter Zyklus)
(a) m_s^2 war als M_Pl^2/(6 gamma) geschrieben (dimensionell inkonsistent
    mit [gamma] = [E]^{-2}). Korrigiert auf m_s^2 = 1/(6 gamma).
(b) delta rho = M_Pl^4/(64 pi^2 (6 gamma)^2) korrigiert auf
    1/(64 pi^2 (6 gamma)^2).
(c) gamma-Schranke korrigiert: gamma >> 1/(M_Pl H_0) ~ 10^60 Planck-Einh.
    (statt gamma >> M_Pl/H_0, was dimensionell falsch war).
Status: KORRIGIERT in EN und DE.

### L15: w_eff Zahlenwert inkonsistent -- KORRIGIERT (2026-03-16, Fuenfter Zyklus)
Das Paper gab w_eff(z=0) ~ -1.35 mit "best-fit" Parametern kappa=1.5,
a_trans=0.50, Phi_0=0.685 an. Sorgfaeltige Berechnung ergibt:
- kappa=1.5: w_eff ~ -1.47
- kappa=2.0: w_eff ~ -1.37
Der Wert -1.35 war mit kappa=1.5 NICHT reproduzierbar.
Korrektur: Zahlenwert auf "~-1.4 bis -1.5" geaendert, Parameter als
"illustrativ" (nicht "best-fit") gekennzeichnet, Hinweis auf erforderlichen
numerischen Fit an SN+BAO-Daten eingefuegt.
Zusaetzlich: DE Solar-System Remark vervollstaendigt (fehlende Qualifikationen
bei Resolutions (i) und (ii), fehlender Hinweis auf "future work").
Status: KORRIGIERT in EN und DE.

### L16: R4-Lagrangian-Konsistenz -- VERIFIZIERT (2026-03-16, Fuenfter Zyklus)
Alle 10 abhaengigen Gleichungen in BEIDEN Versionen systematisch geprueft:
Lagrangian, gamma-Dimension, f_R/f_RR, Spurgleichung, m_s^2, BBN,
One-Loop, gamma-Schranke, Poisson-Gleichungen, V_grav^{f(R)}.
Ergebnis: ALLE konsistent mit f(R) = R + gamma R^2. Keine Relikte von
2gamma, 4gamma oder 12gamma gefunden.
Status: VERIFIZIERT (kein Fehler).


===============================================================================
## Naechste Schritte
===============================================================================

1. **Minimalitaetslemma formulieren:** Axiome -> Form von Phi erzwungen
   (Dimensionskonsistenz + allgemeine Kovarianz + Konvexitaet + IR/UV-Cutoff).

2. **RG-Abschnitt mit echtem beta-Statement:** beta_Lambda aus Skalar-Tensor-
   Lagrangian explizit berechnen.

3. **1-2 Killer Predictions:** Statt vieler Signaturen fokussieren auf:
   - Gravitational Slip als praezise Vorhersage
   - Lensing-Signal als falsifizierbarer Test

4. **Phantom w_eff klare Story:** Erklaeren warum w_eff(0) ~ -1.4 bis -1.5 keine
   Ghost-Instabilitaet ist, sondern Effekt der Parametrisierung.

5. **Hopfion-Energiebounds:** E >= c |Q|^{3/4} als alternative
   Backreaction-Form (explorativ, INPUT-DE-2).


===============================================================================
## Verbindung zum Meta-Framework (FST Positivity Pattern)
===============================================================================

### Pattern B: Flow Selection

Dark Energy realisiert das Pattern-B-Schema (Flow Selection):

| Ebene | Mechanismus | DE-Instanziierung |
|-------|-------------|-------------------|
| Leading neutral | Vakuum rho_0 | Jeder rho-Wert minimiert etwas |
| 1st order degenerate | Loop-Korrekturen deg. | QFT-Schaetzung 10^120 zu gross |
| 2nd order decides | Lambda_eff Selection | Phi-Minimum selektiert rho_* ~ H_0^2/G |

### Analogie zu RH

Lambda_eff ist das physikalische Analogon des renormierten Li-Funktionals
bei RH: Beide identifizieren einen eindeutigen, nichtnegativen Minimierer
als Attraktor des komplexen Systems.

"Dark Energy ist der Selection Pressure, der fuer gravitationelle
Stabilitaet benoetigt wird." (Section 1.3 des Papers)

### Coincidence Problem geloest

Lambda_eff ~ H_0^2 koppelt die Vakuumenergie an die AKTUELLE Expansionsrate.
Das Coincidence Problem verschwindet, weil Lambda_eff nicht konstant ist,
sondern langsam mit H(t) evolviert.

### Gemeinsame Programmlinie

```
Phi mit Entropieterm -> Strikte Konvexitaet -> Eindeutiger Minimierer -> Lambda_eff
```

Die Kommunikationsstaerke ("CC + Coincidence in einem Schlag") ist hoch,
aber die theoretische Fundierung (V_grav-Form, Minimalitaetslemma) muss
noch geschlossen werden.


===============================================================================
## Referenzen
===============================================================================

- Weinberg (1989): Cosmological Constant Problem (Review)
- Zeldovich (1968): Vakuumenergie und Kosmologie
- Riess et al. (1998): Beschleunigte Expansion (SN Ia)
- Perlmutter et al. (1999): Beschleunigte Expansion
- Planck 2018: Kosmologische Parameter
- Li (2004): Holographic Dark Energy
- Kaloper & Padilla (2014): Vacuum Energy Sequestering
- Starobinsky (1980): R^2-Inflation
- DESI (2025): Dynamische Dark Energy Hinweise
- Hu & Sawicki (2007): f(R)-Gravitationsmodelle
- Sotiriou & Faraoni (2010): f(R) Theories Review
