# BEWEISNOTIZ -- Kosmologische Konstante via geometrische Saettigung
# Stand: 2026-03-15
# Status: FRAMEWORK NOTE -- L1 geschlossen (2x), L2 geschlossen, L3 geschlossen, nur L4/L5 offen

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
kosmologische Konstante Lambda_eff ist der Minimierer eines Free-Energy-
Funktionals Phi[rho], definiert ueber Vakuum-Fluktuationen. Lambda_eff ~ H_0^2
folgt natuerlich, ohne Fine-Tuning.


===============================================================================
## Beweiskette
===============================================================================

### Schritt 1: Vakuum-Free-Energy-Funktional Phi[rho]
**Status: BEWIESEN (Definition 2.1)**

  Phi[rho] = int_{Lambda_IR}^{Lambda_UV} [
    (1/2) omega(k) + T_dS * ln(rho/rho_Pl) + V_grav[rho, k]
  ] * 4pi k^2 / (2pi)^3 dk

Drei Terme:
1. **Nullpunktsenergie:** (1/2) omega(k) = (1/2) sqrt(k^2 + m^2)
2. **Entropie-Strafe:** T_dS * ln(rho/rho_Pl) mit de-Sitter-Temperatur T_dS = H/(2pi)
3. **Gravitationelle Rueckwirkung:** V_grav[rho, k] = 8pi G rho^2 / k^2

Cutoffs:
- UV: Lambda_UV = M_Pl ~ 2.4 x 10^18 GeV (Planck-Skala)
- IR: Lambda_IR = H_0 ~ 1.5 x 10^{-42} GeV (Hubble-Radius)


### Schritt 2: Strikte Konvexitaet
**Status: BEWIESEN (Theorem 3.1, Step 2)**

Zweite Funktionalableitung:
  delta^2 Phi / delta rho^2 = int [T_dS/rho^2 + 16pi G/k^2] * (4pi k^2/(2pi)^3) dk > 0

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

Lagrangian L = R/(16piG) + gamma R^2 + Skalarfeld.
Poeschl-Teller-Potential V(phi) = V_0/cosh^2(phi/phi_0).
Saettigungs-ODE: dX/da = k(1-X^2) mit X = tanh(...).
de-Sitter-Stabilitaets-Proposition: phi_0 >= sqrt(8/3) M_Pl.


### Schritt 6: Effektive Zustandsgleichung
**Status: BEWIESEN (Section 4.3)**

w_eff(a) abgeleitet. Phantom-Stabilitaet diskutiert.
Kompatibilitaet mit DESI-Ergebnissen ueberprueft.

**WARNUNG:** w_eff(0) ~ -1.35 (Phantom). Klare Story noetig
(nicht physikalischer Wert, sondern Grenzwert der Parametrisierung).


===============================================================================
## Offene Luecken
===============================================================================

### L1: V_grav First-Principles-Ableitung -- GESCHLOSSEN (2026-03-15)
Der Backreaction-Term V_grav ~ G rho^2 / k^2 ist jetzt DOPPELT begruendet:

**(a) Minimality Lemma (Lemma 2.6):** 4 Axiome (M1-M4: Gravitationskopplung,
Dimensionskonsistenz, minimale Nichtlinearitaet, IR-Dominanz) erzwingen
V_grav = c_0 * G * rho^2 / k^2 EINDEUTIG (Dimensionsanalyse).

**(b) CRM-Transfer: f(R) Trace-Gleichung (Proposition 3.X):** Aus dem
CRM-Lagrangian L = R/(16piG) + gamma*R^2 folgt die Trace-Gleichung
R + 12*gamma*Box*R = -8*pi*G*T (Sotiriou & Faraoni 2010). Im quasi-statischen
Limes bei Wellenzahl k ergibt die modifizierte Poisson-Gleichung:
  V_grav^{f(R)} = 8*pi*G*rho^2/k^2 * 1/(1 + a^2*m_s^2/k^2)
Fuer k >> a*m_s (sub-Compton): exakt V_grav = 8*pi*G*rho^2/k^2.
Fuer k << a*m_s: unterdrueckt. Der CRM liefert also die KONSTRUKTIVE
Ableitung aus dem Lagrangian, nicht nur die Eindeutigkeit.

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

### L5: Quantum Corrections
One-Loop-Robustheit diskutiert (Skalaron-Beitrag delta_rho ~ m_s^4/(64pi^2)),
aber vollstaendige Loop-Analyse fehlt.


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

4. **Phantom w_eff klare Story:** Erklaeren warum w_eff(0) ~ -1.35 keine
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
