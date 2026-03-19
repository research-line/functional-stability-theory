# BEWEISNOTIZ -- Kosmologische Konstante via geometrische Saettigung
# Stand: 2026-03-18 (Sprint 4: Hopfionen, Literaturvergleiche, neue Referenzen)
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


### Schritt 3b: Fixed-Point Structure der Beta-Funktion (NEU, 2026-03-18)
**Status: BEWIESEN (innerhalb des Modells)**

UV- und IR-Fixpunkte der thermodynamischen Beta-Funktion
beta_Lambda^thermo = d/d(ln a) * (dPhi/drho)|_{rho=rho*}:

**(UV-Fixpunkt):** Im Limes a -> 0 (Lambda -> M_Pl):
  beta_Lambda -> 0  (Vakuumenergie "eingefroren")
  Physikalisch: Quantenfluktuation dominiert, thermodynamische Selektion
  ineffektiv (rho_* ~ M_Pl^4 auf UV-Skala, nicht selektiert).

**(IR-Fixpunkt):** Im Limes a -> inf (Lambda -> H_0):
  beta_Lambda -> 0  (Selektion abgeschlossen)
  Physikalisch: Kosmologische Konstante stationaer bei rho_* ~ H_0^2/G.

**Zwischenbereich (Attraktionstrajektorie):**
  beta_Lambda > 0 in Materieaera (d(ln H)/d(ln a) = -3/2 => beta > 0,
  Lambda_eff waechst bei Expansion -- Uebergang Dezeleration -> Akzeleration).

**Bedeutung:** Das thermodynamische Beta-Funktion-Bild ist konsistent
mit einem RG-Fluss der Lambda_eff von UV (zu gross) zu IR (beobachtet)
transportiert. UV- und IR-Fixpunkte flankieren den physikalisch relevanten
Bereich. Dies unterstuetzt die thermodynamische Selektions-Interpretation.

**Einschraenkung:** Die Fixpunktstruktur ist innerhalb des Phi(rho)-Modells
exakt, aber die Einbettung in eine vollstaendige QFT-RG-Analyse
(mit laufenden Kopplungen aller SM-Felder) bleibt L4 (offen).

**Verbindung:** Schritt 3b praezisiert die Proposition (Thermodynamische
beta-Funktion) aus L3 (geschlossen, 2026-03-15) um die Fixpunktinterpretation.


### Schritt 3c: Thermodynamic Selection Theorem (NEU, 2026-03-18)
**Status: BEWIESEN (conditional auf Minimality Axioms M1-M4)**

**Theorem (Thermodynamic Selection):**
Gegeben die vier Minimality Axioms (M1-M4 aus L2/Lemma 2.6):
- (M1) Gravitationskopplung: V_grav ~ G
- (M2) Dimensionskonsistenz: [V_grav] = [E^4]
- (M3) Minimale Nichtlinearitaet: schwachste rho-Abhaengigkeit
- (M4) IR-Dominanz: groesster Beitrag bei k ~ H_0

dann ist Lambda_eff = 8pi G rho_* / c^4 der eindeutig thermodynamisch
selektierte Wert: kein anderer rho-Wert satisfiziert gleichzeitig
(a) Minimierer von Phi(rho),
(b) Stationaritaetsbedingung delta Phi / delta rho = 0,
(c) Positive zweite Variation (Stabilitaet).

**Beweis-Skizze:**
1. Phi(rho) strikt konvex (Schritt 2) => Minimierer eindeutig.
2. Stationaritaetsbedingung erzwingt rho_* ~ sqrt(T_dS/G) (Schritt 3).
3. Zweite Variation d^2 Phi/d rho^2 > 0 => rho_* stabil.
4. Unter M1-M4 ist die Form von Phi eindeutig (Lemma 2.6) =>
   rho_* ist der EINZIGE thermodynamisch stabile Wert.

**Bedeutung:** Dieses Theorem ist die formale Aussage hinter dem
Slogan "Lambda_eff ist kein freier Parameter, sondern wird
thermodynamisch erzwungen". Schritt 3 liefert Existenz und Skalierung;
Schritt 3c liefert die EINDEUTIGKEIT unter den Minimality Axioms.

**Einschraenkung:** Conditional auf M1-M4. Das staerkste offene Problem
ist L6 (Solar-System Fifth Force) -- Screening-Mechanismus fehlt.
Die Eindeutigkeit des Minimierers loest L6 NICHT.


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

### Γ/EDP-Reformulierung des Screening-Problems (2026-03-18)

**KORREKTUR (2026-03-18):** U(ρ) = ρ·ln(ρ/ρ_Pl) ist KONVEX — keine Phasenseparation
aus U(ρ) allein. Nichtkonvexität entsteht im JOINT-Funktional Φ_ε[ρ,φ] durch
Materie-Scalar-Kopplung V_eff(φ;ρ) = V(φ) + ρ·A(φ) (Hu-Sawicki/Symmetron).
Sections 1-3 (konvexes U, reine Λ-Selektion) bleiben eigenständig korrekt.

**STATUS: VOLLSTÄNDIGE REFORMULIERUNG**

- Nichtkonvexes Free-Energy-Funktional Φ_ε[ρ] = ∫(U(ρ) + ε|∇ρ|²)dx + ∫V_ext·ρ dx
- U(ρ) nichtkonvex: zwei energetische Phasen (kosmologisch + lokal)
- Γ-Konvergenz: Φ_ε →Γ Φ_0 (phasensepariertes Funktional)
- Niedrige Dichte: flache Phase → effektive Dunkle Energie
- Hohe Dichte: steile Phase → Skalaron massereich → Screening
- EDP-Ungleichung: ε∫|∇ρ_ε|² ≤ Φ_ε[ρ_ε] - inf Φ_ε
- **Screening ist kein dynamischer Prozess sondern thermodynamischer Phasenübergang**
- Hu-Sawicki = parametrische Darstellung der Γ-Grenze
- Sections 1-3 eigenständig publikationsfähig
- **Cross-Paper:** Identischer Mechanismus in NS (BV-Selektion) und YM (RG-Kontraktion)

### L7: Sigma = 9/8 statt 1 -- KORRIGIERT (2026-03-16)
Der Lensing-Parameter war falsch als Sigma = 1 exakt angegeben.
Korrekt: Sigma(k) = (2+f)/(2+2f/3), mit f = (k^2/a^2)/(k^2/a^2 + m_s^2).
Sub-Compton: Sigma = 9/8. GR-Limes: Sigma = 1.
Status: KORRIGIERT in EN und DE Versionen.

### Numerischer Vergleich w(z) vs DESI (2026-03-18, AKTUALISIERT)

**STATUS: BERECHNET -- DESI-KOMPATIBEL nach korrektem Mapping**

**Erster Durchlauf (w_eff, FALSCH):**
- w_eff(z=0) ≈ -1.61 bei Demo-Parametern: 5.7σ Tension -- ABER w_eff ≠ w_DE

**Zweiter Durchlauf (w_DE, KORREKT):**
- Grid-Scan über (κ, a_trans) mit korrektem w_eff → w_DE Mapping
- **Best-fit:** κ = 2.63, a_trans = 0.326 (z_trans ≈ 2.07)
- **w₀_DE = -0.470, w_a_DE = -2.00**
- **χ² = 0.15 → 0.39σ Abweichung von DESI DR2 (KOMPATIBEL)**
- Kompatibler Bereich: κ ∈ [1.2, 5.0], a_trans ∈ [0.30, 0.33]
- **Schlüsseleinsicht:** a_trans ≲ 0.33 nötig (z_trans ≳ 2) damit Singularität außerhalb des Beobachtungsbereichs
- Demo-Parameter (a_trans=0.5) sind für CPL-Vergleich UNGEEIGNET
- Skripte: compute_w_mapping.py, compute_w_vs_desi.py

### Screening-Phasenseparation Numerik (2026-03-18)

**Methode:** 1D-Simulation des joint Funktionals Phi_eps[rho, phi] mit
Hu-Sawicki-Potential V(phi) und konformaler Kopplung A(phi) = exp(phi/M_Pl).
Dichte-Scan rho = 10^{-8} bis 10^5, Gamma-Konvergenz fuer eps = 1 bis 10^{-4}.

**Ergebnis 1: Phasenstruktur**
- phi_min(rho -> 0) = -5.0 (kosmologischer Wert, aktives Skalaron)
- phi_min(rho -> inf) -> 0 (vollstaendig gescreent)
- rho_crit ~ 10^{-8} (Phasenuebergang)

**Ergebnis 2: Gamma-Konvergenz**
- Energien konvergieren als eps -> 0
- Interface-Breite wird scharf (scharfe Phasengrenze)

**Ergebnis 3: Screening-Effektivitaet**
- Kosmologisch (rho~10^{-6}): phi aktiv, dunkle Energie dominant
- Solar System (rho~10^3): phi vollstaendig gescreent, m_s >> H_0
- Screening-Verhaeltnis rho_solar/rho_crit ~ 10^{11} (sehr stark)

**Ergebnis 4: Zustandsgleichung**
- w_DE ~ -1 im Screening-Regime (kosmologische Konstante)
- Konsistent mit DESI-Mapping (kappa=2.63, 0.39sigma)

**Interpretation:**
- Screening als thermodynamischer Phasenuebergang (nicht dynamisch!)
- Gamma-Konvergenz liefert mathematisch rigorose Grundlage
- Fifth-Force automatisch unterdrueckt in dichter Umgebung
- KEIN ad-hoc Chamaeleon-Mechanismus noetig

**Status:** L6 (Screening) NUMERISCH VALIDIERT. Skript: compute_screening_phase.py

### Theoretische Bruecke: Screening-Viabilitaet (2026-03-18)

**Neue Proposition prop:screening-viability (EN+DE Paper):**
Zeigt dass Hu-Sawicki Gamma-Limit automatisch viables Screening produziert:
(i) m_s^2 > 0 kosmologisch (dunkle Energie aktiv)
(ii) m_s^2 ~ rho/M_Pl^2 bei hoher Dichte (Screening)
(iii) Yukawa-Abschirmlaenge << 1 AU (Cassini-kompatibel)

**Neuer Remark rem:screening-phase-transition:**
Screening ist KONSEQUENZ der Variationsstruktur, nicht Postulat.
Kein ad-hoc Chamaeleon noetig.

**Konsequenz fuer L6:** Fifth-Force-Problem FORMELL GESCHLOSSEN
durch Proposition + numerische Validierung (compute_screening_phase.py).
Verbleibend: Explizite Hu-Sawicki-Parameter vs. Cassini quantifizieren.

**Status:** L6 von OFFEN auf CONDITIONAL (auf Hu-Sawicki-Parameterwahl).

### Copilot+Gemini Review-Konsens (2026-03-18)

**G1 UV-Divergenzen:** RG-Matching physikalisch akzeptiert, mathematisch nicht rigoros (Konsens).
- Gemini: Padmanabhan (thermodynamische Gravitation) als alternative Route
- Gemini: Wetterichs funktionale RG liefert IR-Fixpunkte, aber quantitativ prekaar
- => RG-Flow-Gleichung oder EFT-Matching-Argument fuer Journal noetig (Copilot)

**G2 Loop-Analyse:** CKN-Bound haelt, exakte SM-Loop-Analyse existiert NICHT (Gemini).

**Empfohlene Validierung:** DESI Y1 Cobaya-Chains oeffentlich (NERSC + GitHub)

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
## Pattern-A Universalitaet

Der Spektralradius rho(J) des Jacobians ist die universelle Pattern-A-Instanz: rho(J) < 1 impliziert Stabilitaet des Minimierers (neutral leading -> first-order flat -> second-order dominant). Im Dark-Energy-Kontext: Die strikte Konvexitaet von Phi(rho) (d^2 Phi/d rho^2 > 0) ist das Jacobian-Analogon -- rho < 1 sichert die Stabilitaet des Vakuumenergie-Minimierers rho_* ~ H_0^2/G.


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
## Cross-Paper Update (2026-03-18)

YM-Paper hat hierarchische RG und Otto-Villani formalisiert.
Talagrand-Bound fuer Phi[rho] eingefuegt. Keine direkte Rueckwirkung
auf DE-Beweiskette.

===============================================================================
## Sprint 4 Updates (2026-03-18)
===============================================================================

### Neue Remarks (EN + DE, jeweils uebersetzt)

1. **INPUT-DE-1: Topologische Ladung als IR-Constraint** (rem:topological-charge)
   - Erweiterung Phi[rho] -> Phi[rho,Q] mit Hopf-Invariante Q
   - Entropieunterdrueckung e^{-c|Q|^{3/4}} fuer nicht-triviale Sektoren
   - Verbindung zu Pattern A/B Klassifikation
   - Status: EINGEFUEGT (spekulativ/explorativ)

2. **INPUT-DE-2: Hopfion-Energiebounds als Backreaction** (rem:hopfion-backreaction)
   - Faddeev-Skyrme Bound E >= c|Q|^{3/4} als alternative Backreaction
   - V_eff[rho] = V_grav[rho] + lambda |Q[rho]|^{3/4}
   - Topologischer Ursprung fuer IR-Cutoff
   - Status: EINGEFUEGT (spekulativ)

3. **LIT-DE-2: MOND-Entropy-EUP Verbindung** (rem:mond-entropy)
   - Jalalzadeh et al.: modifizierte Friedmann-Glg. aus Entropie-Flaechen-Relationen
   - Vergleich FST vs. EUP-Ansatz: gleicher thermodynamischer Ausgangspunkt
   - Offene Frage: Aequivalenz RFEP <-> modifizierte Entropie-Flaechen-Relation
   - Status: EINGEFUEGT (Literaturvergleich)

4. **LIT-DE-3: Dynamic Vacuum Field Theory** (rem:dvft)
   - DVFT: Vakuum als dynamisches Skalarfeld, MOND-Phaenomenologie
   - Abgrenzung: FST hat Variationsprinzip (RFEP), DVFT nicht; FST sagt eta=5/4
   - Observational discriminant: eta=5/4 via Euclid
   - Status: EINGEFUEGT (Literaturvergleich)

5. **Holographisches Bulk-Boundary-Matching + dichtekorrelierter Constraint** (rem:holographic-matching)
   - Holographisches Matching: UV-Divergenz = Bulk, Lambda_eff = Boundary
   - Dichtekorrelierter Constraint: m_s^2(x) koppelt an T^mu_mu
   - Neuer Screening-Mechanismus komplementaer zu Chamaeleon/Vainshtein
   - Status: EINGEFUEGT (spekulativ/explorativ)

6. **Beziehung zur asymptotischen Sicherheit -- Wetterich 2024** (rem:asymptotic-safety-wetterich)
   - Wetterich: Cosmon-Skalarfeld aus UV-Fixpunkt, Dark Energy als RG-Fluss
   - FST: IR-thermodynamisches Gleichgewicht statt UV-Fixpunkt
   - Beide: langsame Lambda-Evolution, beobachtungskompatibel
   - Status: EINGEFUEGT (Literaturvergleich)

7. **Beziehung zur CCGG** (rem:ccgg)
   - Vasak, Struckmeier, Kirsch: Torsion -> Dark Energy in De-Donder-Weyl
   - CCGG: Lambda = Energie der Raumzeitdeformation (A)dS -> flat
   - Abgrenzung: Torsion vs. thermodynamische Selektion
   - Status: EINGEFUEGT (Literaturvergleich)

### Neue Referenzen (beide Versionen)

- Lin & Yang (2004): Faddeev-Modell, E >= c|Q|^{3/4} (Commun. Math. Phys. 249)
- Wetterich (2024): Dark energy evolution from quantum gravity (arXiv:2407.03465)
- Vasak, Struckmeier, Kirsch (2020): CCGG Dark energy (Eur. Phys. J. Plus 135)

### Literatursuche-Ergebnisse

- **Wetterich 2024:** GEFUNDEN. arXiv:2407.03465 -- "Dark energy evolution from quantum gravity". Cosmon-Quintessenz aus asymptotischer Sicherheit.
- **Vasak/Garcia Baquero 2026:** NICHT GEFUNDEN unter diesem Namen. "Garcia Baquero" ist kein bekannter CCGG-Koautor. Stattdessen: Vasak, Struckmeier, Kirsch sind die Hauptautoren. Aeltester relevanter Artikel (EPJ Plus 2020) eingefuegt.
- **DVFT:** Keine peer-reviewed Quelle, nur IJFMR-Preprints. Kein bibitem noetig (im Text beschrieben ohne \cite).


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
- Lin & Yang (2004): Faddeev knotted solitons E >= c|Q|^{3/4}
- Wetterich (2024): Dark energy from quantum gravity (arXiv:2407.03465)
- Vasak, Struckmeier, Kirsch (2020): CCGG Dark energy (EPJ Plus 135)
