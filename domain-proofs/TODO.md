# TODO -- Folgebeweise: Beweisketten-Status

Stand: 2026-03-15 (Update: TU JOURNAL-GRADE (DFC-Sprint), NS BV-Refactor+AGC, DE L1-L3 geschlossen, Framework No-Go)

> Siehe: `PUBLIKATIONSSTRATEGIE.md` in jedem Unterordner fuer Details.

---

## Review-Gesamturteil (2026-03-10)

> "Keines der Programme behauptet mehr implizit einen Beweis, wo keiner existiert.
> Das ist der entscheidende qualitative Sprung."

### Priorisierte Publikationsreihenfolge (Reviewer-Empfehlung):
1. **Turbulenz** -- "Clean Win", am geschlossensten
2. **Navier-Stokes** -- Stark als conditional framework (NICHT als Beweis)
3. **Yang-Mills** -- Strong-coupling Gap (Reframing noetig)
4. **Dark Energy** -- Framework Note (Minimalitaetslemma GESCHLOSSEN 2026-03-15, RG-Kern offen)
5. **Unified Framework** -- 6/6 Review-TODOs geschlossen (2026-03-15), Essay/Preprint-ready

### Gemeinsame Programmlinie:
```
F mit Entropieterm -> Strikte Konvexitaet -> Eindeutiger Minimierer -> Physik
```
Zwei Stellen muessen praezise werden:
- **Axiomatik:** Welche minimalen Prinzipien erzwingen F?
- **Nichttrivialitaet:** Wo reicht "KL ist konvex" NICHT?

---

## 1. Navier-Stokes (~95% -> CONDITIONAL, nicht Beweis)

### Erledigt
- [x] Definition H(u|A) -- Attractor-Distanz sauber definiert
- [x] Lemma (Energy-Dissipation Bound) -- Voll bewiesen
- [x] Haupttheorem 3.1 -- 3-Schritt-Widerspruchsbeweis geschlossen
- [x] Measurable Selection -- Standard-Verweis (Kuratowski-Ryll-Nardzewski)
- [x] Discussion + References komplett
- [x] Lemma A (Biot-Savart) -- Vollstaendiger 4-Schritt-Beweis (2026-03-10)
- [x] Lemma B (Enstrophy Balance) -- Vollstaendiger Beweis mit Gronwall-Kette (2026-03-10)
- [x] Attractor-Regularitaet -- Proposition 2.2 mit Bootstrap (2026-03-10)
- [x] Notation vereinheitlicht auf nicht-fett (2026-03-10)

### Review-Kritik (2026-03-10)
- **KRITISCH:** Gronwall-Luecke bei Attractor-Minimierer-Zeitabhaengigkeit (nicht rigoros geschlossen)
- **REFRAMING:** Als "structural reduction" positionieren, NICHT als Beweis
- **TITEL AENDERN:** -> "A Thermodynamic Obstruction to Finite-Time Blow-Up (Conditional Framework)"

### Offen (Review-bedingt)
- [x] **Titel aendern** -- "A Thermodynamic Obstruction to Finite-Time Blow-Up (Conditional Framework)" (2026-03-11)
- [x] **Gronwall-Luecke** explizit als offene Annahme kennzeichnen -- Assumption G (G1+G2) eingefuehrt, Remark zum Status, Limitations erweitert (2026-03-11)
- [x] **Abstract** -- komplett umgeschrieben als "conditional regularity criterion" (2026-03-11)
- [x] **LaTeX kompilieren und PDF pruefen** -- 13 Seiten, keine Fehler (2026-03-11)
- [x] **Letzter Review:** Konsistenz geprueft (2026-03-11): 3 Fixes -- (1) `remark`-Umgebung im Praeambel definiert, (2) hardcodierte Sektionsnummer "Theorem 3.1" -> "Main Result", (3) Widerspruch in Discussion behoben ("without external hypothesis" -> mit Assumption-G-Verweis). PDF rekompiliert, 13 Seiten, keine Fehler.

### Integrierte Ideen aus Framework/Anwendungsfaelle (2026-03-10)
- [x] ~~**Minimax-Reformulierung (Discussion):**~~ Integriert als Subsection "A Game-Theoretic Perspective" -- Regularitaet als Nullsummenspiel, Legendre-Fenchel S >= F.
- [x] ~~**MFG-Verbindung (Remark):**~~ Integriert als Subsection "Connection to Mean Field Game Theory" + neue Referenz Lasry & Lions 2007.
- [x] ~~**Legendre-Fenchel-Perspektive:**~~ In Game-Theoretic Perspective integriert.

### Integrierte Ideen aus CRM (2026-03-10)
- [x] ~~**Saettigungsdynamik (Remark):**~~ Integriert als Remark nach Theorem 3.1 -- Blow-up-Verhinderung als Saettigungsmechanismus, Analogie zu Poeschl-Teller-Potential.

### Integrierte Copilot-Antworten (2026-03-10)
- [x] ~~**Kumulantenentwicklung fuer Enstrophie:**~~ Integriert als Remark "Cumulant viewpoint on enstrophy fluctuations" -- Sub-Gaussianitaet, Kumulantenschranke |kappa_n| <= n! sigma_H^2/2 alpha^{-(n-2)}, Konzentration um Attraktor. Neue Referenzen: Kuksin-Shirikyan 2012, Hairer-Mattingly 2006.

---

## 2. Yang-Mills (~90% -> REFRAMING NOETIG)

### Review-Kritik (2026-03-10)
- **KRITISCH:** KL-Konvexitaet ist trivial -- gibt keinen Mass Gap (jedes Gibbsmass hat das)
- **KRITISCH:** Hessian = L^2-Kruemmung, aber Mass Gap = Transferoperator-Spektrum -- Bruecke fehlt
- **REFRAMING:** "curvature bound <=> Poincare inequality <=> spectral gap (strong coupling)"
- **Motor:** Holley-Stroock + Dobrushin-Zegarlinski, NICHT Free-Energy-Konvexitaet
- **Continuum-Teil:** Radikal kuerzen oder als "conditional programme" labeln

### Erledigt
- [x] Wilson-Gitter-Regularisierung sauber definiert
- [x] Free-Energy Functional F[mu] definiert mit KL-Zerlegung
- [x] Strikte Konvexitaet -- Voll bewiesen (Standard-Argument)
- [x] Hessian/Second Variation definiert
- [x] Transfer-Operator definiert, Spectral Gap = Mass Gap korrekt formuliert
- [x] Discussion mit explizit benannten offenen Annahmen
- [x] 15 Referenzen (Holley-Stroock 1987, Zegarlinski 1992)

### Integrierte Copilot-Antworten (2026-03-10)
- [x] ~~**Theorem 4.1 Step 1 (Poincare):**~~ Holley-Stroock Bounded Perturbation statt Log-Konkavitaet. kappa_link >= kappa_Haar * e^{-2 d_l beta}.
- [x] ~~**Theorem 4.1 Step 2 (Cluster Expansion):**~~ Dobrushin-Zegarlinski Kriterium. Uniforme Poincare kappa(Lambda) >= (1-eta) kappa_link fuer beta < beta_0.
- [x] ~~**Theorem 4.1 Step 3 (Transfer -> Spectral Gap):**~~ Vollstaendiger Beweis: Transferkern-Konstruktion, reversibeler Markov-Operator P auf L^2(pi), Poincare <=> Spektralluecke, exponentielles Clustering. Delta_Lambda = -log(1-kappa) >= kappa.
- [x] ~~**Reflection Positivity (Lemma 3.3):**~~ Vollstaendiger 3-Schritt-Beweis: (1) Reduktion auf Crossing-Kern-Positivitaet, (2) Einzelne Plaquette via Peter-Weyl + nichtnegative Charakterkoeffizienten, (3) Schur-Produktsatz fuer Produkt positiver Kerne. Remark zu Wilson vs. Heat-Kernel Action.
- [x] ~~**Theorem 4.2 (Kontinuumslimit):**~~ Ehrliche konditionale Formulierung mit 3 Annahmen (CL1-CL3). Remark zu Millennium-Problem.
- [x] ~~**Uniformitaet der Poincare-Konstante:**~~ Bestaetigung, dass Analytizitaets-Argument NICHT funktioniert. Bereits in Remark rem:all-beta korrekt adressiert.
- [x] ~~**Limitations aktualisiert:**~~ 3 klare offene Probleme.

### Offen (Polieren)
- [x] **LaTeX kompilieren und PDF pruefen** -- 11 Seiten, keine Fehler (2026-03-11)
- [x] **Letzter Review:** Konsistenz geprueft (2026-03-11): Alle Cross-Referenzen (lem:reflection-positivity, thm:lattice-gap, thm:continuum-gap, rem:all-beta), Konstanten (kappa, beta_0, d_l=6, eta) und Nummerierungen korrekt. Keine Aenderungen noetig.

### Integrierte Ideen aus Framework/Anwendungsfaelle (2026-03-10)
- [x] ~~**Funktionale Positivitaet unter Constraints:**~~ Integriert als Subsection 5.5 "The Mass Gap as Constrained Functional Positivity" -- Weil-Analogie, Positivitaet auf physikalischer Testklasse.
- [x] ~~**Vakuumlandschaft analog Proteinfaltung:**~~ Integriert als Remark "Analogy to protein folding landscapes" -- Gribov-Kopien vs. misfolded states, Jacobian-Analyse.

### Integrierte Ideen aus CRM (2026-03-10)
- [x] ~~**Confinement als Saettigung (Remark):**~~ Integriert als Remark -- Poeschl-Teller-Analogie fuer Confinement, Saettigungsdynamik im starken Feld.

### Offene Ideen
- [ ] **Stabilitaets-Selektions-Axiomatik (S1)-(S4) fuer YM:** Quelle: FST-GP Sec. 5.5. -> Neue Section oder Appendix.
- [ ] **Ordnungsparameter-Dualitaet (Kosten/Antwort):** Drei Beweisrouten (D1-D3). Quelle: FST-GP Sec. 4.3. -> In Discussion als "Alternative proof strategy".
- [ ] **Pfadintegral als Nash-Selektor:** Quelle: FST-I Sec. 3. -> Remark nach Theorem 4.2.
- [ ] **QCD Confinement als Nash-Phasenuebergang:** Quelle: FST-I Sec. 6.4. -> In Discussion ergaenzen.

---

## 3. Turbulenz (~70% -> CLEAN WIN, Prioritaet 1 fuer Publikation)

### Review-Kritik (2026-03-10)
- **STAERKE:** "Am naechsten an einem echten, in sich geschlossenen Resultat"
- **STAERKE:** Joint minimisation ueber Closures ist "genial"
- **UPGRADE A:** ~~(A2) von Axiom zu Theorem~~ **ERLEDIGT 2026-03-15**: Proposition 3.3 via Buckingham Pi, Axiom-Satz auf A1+A3 reduziert
- **UPGRADE B:** Flux-Constraint relaxieren (penalized/soft constraint), K41 als Grenzfall
- **UPGRADE C:** DNS-Test als reproduzierbares Protokoll (F[E] aus Spektren)

### Erledigt
- [x] Littlewood-Paley-Zerlegung sauber definiert
- [x] Energiefluss-Constraint definiert
- [x] Funktional F[E] (KL-Divergenz) sauber definiert
- [x] Intermittenz-Korrekturen motiviert (K62)
- [x] Discussion mit Limitationen
- [x] 10 Referenzen

### Offen (Beweis-Reparatur)
- [x] ~~**KRITISCH -- Theorem 1 Zirkularitaet behoben:**~~ (2026-03-10) Reformuliert mit Admissible Flux Family (Def. 3.1, drei Axiome). Joint-Minimierung ueber (E, Pi) in R^n_+ x A. 5-Schritt-Beweis. Neues Remark zu Replikator-Dynamik/Nash-GG.
- [x] ~~**KRITISCH -- Theorem 2 F-Monotonie reformuliert:**~~ (2026-03-10) Copilot-Analyse: KL-Konvexitaet unter mass-preserving maps gilt NICHT fuer NS-Transfer (trilinear, kein Markov-Operator). Loesung: Explizite "Assumption NL" (entropy-compatible cascade, Def. 3.x), regularisiertes F_epsilon, Flux-Repraesentations-Lemma (Duchon-Robert), 3-Faelle-Remark (Shell-Modelle, stochastische Kaskaden, NS mit Downhill-Flux).
- [x] ~~**Theorem 2 Step 4 (Quantitative Schranke):**~~ (2026-03-10) Copilot-Analyse: F-Monotonie allein liefert keine nu-unabhaengige untere Schranke. Kontrolliert nur Entropie-Dissipation D_F, nicht physikalische Dissipation epsilon_nu. Loesung: Explizite "Assumption ND" (nu-uniforme Energie-Injektion m_0 > 0). H^{-1} als kanonische Sobolev-Norm identifiziert. Schranke: epsilon_* >= c ||f||_{H^{-1}}^2 / M^2 unter (NL)+(ND).

### Integrierte Ideen aus Framework/Anwendungsfaelle (2026-03-10)
- [x] ~~**MEPP als Begruendung des laminar-turbulenten Uebergangs:**~~ Integriert als Discussion-Absatz "Entropy production and the laminar-turbulent transition".
- [x] ~~**Replikator-Gleichung / Gibbs-Minimierung:**~~ Integriert als Remark 3.5 in Theorem-1-Section (Nash-GG, Lyapunov-Funktion).

### Integrierte Ideen aus CRM (2026-03-10)
- [x] ~~**Skalenabhaengige Kopplung (Subsection):**~~ Integriert als neue Subsection "Scale-dependent coupling and the inertial range" -- Poeschl-Teller-Analogie fuer Energietransfer, Saettigungsmechanismus im Inertialbereich.

### Offene Ideen
- [ ] **Energiekaskade als Minimax-Spiel:** Kolmogorov-Spektrum als Sattelpunkt. Quelle: FST-GP Sec. 6.
- [ ] **Dissipative Adaptation (England):** Quelle: FST-III Sec. 2.3. -> Remark in Discussion.
- [ ] **Kumulantenentwicklung fuer Intermittenz:** Phi >= 0 schrankt Momentenhierarchie. Quelle: FST-GP Principle 2.3. -> Section 4.
- [ ] **Koalitionsbruch-Modell:** Quelle: FST-III Sec. 6. -> Remark in Discussion.

---

## 4. Dark Energy (~55% -> FRAMEWORK NOTE, braucht Upgrades)

### Review-Kritik (2026-03-10)
- **STAERKE:** "Kommunikativ extrem stark" (CC + Coincidence in einem Schlag)
- **SOLLBRUCHSTELLE:** ~~V_grav ~ G*rho^2/k^2 ohne First-Principles-Ableitung~~ **GESCHLOSSEN 2026-03-15**: Doppelt begruendet (Minimality Lemma + CRM f(R) Trace-Gleichung)
- **ROTES TUCH:** Phantom w_eff(0) ~ -1.35 -- klare Story noetig
- **UPGRADE A:** ~~Minimalitaetslemma (wenige Axiome -> Form erzwungen)~~ **ERLEDIGT 2026-03-15**: Lemma 2.6 mit 4 Axiomen (M1-M4) eingefuegt, Beweis via Dimensionsanalyse
- **UPGRADE B:** Sauberer RG-Abschnitt mit echtem Beta-Statement
- **UPGRADE C:** 1-2 Killer Predictions statt viele Signaturen

### Erledigt
- [x] Funktional Phi[rho] definiert (3 Terme)
- [x] Cutoffs physikalisch motiviert
- [x] Strikte Konvexitaet gezeigt
- [x] Numerisches Korollar: ~3.8 x 10^{-53} m^{-2}
- [x] Coincidence Problem qualitativ adressiert
- [x] 7 Referenzen (urspruenglich) -> 11 Referenzen nach CRM-Integration

### Erledigt durch CRM-Integration (2026-03-10)
- [x] **Back-Reaction-Term abgeleitet:** Neue Section 3.5 "Connection to scalar-tensor gravity" -- Lagrangian L = R/(16piG) + gamma*R^2 + Skalarfeld, Poeschl-Teller-Potential, Saettigungs-ODE, de-Sitter-Stabilitaets-Proposition.
- [x] **TODO (rot) Prefaktor g_*:** Behoben -- In Subsection 4.3 "Effective equation of state" integriert.
- [x] **TODO (rot) Lambda_eff(a), z_acc:** Behoben -- Subsection 4.3 mit w_eff(a), DESI-Kompatibilitaet.
- [x] **TODO (rot) epsilon vs. DESI/Euclid:** Behoben -- Subsection 5.5 "Falsifiable predictions" (gravitational slip, lensing, Oszillationen, MOND-Skala).
- [x] **BBN-Schutz:** Neue Subsection 4.4 "BBN protection via trace coupling".
- [x] Neue Referenzen: Starobinsky 1980, DESI 2025, Hu-Sawicki 2007, Sotiriou-Faraoni 2010.

### Offen (Theoretische Luecken)
- [x] **Renormierung:** Neue Subsection 5.3 "Renormalisation and the status of the scaling law" -- RG-Matching-Interpretation, nichtlokale IR-Logs, Reinterpretation von V_grav als effektiver nichtlokaler Term. Verbleibt: Explizite Berechnung von beta_Lambda aus dem Skalar-Tensor-Lagrangian. (2026-03-10, via Copilot)
- [x] **Selbstkonsistenz mit Friedmann:** Neue Subsection 3.6 "Self-consistent dynamics" -- Autonomes System (x,y,phi), exakte Reduktion auf Saettigungs-ODE, Gronwall-Fehlerschranke, Parameter-Identifizierbarkeit, de-Sitter-Stabilitaetsbedingung phi_0 >= sqrt(8/3) M_Pl. (2026-03-10, via Copilot)
- [x] **Quantum Corrections:** Neue Subsection 5.4 "One-loop robustness" -- Skalaron-Beitrag delta_rho ~ m_s^4/(64pi^2), Bedingung gamma >> O(M_Pl/H_0), Materie-Loops Planck-unterdrueckt, Naturalness-Diskussion. (2026-03-10, via Copilot)

### Integrierte Ideen aus Framework/Anwendungsfaelle (2026-03-10)
- [x] ~~**MEPP State Selection fuer Lambda:**~~ Integriert als Subsection 5.4 "Thermodynamic selection of the vacuum energy" + neue Referenz Martyushev & Seleznev 2006.

### Integrierte Ideen aus CRM (2026-03-10)
- [x] ~~**CRM-Instanziierung (Skalar-Tensor-Gravitation):**~~ Portiert als Section 3.5 -- Lagrangian, Poeschl-Teller, Saettigungs-ODE, de-Sitter-Stabilitaet.
- [x] ~~**tanh-Saettigung:**~~ In Saettigungs-ODE und Proposition integriert.
- [x] ~~**Effektive Zustandsgleichung + DESI:**~~ Subsection 4.3 mit w_eff, Phantom-Stabilitaet.
- [x] ~~**Falsifizierbare Vorhersagen:**~~ Subsection 5.5 mit 4 konkreten Tests.

### Offene Ideen
- [ ] **Kosmologische Phasenuebergangs-Kette:** Quelle: FST-I Sec. 6.4. -> Remark in Discussion.
- [ ] **Legendre-Fenchel fuer DE:** S >= F. Quelle: FST-GP Sec. 6 (Route B). -> In Discussion.
- [ ] **England-Ungleichung fuer Expansionsrate:** Quelle: FST-II Sec. 3.2. -> Explorativ, niedrige Prioritaet.

---

## Priorisierung (aktualisiert 2026-03-10, nach Review-Auswertung)

### Neue Reihenfolge (Reviewer-Empfehlung):
1. **Turbulenz** -- "Clean Win". Joint-Minimierung ist genuiner Beitrag. Upgrades A-C noetig.
2. **Navier-Stokes** -- Stark als conditional framework. Gronwall-Luecke offen benennen.
3. **Yang-Mills** -- Reframing auf Poincare/Dobrushin. Continuum-Teil kuerzen.
4. **Dark Energy** -- Framework Note. Minimalitaetslemma + RG-Kern fehlen.

### Gemeinsame Schwachstelle aller Papers:
> "Wo reicht 'KL ist konvex' NICHT? Wo wird echte Struktur gebraucht?"
- Turbulenz: Joint-Minimierung ist das Neue (nicht KL allein)
- NS: Attractor-Stabilitaet ist das Neue (nicht Konvexitaet allein)
- YM: Poincare/Dobrushin ist der Motor (nicht Free-Energy-Konvexitaet)
- DE: RG-Matching ist die Physik (nicht Funktional-Konvexitaet)

**Alle 12 Copilot-Anfragen sind geschlossen.** Keine offenen Anfragen mehr.

---

## 5. Externe Evaluierung -- Review-Inputs (2026-03-11)

Basierend auf umfassender externer Evaluierung des FST/CRM-Gesamtprogramms.
SotA-Abgleich, kritische Luecken, methodische Bruecken fuer alle Folgebeweise.

### 5.1 Navier-Stokes -- Neue Aufgaben

- [ ] **EXT-NS-1: Gronwall-Luecke via Lipschitz-Regularitaet schliessen**
  Entscheidend: Exakte Charakterisierung der Lipschitz-Regularitaet der Projektionen
  auf den Leray-Hopf-Attraktor. Zerlegung H = P_N + Q_N (Stokes-Operator).
  Null-Lipschitz-Abweichung (dev_m(A)=0) impliziert: hochfrequente Moden werden
  vollstaendig durch niederfrequente determiniert ("versklavt").
  -> Lit: Barbu/Cannone (2016) -- Log-Lipschitz-Regularitaet von Projektionen
  -> Lit: Eden, Foias, Nicolaenko -- Inertial Manifolds (klassisch)
  -> Schliesst die als "KRITISCH" markierte Gronwall-Luecke.

- [ ] **EXT-NS-2: Log-Dirichlet-Quotienten in KL-Joint-Minimierung integrieren**
  Die Log-Dirichlet-Quotienten fuer Differenzen von Loesungen auf dem Attraktor
  muessen in das FST-TU KL-Joint-Minimierungs-Verfahren integriert werden.
  -> Stabilisiert Lyapunov-Exponenten.
  -> Transformiert Gronwall-Abschaetzung von exponentiellem zu algebraischem Wachstum.
  -> Verbietet topologisch jede Singularitaet in endlicher Zeit.

- [ ] **EXT-NS-3: Zirkularitaet durchbrechen (Regularitaet <-> Attraktor)**
  Rote Flagge: Standardliteratur setzt Glaette VORAUS um Existenz/endliche Dimension
  des starken Attraktors zu beweisen. FST-NS muss diese Zirkularitaet explizit
  adressieren und aufloesen.
  -> Strikte Kontraktion von Q_N unter KL-Metrik beweisen.
  -> Lit: Zelik (2022) -- Attractors for Dissipative PDEs (modernstes Kompendium)

### Neue Referenzen (Navier-Stokes)
- [ ] Zelik (2022): Trajektorien-Attraktoren, Squeezing, Mane-Projektionen
- [ ] Barbu/Cannone (2016): Log-Lipschitz-Regularitaet, Zero Lipschitz Deviation
- [ ] Eden, Foias, Nicolaenko: Inertial Manifolds (Klassiker)
- [ ] "Geometric Redefinition of Turbulence" (2025): Fraktale Dimension schliesst
  Singularitaeten topologisch aus. Direkte Validierung des FST-NS-Ansatzes.

### 5.2 Yang-Mills -- Neue Aufgaben

- [ ] **EXT-YM-1: Kontinuumslimes-Problem explizit adressieren (KRITISCH)**
  Rote Flagge: Strong-Coupling-Massenluecke resultiert schlicht aus Gitter-Diskretheit.
  Im Limes beta->infinity (a->0) divergiert die Korrelationslaenge.
  Holley-Stroock-Konstanten wachsen exponentiell mit Variationspotential-Amplitude.
  -> Variationelle freie Energie explodiert im Weak-Coupling-Limit.
  -> Poincare-Ungleichungen werden wertlos (Massenluecke -> 0).
  -> LOESUNG: Hierarchische Block-Spin-RG INNERHALB des LSI-Frameworks.
     Skalenaufgeloeste Zerlegung der freien Energie statt globales Holley-Stroock.

- [ ] **EXT-YM-2: Reflexionspositivitaet im variationellen Ansatz erhalten**
  Ohne Reflexionspositivitaet keine gueltige Wightman-Theorie (Osterwalder-Schrader).
  -> Pruefen ob variationeller Ansatz RP bricht.
  -> Lemma 3.3 (Reflection Positivity) existiert, aber nur fuer Strong-Coupling.
  -> Beweis muss auf Weak-Coupling-Regime erweitert werden.

- [ ] **EXT-YM-3: LSI-Konstanten unter Skalentransformation binden**
  Beweisen: Holley-Stroock-Konstanten bleiben unter Renormierungsgruppen-Fluss
  invariant oder polynomiell beschraenkt (nicht exponentiell divergent).
  -> Birkhoff-projektive Kontraktion + KL-Minimierung zwischen RG-Skalen
     transformiert exponentielles in polynomielles Skalierungsproblem.
  -> Lit: "Emergence as Regime Formation" (2025) -- Polynomielle untere Schranken
     fuer SU(2)-Gittermassenluecke via Dobrushin-Shlosman + LSI.

### Neue Referenzen (Yang-Mills)
- [ ] "Emergence as Regime Formation" (2025/Recent): Exakte polynomielle untere
  Schranken via Dobrushin-Shlosman + LSI. Blueprint fuer FST-YM-Uebergang.
- [ ] Zegarlinski/Stroock (1992): Grundstein Log-Sobolev auf dem Gitter.
- [ ] Osterwalder-Schrader via Zegarlinski: Reflexionspositivitaet im Kontinuum.

### 5.3 Turbulenz -- Neue Aufgaben

- [ ] **EXT-TU-1: Energiekaskade als Minimax via Inertial Manifolds formalisieren**
  Null-Lipschitz-Abweichung auf dem Attraktor -> deterministische Versklavung
  der hochfrequenten Moden. Direkte Verbindung zur bestehenden offenen Idee
  "Energiekaskade als Minimax-Spiel".
  -> Kolmogorov-Spektrum als Sattelpunkt der skalenspezifischen KL-Divergenz.
  -> Lit: Eden, Foias, Nicolaenko -- Inertial Manifolds

### 5.4 Uebergreifend -- Strukturelle Isomorphie (Review-Einsicht)

> "Die Convexity Trap (RH) und die Gronwall-Luecke (NS) haben exakt dieselbe
> strukturelle mathematische Wurzel: Den Sprung von globalen makroskopischen
> Integralen zu lokaler Positivitaet/Regularitaet."

- [ ] **EXT-CROSS-1: SOS-Zertifikat als universelles Werkzeug**
  Wenn ein algebraisch geschlossenes SOS-Zertifikat fuer indefinite quadratische
  Formen (exakt 2 negative Eigenwerte) auf Hilbertraeumen konstruiert werden kann:
  -> Loest Convexity Trap (RH): Li-Koeffizienten-Positivitaet
  -> Loest Gronwall-Luecke (NS): Hochfrequente Moden-Kontrolle
  -> Potenziell auch YM-Kontinuumslimes
  -> "Universelles mathematisches Skalpell" (Review-Formulierung)
  -> Lit: Parrilo/Blekherman (SOS), Yakubovich (S-Lemma), Lasserre-Hierarchien

### Strategische Empfehlung (Review)

> "Der absolute theoretische Flaschenhals ist die Funktionalanalysis
> unendlich-dimensionaler Operatoren: SOS-Relaxierung in Hilbertraeumen
> und die unendlichdimensionale S-Procedure."

> Fokus YM: "Log-Sobolev-Ungleichungen unter RG-Fluss. Beweis nur haltbar wenn
> Holley-Stroock-Konstanten unter Skalentransformation invariant/polynomiell bleiben."

> Fokus NS: "Inertial Manifolds + Zero Lipschitz Deviation. Strikte Kontraktion
> von Q_N unter KL-Metrik beweisen um Zirkularitaet zu durchbrechen."

## 6. Knoten-Solitonen (Hopfionen) -- Neue Inputs (2026-03-11)

Basierend auf externer Analyse (INPUT.txt): Knoten-Solitonen-Mathematik (Hopfionen,
knotted fields, Faddeev-Niemi/Ackerman-Smalyukh) als strukturelle Ergaenzung.

### 6.1 CRM Paper III -- DIREKT INTEGRIERT (2026-03-11)

- [x] **Soliton-Ensemble als 6. UV-Completion-Kandidat** in Sec. 3.2 "UV-Completion Candidates"
  - Mean-Field-Ableitung: p(1-p) Kinetik -> X=2p-1 -> dX/da=k(1-X²)
  - Kapazitaetsauslastungs-Interpretation: X = C(a)/C_max
  - 4 neue Referenzen: Faddeev-Niemi 1997, Ackerman-Smalyukh 2017, Battye-Sutcliffe 1998, Vilenkin-Shellard 1994
- [x] **Saturation Universality Conjecture** erweitert um Kapazitaets-Remark (Sec. 4.4)
  - Drei Bedingungen aequivalent zu: endliches C_max, katalytisches Wachstum, monotoner Kontrollparameter
- [x] LaTeX kompiliert: 22 Seiten (vorher 21), keine Fehler

### 6.2 Navier-Stokes -- Offene TODOs

- [ ] **INPUT-NS-1: Helizitaets-Stratifizierung des Attraktors**
  Zerlegung A = Union_h A_h nach Helizitaets-Schichten.
  Stratifizierte Projektion v_h(t) statt nackter arg-min stabilisiert Minimierer-Auswahl.
  Koennte Assumption G (G1+G2) auf kontrollierter Auswahl begruenden.
  -> Lit: Moffatt (1969, Helicity in turbulence), Arnold & Khesin (1998)
  -> Substantielle theoretische Arbeit noetig. Mittlere Prioritaet.

- [ ] **INPUT-NS-2: Energie-Helizitaets-Ungleichungen als zweite Budget-Schranke**
  Topologie erzwingt Mindest-Enstrophie/Gradientenenergie fuer verknotete Konfigurationen.
  E >= c |Hopf|^{3/4} (Faddeev-Skyrme-Typ). Ergaenzt H(u|A)-Budget um topologische Koerzivitaet.
  -> Nicht 1:1 uebertragbar auf NS, aber Philosophie wertvoll. Niedrige Prioritaet.

- [ ] **INPUT-NS-3: Knot-Complexity als Attraktor-Regularisator**
  Funktionale Komplexitaet K(u) (Helizitaetsdichte-Spektrum, Linking-Statistik).
  Wenn K auf dem Attraktor gleichmaessig beschraenkt -> kontrollierte Minimierer-Variation.
  -> Reine Forschungsidee. Niedrige Prioritaet.

### 6.3 Turbulenz -- Offene TODOs

- [ ] **INPUT-TU-1: Vortex-Tubes als knotted Quasipartikel der Kaskade**
  Knottedness erhoeht lokale Dehnung/Curvature -> beeinflusst Transfer T_j.
  Topologischer Strafterm: F_lambda[E, structure] = F[E] + lambda E[knot complexity].
  Pruefen ob K41 robust bleibt waehrend Hessian-Fluktuationen "realistischer" werden.
  -> Passt zur bestehenden offenen Idee "Energiekaskade als Minimax-Spiel".

- [ ] **INPUT-TU-2: Entropy-Kompatibilitaet als "downhill in knot-space"**
  Annahme (NL) mikroskopisch begruenden: Rekonnexionen/Entknotungen sind dissipativ
  und bauen komplexe Konfigurationen ab -> "downhill"-Prinzip geometrisch statt
  probabilistisch. -> Remark in Discussion.

### 6.4 Yang-Mills -- Offene TODOs

- [ ] **INPUT-YM-1: Topologische Sektoren als "no flat directions"**
  Jede nichttriviale topologische Deformation (Aenderung robuster Observable)
  erzwingt positive Free-Energy-Kruemmung -> konzeptionelle Stuetze der "no flat
  directions"-Story. -> Remark nach Theorem 4.1 oder in Discussion.

- [ ] **INPUT-YM-2: Knotted Flux Tubes als Kandidaten fuer erste massive Stufe**
  Knotted flux tubes = Glueball-artige Zustaende. Liefert Interpretation
  WARUM der Gap thermodynamisch stabil ist und WELCHE Observablen O optimal koppeln.
  -> Niedrige Prioritaet, eher interpretativ.

### 6.5 Dark Energy -- Offene TODOs

- [ ] **INPUT-DE-1: Topologische Ladung als zusaetzlicher IR-Constraint**
  Phi[rho] erweitern zu Phi[rho, Q] mit Hopf-Invariante Q als Sektor-Constraint.
  Minimierung bei festem Q oder Sektor-Ensemble mit Entropie-Strafe.
  -> Passt zu Pattern A/B im Unified Framework. Mittlere Prioritaet.

- [ ] **INPUT-DE-2: Hopfion-Energiebounds als Backreaction-Form**
  E >= c |Q|^{3/4} als nichtquadratische Alternative zu V_grav ~ rho².
  Backreaction-Kostenfunktion aus Topologie + Steifigkeit statt nur Self-Gravity.
  -> Spekulativ, aber konzeptionell interessant. Niedrige Prioritaet.

### 6.6 RH Papers -- Dokumentiert (FROZEN)

> Alle RH Papers ausser RH3 sind eingefroren. Folgende Punkte nur als Referenz:

- [ ] **INPUT-RH-1: Holonomy-Monoid als Linking-Struktur (RH2)**
  Cocycle "no inverses, no cancellation" ist strukturell nah an topologischer
  Verknotung. Abbildung des Monoids in *-Algebra mit positiver Spur koennte
  A8/A9 als Kernel-Positivitaet liefern. -> Nur wenn RH2 reaktiviert wird.

- [ ] **INPUT-RH-2: Chern-Simons/Hopf-Term fuer positiven Fredholm-Kernel (RH2)**
  A_n als Quadratform in positivem Kernel: A_n = <f_n, K f_n> mit K >= 0.
  -> Nur wenn konkreter Kernel gebaut wird, nicht als Metapher.

### Referenzen (Knoten-Solitonen)

| Referenz | Relevanz | Bereits zitiert |
|----------|----------|----------------|
| Faddeev & Niemi (1997), Nature 387 | Grundlage Hopfionen | CRM Paper III |
| Ackerman & Smalyukh (2017), PRX 7 | Experimentelle Realisierung | CRM Paper III |
| Battye & Sutcliffe (1998), PRL 81 | Numerische Stabilitaet | CRM Paper III |
| Vilenkin & Shellard (1994), Cambridge UP | Kosmologische Defekte | CRM Paper III |
| Moffatt (1969) | Helizitaet in Turbulenz | -- (NS, TU) |
| Arnold & Khesin (1998) | Topologische Hydrodynamik | -- (NS) |
| Skyrme (1961) | Skyrmionen als Teilchen | -- (YM) |

---

## 8. RH-Durchbruch-Integration (2026-03-14)

> Die RH-Analyse hat eine universelle Architektur aufgedeckt:
> **Leading neutral -> 1. Ordnung degeneriert -> 2. Ordnung entscheidet**
> Dies ist das "Second-Order Resolvent Dominance Pattern" und gilt fuer alle 5 Probleme.

### 8.1 Universelles Resolvent-Muster (alle Papers)

- [x] **RH-CROSS-1: Second-Order Resolvent Pattern als universelle Architektur dokumentieren** -- ERLEDIGT (2026-03-14)
  Remark "Universal Resolvent Pattern" mit 5x4-Tabelle in alle 4 Folgebeweis-Papers eingefuegt
  (rem:universal-resolvent in YM, NS, TU, DE). Framework-Pattern-A Definition auf
  "Second-Order Resolvent Dominance" aktualisiert. Alle Papers kompiliert.

### 8.2 Navier-Stokes -- RH-Integration

- [x] **RH-NS-1: Resolvent-Muster im Attraktor-Minimierer** -- ERLEDIGT (2026-03-14)
  Remark rem:resolvent-minimiser nach Theorem 3 eingefuegt: rank-2 instabiles
  Residuum, Hellmann-Feynman-Zerlegung, H-Metrik Kontraktion als 2. Ordnung.

- [x] **RH-NS-2: Lipschitz-Kontrolle via Resolvent-Architektur** -- ERLEDIGT (2026-03-14)
  Assumption-G-Remark (rem:assumption-G) erweitert um "Resolvent-architectural
  reformulation": Sektor-Labels, projective Familie P_N, Sektor-Kreuzungs-Argument,
  epsilon-Netz-Schranke fuer endliche fraktale Dimension.

### 8.3 Yang-Mills -- RH-Integration

- [x] **RH-YM-1: Topologische Sektoren als Shift-Parity-Analog** -- ERLEDIGT (2026-03-14)
  Remark rem:shift-parity-YM nach Theorem 4.1 eingefuegt: Instanton-Sektoren als
  Shift-Parity-Analog, Kirk zero-free-strip als YM-Shift-Parity.
  Limitation 1 um topologische Sektor-Argumentation und Pattern-A-Referenz erweitert.

- [x] **RH-YM-2: Transfer-Operator Second-Order Resolvent Pattern** -- ERLEDIGT (2026-03-14)
  Remark rem:transfer-resolvent eingefuegt: Transfer-Operator 3-Level-Hierarchie
  (neutraler Vakuum-Eigenwert, degenerierende Poincare-Konstante, zweite-Ordnung
  Dobrushin-Zegarlinski Resummation als Gap-Mechanismus).

### 8.4 Turbulenz -- RH-Integration

- [x] **RH-TU-1: Joint-Minimierung als Second-Order Selektion** -- ERLEDIGT (2026-03-14)
  Step 5b "Second-order resolvent interpretation" im Theorem-K41-Beweis eingefuegt:
  rank-1 Hessian-Entartung, K41 als einziger stabiler Punkt via 2. Ordnung,
  Analogie zu E_sin/E_cos-Splitting.

- [x] **RH-TU-2: Resolvent-Spektrum der Kaskade** -- ERLEDIGT (2026-03-14)
  Neue Subsection 3.3 "Resolvent Dominance in the Cascade" eingefuegt:
  rank-0 Muster im Inertialbereich, rank-2 Korrektur am Dissipations-Cutoff,
  intermittente Bursts als transiente 1.-Ordnung-Verletzungen.

### 8.5 Dark Energy -- RH-Integration

- [x] **RH-DE-1: Kosmologische Konstante als Resolvent-Vakuum-Problem** -- ERLEDIGT (2026-03-14)
  Neue Subsection "The Resolvent Vacuum Problem" (sec:resolvent-vacuum) vor
  Discussion eingefuegt: 3-Level-Hierarchie (Leading neutral M_Pl^4,
  1. Ordnung Loop-Entartung, 2. Ordnung topologische Selektion rho_CC)
  mit Universal-Resolvent-Remark.

- [x] **RH-DE-2: Finsler-Diskriminante via Resolvent schaerfen** -- ERLEDIGT (2026-03-14)
  Predictions-Sektion um neuen Punkt 5 "Resolvent-vacuum prediction for eta"
  erweitert: eta=5/4 als Resolvent-Konsequenz, 3-Wege-Diskriminante
  CRM vs. Finsler vs. LCDM, DESI DR1 2027 Sensitivitaet.

### 8.6 Framework -- RH-Integration

- [x] **RH-FW-1: Pattern A generalisieren** -- ERLEDIGT (2026-03-14)
  FRAMEWORK_ABSTRACT.md war bereits auf "Second-Order Resolvent Dominance"
  aktualisiert (Section 3). FST_Unified_Framework.tex Pattern-A-Definition
  von "Finite Negativity as Gauge-Signal" auf "Second-Order Resolvent Dominance"
  umgestellt.

- [x] **RH-FW-2: FST_Unified_Framework.tex Zeilen 92-94 aktualisieren** -- ERLEDIGT (2026-03-14)
  Pattern A Punkt b) von "Finite Negativity as Gauge-Signal" auf
  "Second-Order Resolvent Dominance" aktualisiert. Kompiliert fehlerfrei.

- [x] **RH-FW-3: Euler-Maclaurin als universales Werkzeug** -- ERLEDIGT (2026-03-14)
  Neuer Paragraph nach der Euler-Maclaurin Endpoint Analysis Subsection eingefuegt:
  EM als universales Diagnostik-Tool, diskret-zu-kontinuierlich Uebergang in allen
  5 Problemen, Bernoulli-Korrekturen als endlich-rank Gauge-Shift.

---

## Quellen-Uebersicht

| Idee | NS | YM | DE | TU | Quelle |
|------|:--:|:--:|:--:|:--:|--------|
| Phi-Funktional (>= 0, KL, Kumulanten) | x | x | x | x | FST-GP Sec. 1-2 |
| CRM-Instanziierung (I_CRM) | | | x | | FST-GP Sec. 3 |
| Universelle Bruecke (Cost -> Response) | x | x | x | x | FST-GP Sec. 5 |
| Weil-Positivitaet / Funkt. Positivitaet | | x | | | FST-GP Sec. 6 |
| Minimax/Nash-Selektion | x | x | | x | FST-GP Sec. 6 |
| Legendre-Fenchel (S >= F) | x | | x | x | FST-GP Sec. 6 |
| K1-K5 Matrix (Pattern A/B) | x | x | x | x | FST-Meta Sec. 1-2 |
| SSA (S1)-(S4) | x | x | | | FST-GP Sec. 5.5 |
| Entropieproduktions-Funktional | x | x | x | x | FST-I Sec. 5 |
| MFG (HJB+FP) | x | | | x | FST-I Sec. 2 |
| Symmetriebrechung als Nash-Phasenuebergang | | x | x | | FST-I Sec. 6 |
| Replikator-Gleichung | x | | | x | FST-II Sec. 5.2-5.3 |
| Gibbs-Minimierung als Nash-GG | | | | x | FST-II Sec. 4.2 |
| Lyapunov dH/dt = -(E-E*)/T | x | | x | x | FST-II Sec. 5.4 |
| MEPP State Selection | | | x | x | FST-III Sec. 2.2 |
| Dissipative Adaptation | | | | x | FST-III Sec. 2.3 |
| Nash-Frustration / Jacobian | x | x | | | FST-III Sec. 3.3-3.4 |
| Proteinfaltung -> YM-Vakuum | | x | | | FST-III Sec. 3.2 |
| Skalar-Tensor-Gravitation (Lagrangian) | | | x | | CRM Paper 1-4 |
| Poeschl-Teller / Saettigung | x | x | x | x | CRM Paper 1-4 |
| BBN-Schutz (Trace Coupling) | | | x | | CRM Paper 3 |
| Falsifizierbare Vorhersagen (Slip, Lensing) | | | x | | CRM Paper 4 |
| Skalenabhaengige Kopplung | | | | x | CRM Paper 1-2 |
| Knoten-Solitonen / Hopfionen | x | x | x | x | INPUT.txt (2026-03-11) |
| Kapazitaetsauslastung (C/C_max) | | | x | | INPUT.txt -> CRM Paper 3 |
| Kirk O(4)-Restauration / Cluster Expansion | | x | | | LitReview (2026-03-12) |
| TMT Praegeometrie / Glue-Lattice Rigidity | | x | | | LitReview (2026-03-12) |
| Finsler-Gravitation (Konkurrent) | | | x | | LitReview (2026-03-12) |
| MOND-Entropy-EUP Verbindung | | | x | | LitReview (2026-03-12) |
| Hunt Informational Persistence | | | x | | LitReview (2026-03-12) |
| Hyperviskose Inertialmannigfaltigkeiten | x | | | | LitReview (2026-03-12) |

---

## 7. Literaturrecherche-Review -- Neue Inputs (2026-03-12)

Basierend auf erschoepfender externer Literaturrecherche und kritischer Analyse
des FST-Programms im Kontext des Forschungsstandes 2025/2026.

### 7.1 Yang-Mills -- Neue Erkenntnisse & TODOs

**Direkt integriert (2026-03-12):**
- [x] **Kirk (2026) + TMT (Baquero 2026) in Limitations referenziert**
  Kirk's O(4)-Restauration via Dobrushin-Shlosman complete analyticity +
  Pro-Polymer Cluster Expansion füllt CL1-CL3 direkt fuer SU(2).
  TMT als praegeometrischer Alternativansatz erwaehnt.
  3 neue Bibitems: Kirk2026, Kirk2026b, GarciaBaquero2026.

**Neue offene TODOs:**
- [~] **LIT-YM-1: Kirk's Methodik auf FST-Beweis abbilden (P1)** -- TEILWEISE ERLEDIGT (2026-03-12)
  Kirk nutzt Dobrushin-Shlosman complete analyticity -- genau die Methodik,
  die FST-YM in Theorem 4.1 (Holley-Stroock + Dobrushin-Zegarlinski) verwendet.
  -> Pruefen ob Kirk's Pro-Polymer-Abschaetzungen die FST-Konstanten direkt
     uebernehmen koennen (kappa_link, eta, beta_0).
  -> Falls ja: CL1-CL3 koennen als bewiesen markiert werden (zumindest fuer SU(2)).
  -> Lit: Kirk (2026), Zenodo DOI 10.5281/zenodo.14894820
  **Erledigt:** Remark rem:kirk-mapping mit Vergleichstabelle FST vs. Kirk eingefuegt.
  Numerische Verifikation: kappa_Haar=3, d_l=6, beta_0~0.056, eta(beta=2.2)=39.6>>1.
  Kirk's Pro-Polymer-Technik ueberbrueckt exakt die Luecke jenseits beta_0.
  **Offen:** Detailvergleich mit Kirk's tatsaechlichen Pro-Polymer-Abschaetzungen
  (erfordert Lektuere des vollstaendigen Papers).

- [ ] **LIT-YM-2: Vergleichstabelle FST vs. Kirk vs. TMT erstellen (P2)**
  3-Wege-Vergleich: Ontologischer Status des Gitters, UV-Problem-Loesung,
  Mechanismus der Massenluecke, Limitierungen.
  -> In Discussion als Tabelle oder als Appendix.

- [ ] **LIT-YM-3: Kirk's zero-free strip fuer RG-Fluss nutzen (P1)**
  Kirk's Companion-Preprint zeigt: Keine Phasenuebergaenge 1. Ordnung entlang
  des RG-Flusses. Schliesst Limitation 1 (Gap-scaling under renormalisation)
  zumindest fuer SU(2) partiell.
  -> In Limitation 1 referenzieren wenn Kirk peer-reviewed wird.

### 7.2 Dark Energy -- Neue Erkenntnisse & TODOs

**Direkt integriert (2026-03-12):**
- [x] **Neue Subsection "Competing approaches and recent developments"**
  Finsler-Gravitation als starker Konkurrent, MOND-Entropy-EUP-Stuetzung,
  Hunt Informational Persistence. 3 neue Bibitems.

**Neue offene TODOs:**
- [~] **LIT-DE-1: Finsler-Gravitation als Diskriminante nutzen (P1)** -- TEILWEISE ERLEDIGT (2026-03-12)
  Finsler-Modelle sagen eta != 5/4 voraus, CRM sagt eta = 5/4.
  -> Diesen Diskriminanten als "smoking gun test" in Predictions schaerfen.
  -> Euclid DR1 (2027) koennte entscheiden.
  **Erledigt:** Vollstaendige Herleitung eta=5/4 aus modifizierten Poisson-Gleichungen
  fuer f(R)=R+gamma*R^2 Skalar-Tensor-Theorie (quasi-statische Sub-Horizont-Naeherung).
  Crossover-Skala k_cross = a*m_s identifiziert. Finsler-Diskriminante als
  "model-discriminating observable" in Predictions-Sektion eingefuegt.
  Neue Subsection "Competing approaches and recent developments" hinzugefuegt.
  **Offen:** Quantitative Euclid-Sensitivitaetsanalyse (sigma(eta) vs. eta_CRM vs. eta_Finsler).

- [ ] **LIT-DE-2: MOND-Entropy-EUP Verbindung vertiefen (P2)**
  Jalalzadeh et al. leiten modifizierte Friedmann-Gleichungen aus
  verallgemeinerten Entropie-Flaechen-Relationen ab.
  -> Pruefen ob CRM-Friedmann als Spezialfall ableitbar ist.
  -> Stuetzt thermodynamische Ableitung von Lambda_eff.

- [ ] **LIT-DE-3: Dynamic Vacuum Field Theory als verwandten Ansatz einordnen (P2)**
  DVFT modelliert Vakuum als dynamisches Skalarfeld und leitet MOND ab.
  -> Abgrenzung zu CRM in Discussion. Gemeinsamkeiten: Skalarfeld + MOND.
  -> Unterschied: DVFT hat kein Entropie-/Thermodynamik-Argument.

### 7.3 Navier-Stokes -- Neue Erkenntnisse & TODOs

**Direkt integriert (2026-03-12):**
- [x] **Hyperviskose Modelle + Spectral Gap in Assumption-G-Remark**
  Inertialmannigfaltigkeiten fuer beta >= 5/4 bewiesen (Mallet-Paret/Sell).
  Klassischer Fall beta=1 offen wegen zu langsamen Eigenwertwachstums.
  Perfect-Slip H2-Resultate erwaehnt. Neues Bibitem MalletParetSell1988.

**Neue offene TODOs:**
- [ ] **LIT-NS-1: MFG-Wasserstein-Formulierung fuer NS vertiefen (P2)**
  Review identifiziert: Neueste MFG-Theorie arbeitet mit PDEs im
  Wasserstein-Raum und Master-Gleichungen (U Michigan 2026).
  -> Bestehendes MFG-Subsection (Sec. 5) koennte davon profitieren.
  -> Strukturelle Parallelen: Fokker-Planck ~ Advection-Diffusion.

- [ ] **LIT-NS-2: Zelik (2022) Trajektorien-Attraktoren einarbeiten (P1)**
  Zelik's Kompendium ist das modernste zu Attraktoren dissipativer PDEs.
  -> Fuer Assumption-G-Diskussion: Squeezing, Mane-Projektionen.
  -> Bereits in EXT-NS-3 vorgemerkt, aber Lit-Verweis noch nicht im Paper.

### 7.4 Turbulenz -- Neue Erkenntnisse

Keine neuen Erkenntnisse ueber bestehende TODOs hinaus.
Review bestaetigt "Clean Win"-Status und Joint-Minimierung als genuinen Beitrag.
Onsager-Verbindung und Duchon-Robert-Defektmass bereits im Paper.

### 7.5 Uebergreifende strategische Einsicht (Review)

> "Geigers thermodynamische Intuition, wonach die Massenluecke eine durch
> Zwangsbedingungen erzwungene globale Positivitaet darstellt (K1-Problem),
> ist konzeptionell brillant. Ohne die schweren analytischen Geschuetze
> der Clusterentwicklung (Kirk) oder einen radikalen ontologischen Wechsel
> zur diskreten Praegeometrie (Baquero) bleibt der Beweis im Rahmen der
> Standard-QFT jedoch unvollstaendig."

> "Das FST-Programm etabliert sich als brillantes, ordnendes Meta-Framework.
> Um den Anspruch auf Letztgueltigkeit zu erheben, muessen die aufgezeigten
> analytischen Luecken durch Integration spezifischer Methoden aus
> Gittereichtheorien, Machine Learning Potenzialen und stochastischen
> Differentialspielen geschlossen werden."

---

## Review-Findings 2026-03-14

> Vollstaendiger Review aller 7 Papers. Fokus auf 3 neue Papers (Hodge, BSD, P-vs-NP),
> Schnellscan der 4 bereits mehrfach reviewten Papers (YM, NS, TU, DE).
> Triviale LaTeX-Fixes direkt umgesetzt, substanzielle Punkte gesammelt.

### Triviale Fixes (direkt umgesetzt)

- [x] **Hodge, Zeile 137:** Lefschetz-Operator-Definition praezisiert: `$\alpha \mapsto L \cup \alpha$` -> `$\Lambda_L(\alpha) = L \cup \alpha$`
- [x] **P-vs-NP, Kommentare:** Doppelte Section-Nummerierung behoben (Sec 8 war doppelt). Kommentare 8->9, 9->10, 10->11 korrigiert.
- [x] **Hodge + P-vs-NP:** Kompiliert, PDFs aktualisiert (7 bzw. 8 Seiten, keine Fehler).

### 9.1 Hodge Paper -- Substanzielle Findings

- [ ] **HOD-1: Proposition 2.3 Beweis ist logisch unvollstaendig (MITTEL)**
  Der Beweis von "APC implies HC" sagt: "the question reduces to whether (AP2) and (AP3)
  are also satisfied." Aber das ist genau das, was gezeigt werden muesste -- es wird
  nichts reduziert, sondern nur umformuliert. Die Richtung "APC => HC" ist trivial
  (APC sagt alle AP-Klassen sind algebraisch, also sind alle Hodge-Klassen die AP
  sind auch algebraisch), aber der Beweis argumentiert umgekehrt und erzeugt so einen
  zirkulaeren Eindruck. Umformulieren: "If APC holds and if every Hodge class is
  arithmetically positive, then HC follows" -- das macht die zwei separaten Annahmen klar.

- [ ] **HOD-2: Lefschetz-Operator Notation inkonsistent (NIEDRIG)**
  $\Lambda_L$ wird als Cup-Product mit $L$ definiert (Zeile 137), aber in der
  Standardliteratur (Voisin, Griffiths-Harris) ist $\Lambda_L$ der duale/adjungierte
  Operator (Kontraktion), waehrend das Cup-Product einfach als $L$ oder $L \wedge$
  notiert wird. Das koennte Reviewer verwirren. -> Entweder $L_*$ fuer den Cup-Operator
  verwenden oder eine Fussnote zur Konvention hinzufuegen.

- [ ] **HOD-3: 7 Bibliographie-Eintraege nicht zitiert (NIEDRIG)**
  Grothendieck1969, Milne2013, Lewis1999, GrossZagier1986, Deligne1974, Cattani2014,
  Peters2008 -- keine \cite-Referenz im Text. Entweder als "Further Reading" markieren
  oder \nocite{*} verwenden, oder besser: relevante Stellen im Text mit \cite versehen.

- [ ] **HOD-4: Proof Sketch fuer Theorem 5.1 (Abelian Varieties) hat Luecke in Step 5 (HOCH)**
  Step 5 sagt: "The general case reduces to showing that the arithmetic positivity
  conditions on an abelian variety force the class to lie in the Tate module's
  Galois-fixed part, which is algebraic by the Tate conjecture for abelian varieties
  over finite fields (proved by Faltings)." Aber:
  (a) Der Uebergang von Hodge-Klassen ueber C zu Tate-Klassen ueber endlichen Koerpern
      ist selbst eine nicht-triviale Reduktion (benoetigt gute Reduktion + Vergleichssaetze).
  (b) "Galois-fixed part" der Tate-Vermutung betrifft l-adische Kohomologie, nicht
      direkt die Hodge-Klassen. Der Zusammenhang ist der Andre/Abdulali-Pfad, aber
      das wird nicht praezisiert.
  (c) Die Behauptung "The equivalence follows" nach Steps 1-5 ist zu stark --
      Step 5 ist ein Programm, kein Beweis.
  -> Als "Proof outline" oder "Strategy" labeln, nicht als "Proof".

- [ ] **HOD-5: AP1 Vorzeichen-Konvention (MITTEL)**
  AP1 verlangt $Q_L(\alpha_0, \alpha_0) \geq 0$ fuer primitive $(p,p)$-Klassen.
  Aber die Hodge-Riemann-Relationen sagen $(-1)^p Q_L(\alpha_0, \bar\alpha_0) > 0$
  fuer $\alpha_0 \in P^{p,p}$. Da $\alpha_0$ reell (rational) ist, gilt $\bar\alpha_0 = \alpha_0$,
  also $(-1)^p Q_L(\alpha_0, \alpha_0) > 0$, was fuer ungerades $p$ gerade
  $Q_L(\alpha_0, \alpha_0) < 0$ ergibt. Die AP1-Bedingung muesste also
  $(-1)^p Q_L(\alpha_0, \alpha_0) \geq 0$ lauten, nicht $Q_L \geq 0$.
  -> **KRITISCH:** Vorzeichen-Fehler in der zentralen Definition. Muss in Def 3.1,
  Eq (6), und allen Folge-Propositionen korrigiert werden.

- [ ] **HOD-6: Functoriality Proofs zu skizzenhaft (MITTEL)**
  Die Paper-eigene Limitation (Section 9, Punkt 2) gibt zu, dass die Proofs in
  Section 5 Skizzen sind. Das ist akzeptabel fuer ein Preprint, aber fuer ein
  Publikation muessten Pullback (Prop 5.1) und Pushforward (Prop 5.2) sauberer
  behandelt werden. Insbesondere: die Behauptung in Prop 5.1 dass
  $Q_{f^*L_X}(f^*\alpha_0, f^*\alpha_0) = \deg(f) \cdot Q_{L_X}(\alpha_0, \alpha_0)$
  gilt nur fuer generisch endliche Morphismen, nicht allgemein.

### 9.2 BSD Paper -- Substanzielle Findings

- [ ] **BSD-1: Theorem 4.1 ist zirkulaer aufgebaut (HOCH)**
  Das "BSD Normal-Form Theorem" (Theorem 4.1) setzt Axiome I-IV voraus und folgert
  BSD. Aber Axiom I (Sha-Endlichkeit) ist fuer Rang >= 2 offen und TEIL der BSD-Vermutung.
  Axiom IV (Euler-System) ist in voller Allgemeinheit ebenfalls unbewiesen. Somit
  reduziert das Theorem BSD auf "offene Probleme die BSD implizieren", was als
  Reformulierung wertvoll aber als Theorem ueberzogen ist. -> Als "Conditional
  Theorem" oder "Framework Theorem" labeln.

- [ ] **BSD-2: Step 2 der Proof-Strategie ist nicht bewiesen (HOCH)**
  "Composition forces the order" behauptet: Euler-System + L-Funktion => ord >= rank.
  Die Umkehrung "ord > rank => positivity violation" ist ebenfalls nur skizziert.
  Das Argument "the leading term vanishes faster than the regulator, violating
  positivity" ist eine Intuition, kein Beweis. -> Praezise Formulierung noetig.

- [ ] **BSD-3: Exclusion E2 hat logische Luecke (MITTEL)**
  E2 sagt: rank > ord impliziert R_E > 0 aber L^{(s)}(E,1) * (positive factors) = 0
  fuer s < r, "contradicting R_E > 0". Aber der Widerspruch folgt nur, wenn man
  die BSD-Formel bereits als wahr voraussetzt -- genau das, was bewiesen werden soll.
  -> Umformulieren als "under the assumption that the leading-term formula holds".

- [ ] **BSD-4: CRM-Referenz \cite{Geiger2026e} ist zirkulaer (NIEDRIG)**
  Das BSD-Paper zitiert das CRM Paper V als "Companion paper" fuer die Axiom-Analogien.
  Das ist in Ordnung fuer ein internes Forschungsprogramm, aber fuer ein externes
  Publikum wirkt die Selbstreferenz auf ein unpubliziertes "companion paper" schwach.
  -> Bei Publikation entweder Zenodo-DOI angeben oder den Verweis kuerzen.

- [ ] **BSD-5: Section 7 (Iwasawa) vermischt bewiesene und offene Ergebnisse (MITTEL)**
  Die Iwasawa-Hauptvermutung wird als "proven in many cases" zitiert, ohne zu
  praezisieren, in welchen Faellen sie fuer die Anwendung auf BSD ausreicht.
  Skinner-Urban gilt fuer gewoehnliche Primzahlen, Kato fuer alle p aber nur
  in eine Richtung (Teilbarkeit). -> Praeziser formulieren.

### 9.3 P-vs-NP Paper -- Substanzielle Findings

- [ ] **PNP-1: Korollar 4.3 (Witness Entropy Gap) hat eine vom Paper selbst erkannte Luecke (HOCH)**
  Das Paper erklaert korrekt in Remark 4.4, dass die Compression-Schranke fuer den
  *spezifischen* Witness gilt, den A produziert, waehrend die High-Entropy-Schranke
  fuer *alle* Witnesses gilt. Fuer Instanzen mit vielen Witnesses koennte A einen
  low-complexity Witness waehlen. Die Einschraenkung auf unique-SAT-Instanzen
  schliesst die Luecke theoretisch, aber die Konstruktion von dichten Familien
  solcher Instanzen (Strategy B1) ist nicht bewiesen.
  -> Das ist der Kernpunkt des Papers und richtig dargestellt, aber die Praesentation
  suggeriert mehr, als tatsaechlich bewiesen ist.

- [ ] **PNP-2: Proposition 3.1 (Barrier Immunity) ist zu optimistisch (HOCH)**
  (S1) behauptet: "K(x) does not change (up to O(1)) when oracles are added".
  Aber resource-bounded Kolmogorov-Komplexitaet K^t(x) haengt sehr wohl von Orakeln ab.
  Da das Paper spaeter K^t (Section 7) verwendet, nicht K, ist die Barriere-Immunität
  fuer den ressource-beschraenkten Fall NICHT etabliert. Insbesondere:
  - K^t relativiert (Baker-Gill-Solovay gilt analog fuer resource-bounded Kolmogorov).
  - Die Entropy Hardness Hypothesis (EH) in Theorem 7.2 ist eine Aussage ueber K^t,
    nicht ueber K, und unterliegt daher den Barrieren.
  -> Klar trennen: K (uncomputable) bypassed Barrieren, K^t (computable) nicht.

- [ ] **PNP-3: Entropic Separation Conjecture ist aequivalent zu P!=NP (MITTEL)**
  Das Paper gibt selbst zu (Sec 9, Punkt 3): "The No-Go theorem is conditional on
  the ESC, which is itself a reformulation of P != NP, not a proof." Das ist korrekt
  und ehrlich. Aber die Praesentation in Abstract und Introduction suggeriert mehr
  Fortschritt als tatsaechlich erreicht wird. -> In Abstract deutlicher machen, dass
  es sich um eine Reformulierung, nicht um einen Beweis handelt.

- [ ] **PNP-4: Theorem 4.2 (Existence of high-entropy witnesses) Proof Sketch unvollstaendig (MITTEL)**
  Der Proof Sketch konstruiert SAT-Instanzen mit unique satisfying assignments via
  "encoding a Kolmogorov-random string r". Aber:
  (a) Die Reduktion von r zu einer SAT-Instanz x_n erhaelt K(r|x_n) >= m - O(log n)
      nur, wenn die Reduktion injektiv ist und die SAT-Formel x_n "nicht mehr
      Information enthaelt als r". Das muss praezisiert werden.
  (b) Die Behauptung "every valid witness w_n satisfies K(w_n|x_n) >= |w_n| - O(log n)"
      gilt nur fuer den unique witness. Fuer allgemeine NP-Instanzen (mit vielen
      Witnesses) ist die Aussage falsch.
  -> "Proof sketch" praeziser machen oder Verweis auf Literatur (Li-Vitanyi).

- [ ] **PNP-5: AP4 (Monotone Stability) ist fragwuerdig als Axiom (NIEDRIG)**
  Die Behauptung, dass SAT "extreme sensitivity" zeigt (Zeile 346), ist zu pauschal.
  Fuer viele SAT-Instanzen ist die Satisfiability stabil unter kleinen Aenderungen
  (weit weg vom Phase Transition). Die Phase-Transition-Argumentation gilt nur fuer
  random k-SAT nahe der Schwelle, nicht allgemein. -> Differenzierter formulieren.

### 9.4 Alte Papers -- Schnellscan-Findings

- [ ] **YM: Kirk-Referenzen (Kirk2026, Kirk2026b) verifizieren**
  Die Referenzen sind als "Zenodo preprint, February/March 2026" angegeben.
  Pruefen ob tatsaechlich auf Zenodo verfuegbar (DOI 10.5281/zenodo.14894820).
  Falls nicht: als "private communication" oder "in preparation" labeln.

- [ ] **NS: \cite{Temam2001} doppelt im BibTeX** -- Nicht doppelt, aber nur 12 Referenzen
  insgesamt. Fuer ein 13-seitiges Paper duenn. -> Bei Publikation erweitern.

- [ ] **TU: \cite{FST-Unified2026} und \cite{FST-RH3} sind interne Referenzen**
  Beide Papers sind unpubliziert. Fuer externe Publikation problematisch.
  -> Zenodo-DOIs eintragen oder als "in preparation" kennzeichnen.

- [ ] **DE: Gleiches Problem mit internen Referenzen.**
  -> Wie TU behandeln.

- [ ] **ALLE 4 ALTEN: Pattern-A / Resolvent-Remarks prueferisch ueberarbeiten**
  Die am 2026-03-14 integrierten "Universal Resolvent Pattern" Remarks (rem:universal-resolvent)
  enthalten eine 5x4-Tabelle, die in allen 4 Papers identisch ist. Das ist methodisch
  sauber, aber bei einer externen Einreichung koennte die exakte Duplikation als
  self-plagiarism gewertet werden. -> Jedes Paper sollte nur die eigene Zeile +
  Kontext enthalten, nicht die volle Tabelle.

### 9.5 Uebergreifende Bewertung der 3 neuen Papers

> **Hodge:** Konzeptionell stark. Die Idee der "arithmetic positivity" als Vereinigung
> von Hodge-Riemann-Positivitaet, Absolutheit und Lefschetz-Kompatibilitaet ist originell.
> **ABER:** Vorzeichen-Fehler in AP1 (HOD-5) ist KRITISCH und muss vor Publikation
> korrigiert werden. Step 5 im Abelian-Varieties-Beweis ist ein Programm, kein Beweis.

> **BSD:** Gut strukturierte Reformulierung. Die Parallele Gross-Zagier = Positivity Identity
> ist erhellend. **ABER:** Theorem 4.1 ist zirkulaer (setzt offene Vermutungen als Axiome
> voraus). Die "Normal-Form"-Sprache suggeriert mehr als gezeigt wird.

> **P-vs-NP:** Ambitioniertestes der drei Papers. Barrier-Analyse ist solide.
> **ABER:** Der zentrale Punkt -- die Uniformity Bridge -- ist ungeloest, und die
> Barrier-Immunität-Behauptung (Prop 3.1) gilt streng genommen nur fuer unbeschraenktes K,
> nicht fuer K^t. Das Paper ist am ehrlichsten ueber seine Limitationen.

### Publikationsreihenfolge (erweitert um 3 neue Papers):

1. **Turbulenz** -- Clean Win (unveraendert)
2. **Navier-Stokes** -- Conditional Framework (unveraendert)
3. **Hodge** -- Stark nach AP1-Vorzeichen-Fix und Abelian-Reframing
4. **BSD** -- Wertvoll als Reformulierung, Axiom-Labelling noetig
5. **Yang-Mills** -- Reframing (unveraendert)
6. **P-vs-NP** -- Ueberzeugendste Barrier-Analyse, aber groesste offene Luecke
7. **Dark Energy** -- Framework Note (unveraendert)
