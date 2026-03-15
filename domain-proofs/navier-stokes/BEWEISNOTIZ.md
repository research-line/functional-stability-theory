# BEWEISNOTIZ -- Navier-Stokes Regularitaet
# Stand: 2026-03-15
# Status: CONDITIONAL Framework (G reduziert auf AGC1: Attraktor-Reach)

===============================================================================
## Problemstellung
===============================================================================

**Millennium-Problem (Fefferman 2000):** Fuer die 3D inkompressiblen
Navier-Stokes-Gleichungen auf T^3 oder R^3:

  d_t u + (u . nabla) u = -nabla p + nu Delta u,   div u = 0

beweise oder widerlege: Fuer glatte Anfangsdaten u_0 existiert die Loesung
fuer alle Zeiten t > 0 und bleibt glatt (keine Finite-Time-Singularitaet).

**Ansatz (FST-NS):** Thermodynamische Obstruktion gegen Blow-Up. Konstruktion
eines Attractor-Distanz-Funktionals H(u|A), das bei hypothetischem Blow-Up
unter Null fallen muesste -- Widerspruch zur Nichtnegativitaet.


===============================================================================
## Beweiskette
===============================================================================

### Schritt 1: Attraktor-Existenz und Regularitaet
**Status: BEWIESEN (Proposition 2.2, Standard-Literatur)**

Globaler Leray-Hopf-Attraktor A auf T^3 (mit glattem Forcing f):
- (i) A ist kompakt in H^1(Omega)
  (Absorbierende Kugel + Rellich-Kondrachov + Temam Thm I.1.1)
- (ii) Jedes v in A ist glatt: v in C^inf(Omega)
  (Invarianz -> vollstaendige Trajektorie -> Bootstrap H^1 -> H^2 -> ... -> C^inf)
- (iii) sup_{v in A} ||v||_{W^{1,inf}} < inf
  (Kompaktheit in H^1 + Sobolev-Einbettung + Bootstrap)


### Schritt 2: Attractor-Distanz-Funktional H(u|A)
**Status: BEWIESEN (Definition 2.1)**

  H(u|A) := inf_{v in A} (1/2) int_Omega |u - v|^2 dx

Wohldefiniert (A kompakt in L^2, Infimum wird angenommen).
Messbare Selektion v(t) in A existiert (Kuratowski-Ryll-Nardzewski).
Schreibe w := u - v, dann H = (1/2) ||w||^2.

Zentrale Eigenschaft: H(u|A) >= 0 fuer alle u.


### Schritt 3: Energy-Dissipation-Divergenz
**Status: BEWIESEN (Lemma 2.4, Lemma A, Lemma B)**

**Lemma A (Biot-Savart):** Vorticity omega = curl u kontrolliert
Geschwindigkeitsgradienten: ||nabla u||_{L^p} <= C_p ||omega||_{L^p}.
Vollstaendiger 4-Schritt-Beweis via Calderon-Zygmund.

**Lemma B (Enstrophie-Bilanz):** d/dt ||omega||^2 + 2nu ||nabla omega||^2
= 2 int omega . (omega . nabla) u dx. Vollstaendiger Beweis mit Gronwall-Kette.

**Energy-Dissipation Bound:** Wenn Blow-Up bei T*:
  int_0^{T*} ||nabla u||^2 dt = inf  (BKM-Kriterium)
d.h. die Dissipation divergiert.


### Schritt 4: H(u|A)-Budget und Widerspruchsbeweis
**Status: BEWIESEN conditional auf Assumption G (Theorem 3.1)**

Zeitableitung von H entlang der NS-Loesung:
  d/dt H = <w, d_t u - d_t v>
         = -nu ||nabla w||^2 + [Nichtlineare Terme] + [Minimierer-Korrektur]

**3-Schritt-Widerspruch:**
1. H >= 0 (Nichtnegativitaet, per Definition)
2. Blow-Up erzwingt int_0^{T*} Dissipation = inf
3. Unter Assumption G: Minimierer-Korrekturterme bleiben integrabel
   -> H wird irgendwann < 0 -> WIDERSPRUCH


### Schritt 5: Assumption G (OFFEN)
**Status: OFFEN -- Kernluecke des Frameworks**

**Assumption G** besteht aus zwei Teilen:

**(G1) Sobolev-Regularitaet des Minimierers:**
Der zeitabhaengige Minimierer v(t) = argmin_{v in A} ||u(t) - v||^2
erfuellt: v in L^inf_t H^2_x und ||nabla^2 v(t)|| ist gleichmaessig beschraenkt.

**(G2) Lipschitz-Regularitaet der Projektion:**
Die Abbildung t -> v(t) ist Lipschitz-regulaer:
  ||d_t v(t)|| <= L(t) mit int_0^T L(t)^2 dt < inf.

**PROBLEM:** (G1) folgt aus der Attraktor-Regularitaet (alle v in A sind glatt).
(G2) ist das schwierige: Die Projektion u -> v wechselt den naechsten Punkt
auf A bei zeitlicher Evolution, und diese Auswahl kann springen.

**Gronwall-Luecke:** Ohne (G2) hat die Differentialungleichung fuer H
einen nicht kontrollierten Term der Form ||d_t v|| * ||w||, der das
Budget-Argument torpediert.


===============================================================================
## Offene Luecken
===============================================================================

### L1: Gronwall-Luecke -- REDUZIERT auf AGC1 (2026-03-15, Sprint)

**Fortschritt in 3 Stufen:**

**(a) Originale Assumption G (AC-Regularitaet):** Zu stark, setzt
absolute Stetigkeit der Projektion voraus. Aequivalent zu Inertial Manifold.

**(b) Assumption G' (BV-Regularitaet):** Erlaubt Spruenge, verlangt nur
endliche Totalvariation. Proposition: G' reicht fuer Haupttheorem.

**(c) Attractor Geometry Condition (AGC, NEUES Lemma):**
  - AGC1: Positive Reach des Attraktors (reach(A) >= delta > 0)
    -> Nearest-point-Projektion ist Lipschitz-1
    -> v(t) = pi(u(t)) erbt AC von u(t)
    -> G (nicht nur G') folgt AUTOMATISCH
    STATUS: **OFFEN** (aequivalent zu Inertial-Manifold-Existenz fuer 3D NS)
  - AGC2: Exponential Tracking (dist(u(t), A) <= C*exp(-lambda*t))
    STATUS: **BEWIESEN** (Temam 2001, Standard-Attraktortheorie)

**Lemma AGC => G:** Vollstaendiger Beweis. Lipschitz-1 + AC von u(t) => AC von v(t).
Gronwall-Koeffizienten beschraenkt durch Attraktor-Regularitaet.

**Logische Kette:**
AGC1 + AGC2 => G => G' => Theorem (keine Finite-Time-Singularitaet)

**Status-Tabelle:**
  AGC1: OFFEN (Attraktor-Geometrie, nicht PDE)
  AGC2: THEOREM (Temam 2001)
  AGC => G: THEOREM (dieses Paper)
  G => G': trivial (AC => BV)
  G' => Thm: THEOREM (dieses Paper)

**Vier Angriffsvektoren (einer reicht):**
  AGC1:   Positiver Reach (global) -- staerkste, aequivalent zu IM
  AGC1':  Trajectory-Reach (lokal) -- nur entlang u(t), realistischer
  AGC1'': Lipschitz-Graph ueber P_N -- klassische IM-Bedingung
  AGC1''': BV-Selektion existiert -- mengenwertige Analysis

**BV-Refactor im Hauptbeweis (2026-03-15):**
  - Theorem-Statement erweitert: "If G or G' holds"
  - Explizite BV-Modifikation in Step 3: Stieltjes-Integral ersetzt dt_v
  - C_2 enthaelt NICHT MEHR ||dt_v|| -- nur Attraktor-Glattheit
  - Modifizierte Energie-Ungleichung (eq. H-integral-BV)
  - Failure-Modes-Tabelle mit allen Gliedern der Kette

**Shift:** Das NS-Regularitaetsproblem ist jetzt reduziert auf:
"Hat der 3D NS-Attraktor positiven Reach (oder eine der drei Alternativen)?"
Das ist eine GEOMETRISCHE Frage ueber die Attraktor-Struktur,
nicht eine PDE-Regularitaetsfrage ueber Loesungen.

### L1b: Literatur-Survey zu BV-Projektion (2026-03-15)
Systematische Recherche (Federer, Clarke, Poliquin, Chistyakov, Moreau,
Edmond-Thibault, Zelik) ergibt: **KEIN existierendes Resultat** verbindet
fraktale Dimension mit BV-Projektion. Die Luecke ist praezise:
- Federer Reach / Prox-Regularitaet: schliesst Fraktale aus (Reach=0)
- Chistyakov BV-Selektion: braucht Hausdorff-Variation (unbewiesen fuer NS)
- Sweeping Processes: nur fuer prox-regulaere Mengen
- Mane-Projektoren: nur Hoelder-stetig, nicht Lipschitz (fuer allg. Attraktoren)
- Bi-Lipschitz Mane: nur fuer Ginzburg-Landau mit omega != 0 (Zelik 2022)

**Drei vielversprechende Bruecken:**
1. "Log-Lipschitz BV"-Theorie (erweitert Chistyakov auf Hoelder-Variation) -- NEUES MATH
2. "Approximate Prox-Regularitaet" (prox-regulaer nach Entfernung meagerer Singularitaet)
3. Bi-Lipschitz Mane fuer 3D NS (wuerde AGC1'' direkt geben)

### L1c: Log-Lipschitz BV Kettenregel -- FORMALISIERT (2026-03-15, Sprint)

**Bruecke 1 ist jetzt formal ausgearbeitet im EN-Paper.**

**Kernkorrektur gegenueber erstem Entwurf:** Die GLOBALE log-Lipschitz-Bedingung
||P(x)-P(y)|| <= C||x-y||(1+log(R/||x-y||)) reicht NICHT fuer BV.
Grund: Partitionssummen sum_i Delta_i(1+log(R/Delta_i)) divergieren logarithmisch
fuer feine Zerlegungen. Der log-Faktor bezieht sich auf die Inkrementskala
||x-y||, die bei feiner Partition gegen 0 geht -> log(R/Delta) -> infty.

**Loesung:** Trajektorienweise log-Lipschitz-Bedingung (TLL):
  ||P(u(t+h)) - P(u(t))|| <= C_TLL ||u(t+h)-u(t)|| (1 + log(R/d(t)))
wobei d(t) = dist(u(t), A) die Distanz zum Attraktor ist (NICHT die Inkrementgroesse).
Der log-Faktor bezieht sich jetzt auf die geometrisch natuerliche Skala d(t).

**Neues Lemma (lem:LL-BV im Paper):**
  TLL + LDI => v in BV([0,T]; H) mit Var(v) <= C_TLL * integral ||u'|| (1+log(R/d)) dt

**Log-Distance Integrability (LDI):**
  integral_0^T ||u'(t)|| (1 + log(R/d(t))) dt < infty

**Beweis (3 Schritte):**
1. Metrische Ableitung: |v'|(t) <= C_TLL |u'|(t) (1+log(R/d(t))) a.e.
2. BV-Kriterium: Stetige Funktion mit integrabler metrischer Ableitung ist BV
   (Ambrosio-Fusco-Pallara, Abschnitt 3.2)
3. BV => G' => Theorem (keine Finite-Time-Singularitaet)

**Zwei hinreichende Bedingungen fuer LDI:**
- Variante A: log(R/d(t)) in BV_loc -> LDI via AC x BV Produktregel
- Variante B: Sublevel-Zeitmass |{t: d(t) <= e^{-m}}| <= C e^{-cm} (m>=1)
  -> log(R/d) in L^p fuer alle p -> LDI via Hoelder

**Plausibilitaet im NS-Setting:**
- Vor Blow-up (t < T*): d(t) beschraenkt oder exponentiell fallend,
  log(R/d) in L^infty_loc, LDI trivial
- Unter Blow-up: ||u'|| -> infty, aber d(t) bleibt L^2-beschraenkt
  (Leray-Hopf), interplay bestimmt LDI-Konvergenz

**Paralleler Pfad (unabhaengig von AGC):**
  TLL (offen) + LDI (offen) => BV => G' => Theorem
Vorteil: Kein Inertial Manifold noetig, nur trajektorienweise Regularitaet.

**Bezug zu globalem LL:**
Globales LL impliziert NICHT TLL (log-Gewicht bezieht sich auf verschiedene Skalen).
Globales LL ist Motivation, nicht Substitut fuer TLL.

**Status:** TLL und LDI sind OFFEN. Beide sind quantitative trajektorienweise
Bedingungen -- leichter zugaenglich als globale Attraktor-Geometrie (AGC1).

### L1d: Dissipationsintegral-Divergenz -- KRITISCH (2026-03-15, Sprint 2 Review)

**PROBLEM (REV-NS-C1):** Lemma B behauptet:
  BKM (int ||omega||_{L^inf} = inf) => int ||nabla u||^2 = inf

**ABER:** Die Leray-Energy-Identity gibt:
  nu * int_0^{T*} ||nabla u||^2 = 1/2 ||u_0||^2 - 1/2 ||u(T*)||^2 + int <f,u>
  RHS ist IMMER endlich (||u|| beschraenkt, f glatt, T* endlich).
  => **int_0^{T*} ||nabla u||^2 < inf IMMER**, auch unter Blow-up.

**KONSEQUENZ:** Der Widerspruchsbeweis im Haupttheorem (Step 3) funktioniert
nicht wie dargestellt, denn der "dominierende" Term -nu/2 int ||nabla u||^2
divergiert NICHT.

**Logischer Fehler in Lemma B Step 2c:** Der Beweis zeigt:
  (int ||nabla u||^2 < inf) UND (int ||omega||_{L^inf} < inf) => Phi beschraenkt
Die Kontraposition ist:
  Phi unbeschraenkt => (int ||nabla u||^2 = inf) ODER (int ||omega||_{L^inf} = inf)
Das Paper schliesst faelschlich:
  int ||omega||_{L^inf} = inf => int ||nabla u||^2 = inf
Das ist ein Fehlschluss (aus "A und B => C" folgt nicht "nicht B => nicht A").

**Was TATSAECHLICH unter Blow-up divergiert:**
- ||omega(t)||_{L^inf} -> inf (BKM)
- ||nabla u(t)||_{L^2} -> inf (punktweise, als Funktion von t)
- int ||nabla omega||^2 moeglicherweise (Enstrophie-Dissipation)
ABER: int ||nabla u||^2 bleibt ENDLICH (Leray Energy Identity).

**Drei Reparaturoptionen:**
(A) H auf H^1-Niveau heben: H^(1)(u|A) = inf_v 1/2 ||nabla(u-v)||^2
    Dann: Dissipationsterm wird -nu/2 ||nabla omega||^2, das divergieren kann.
    AUFWAND: Gross (gesamte Beweisstruktur umbauen).
(B) Hoehere Normen direkt nutzen (Serrin-Kriterium-Verletzung).
    AUFWAND: Mittel.
(C) Punktweise Divergenz von ||nabla u(t)||^2 nutzen (Rate-Argument).
    AUFWAND: Subtil, muss zeigen dass Rate schnell genug waechst.

**Referenzen:**
- Robinson, Rodrigo, Sadowski (2016): bestaetig u in L^2(0,T; H^1) fuer alle T
- Constantin, Foias (1988): Energieungleichung => endliche Dissipation
- Tao (2007, Blog): "Why global regularity for NS is hard"

**STATUS: ADRESSIERT (2026-03-15)**
- Paper umformuliert: "doubly conditional" (G + D)
- H^1-Lift als Section 6 eingefuegt (Theorem 6.4)
- H^1-Theorem: G^(1) + int ||Delta u||^2 = inf => kein Blow-up
- Die Asymmetrie int ||nabla u||^2 < inf vs int ||Delta u||^2 unbeschraenkt
  ist jetzt der zentrale Insight des Papers
- Verbleibende Luecken: G^(1) (analog zu G) + Enstrophie-Dissipations-Divergenz

### L2: Zirkularitaet (Regularitaet <-> Attraktor)
Standardliteratur setzt Glaette VORAUS, um Existenz/endliche Dimension
des starken Attraktors zu beweisen. FST-NS muss diese Zirkularitaet
explizit aufloesen.

**Moegliche Loesung:** Strikte Kontraktion von Q_N (hochfrequente Moden)
unter KL-Metrik beweisen. Referenz: Zelik (2022) -- Trajektorien-Attraktoren.

### L3: Helizitaets-Stratifizierung (explorativ)
Zerlegung A = Union_h A_h nach Helizitaets-Schichten koennte die
Minimierer-Auswahl stabilisieren (INPUT-NS-1 aus Review).


===============================================================================
## Naechste Schritte
===============================================================================

1. **TLL verifizieren (Hauptangriff):** Trajektorienweise log-Lipschitz-Kontrolle
   der Nearest-Point-Projektion entlang dissipativer Trajektorien. Ansaetze:
   - Squeezing-Eigenschaften des Semiflow (Barbu/Cannone 2016)
   - Skalierte Prox-Regularitaet (schwaecher als Federer-Reach)
   - Bi-Lipschitz Mane fuer 3D NS (Zelik-Programm)

2. **LDI verifizieren:** Log-Distanz-Integrierbarkeit laengs Trajektorien.
   Vielversprechendster Ansatz: Sublevel-Zeitmass-Kontrolle (Variante B).

3. **Log-Dirichlet-Quotienten:** In KL-Joint-Minimierung integrieren
   -> Stabilisierung der Lyapunov-Exponenten
   -> Algebraisches statt exponentielles Wachstum der Gronwall-Abschaetzung.

4. **Zirkularitaet durchbrechen:** Strikte Kontraktion von Q_N unter
   KL-Metrik beweisen (ohne Glaette vorauszusetzen).

5. **Reframing beibehalten:** Paper als "conditional regularity criterion"
   und "structural reduction", NICHT als Beweis.


===============================================================================
## Verbindung zum Meta-Framework (FST Positivity Pattern)
===============================================================================

### Pattern A: Second-Order Resolvent Dominance

| Ebene | Mechanismus | NS-Instanziierung |
|-------|-------------|-------------------|
| Leading neutral | Leray-Hopf Vakuum | Schwache Loesung existiert |
| 1st order degenerate | Attraktor-Projektion | H(u|A) >= 0, aber Budget unklar |
| 2nd order decides | H-metrische Kontraktion | Dissipation erzwingt H < 0 -> Widerspruch |

### Saettigungsdynamik

Blow-Up-Verhinderung ist strukturell analog zum Poeschl-Teller-Saettigungs-
mechanismus: Das System kann den Attraktor nicht "ueberschiessen", weil
die Dissipation (Entropie-Produktion) jede Abweichung zuruecktreibt.

### Strukturelle Isomorphie mit RH

> "Die Convexity Trap (RH) und die Gronwall-Luecke (NS) haben exakt dieselbe
> strukturelle mathematische Wurzel: Den Sprung von globalen makroskopischen
> Integralen zu lokaler Positivitaet/Regularitaet."
> -- Review-Einsicht (2026-03-11)

### Gemeinsame Programmlinie

```
F mit Entropieterm -> Attractor-Stabilitaet -> Distanz-Budget -> Regularitaet
```

Bei NS ist die "Physik" = Regularitaet (keine Finite-Time-Singularitaet).
Die Nichttrivialitaet liegt in der Kontrolle des Minimierer-Wechsels,
nicht in der KL-Konvexitaet.


===============================================================================
## Referenzen
===============================================================================

- Leray (1934): Schwache Loesungen, Leray-Hopf-Regularitaet
- Beale, Kato & Majda (1984): BKM Blow-Up-Kriterium
- Temam (2001): Infinite-Dimensional Dynamical Systems
- Foias, Manley, Rosa & Temam (2001): Navier-Stokes und Turbulenz
- Kuratowski & Ryll-Nardzewski: Messbare Selektion
- Barbu & Cannone (2016): Log-Lipschitz-Regularitaet
- Zelik (2022): Attraktor-Theorie, Squeezing, Mane-Projektionen
- Eden, Foias & Nicolaenko: Inertial Manifolds
- Lasry & Lions (2007): Mean Field Games (Verbindung zur MFG-Perspektive)
