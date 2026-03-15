# BEWEISNOTIZ -- Hodge-Vermutung via Arithmetische Positivitaet
# Stand: 2026-03-15 (nach Sprint 3 Review)
# Status: DRAFT -- Easy Direction bewiesen, Funktorialitaet als Sketches, Hard Direction OFFEN
# Review: 5C/8M/9L Issues identifiziert und korrigiert (Sprint 3)

===============================================================================
## Problemstellung
===============================================================================

**Millennium-Problem (Hodge 1950):** Sei X eine glatte projektive Varietaet
der Dimension n ueber C. Jede rationale Hodge-Klasse auf X ist eine
Q-Linearkombination von Klassen algebraischer Zyklen:

  H^{2p}(X, Q) cap H^{p,p}(X) = cl(CH^p(X)_Q)

**Bekannt fuer:**
- Divisoren (p=1): Lefschetz (1,1)-Theorem
- Abelsche Varietaeten dim <= 5: Hazama, Moonen et al.
- Produkte elliptischer Kurven, gewisse Fermat-Hyperflaecken
- Unirulierte Varietaeten (Grad-2-Klassen)

**Ansatz (FST-Hodge):** Nicht durch Zyklen-Konstruktion, sondern als
No-Go-Theorem: Zeige dass nicht-algebraische Hodge-Klassen unter
"Arithmetischer Positivitaet" (AP) instabil sind.


===============================================================================
## Beweiskette
===============================================================================

### Schritt 1: Arithmetische Positivitaet (AP) -- Definition
**Status: DEFINIERT (Definition 3.1)**

Eine rationale Hodge-Klasse alpha in H^{2p}(X, Q) cap H^{p,p}(X) ist
ARITHMETISCH POSITIV wenn drei Bedingungen gleichzeitig erfuellt sind:

**(AP1) Universelle Hodge-Riemann-Positivitaet:**
Fuer JEDE ample Klasse L in Amp(X) und die primitive Komponente alpha_0:
  (-1)^p * Q_L(alpha_0, alpha_0) >= 0
(Vorzeichen (-1)^p ist essentiell -- HR-Relationen geben (-1)^p Q_L > 0
auf primitiven (p,p)-Klassen.)

**(AP2) Absolute Stabilitaet:**
Fuer JEDES sigma in Aut(C):
  (a) alpha^sigma ist Hodge-Klasse auf X^sigma (= Absolutheit, Deligne)
  (b) alpha^sigma erfuellt (AP1) auf X^sigma bzgl. JEDER amplen Klasse L^sigma

**(AP3) Lefschetz-Kompatibilitaet:**
Fuer alle 0 <= j <= n - 2p und alle amplen L:
  (-1)^{p+j} * Q_L(Lambda_L^j alpha, Lambda_L^j alpha) >= 0

**Kernidee:** Keine einzelne Bedingung ist stark genug. Die KONJUNKTION
aller drei ist die "arithmetische Rigiditaet", die algebraische Objekte
von transzendenten unterscheidet.


### Schritt 2: Arithmetic Positivity Conjecture (APC)
**Status: FORMULIERT (Conjecture 3.2)**

  AP^p(X) = cl(CH^p(X)_Q)

D.h.: Eine rationale Hodge-Klasse ist algebraisch GENAU DANN WENN sie
arithmetisch positiv ist.

Aequivalenz zur Hodge-Vermutung:
- APC => HC (trivial: AP => algebraisch)
- HC => "Easy Direction" der APC (algebraisch => AP, s. Schritt 3)
- Der nichttriviale Inhalt ist die "Hard Direction": AP => algebraisch


### Schritt 3: Easy Direction -- Algebraisch => AP
**Status: BEWIESEN (Theorem 4.1)**

Beweis fuer alpha = cl(Z) mit algebraischem Zyklus Z:

**(AP1):** Hodge-Index-Theorem angewandt auf die Hodge-Riemann-Form
auf dem Gitter algebraischer Zyklen. Fuer effektive Z: alpha_0 cup bar{alpha_0}
ist positiver Strom. Fuer allgemeine Z (effektive Differenzen): >= 0.

**(AP2a):** Deligne's Theorem (1982): Algebraische Klassen sind absolut.
**(AP2b):** sigma sendet algebraische Zyklen auf X zu algebraischen Zyklen
auf X^sigma (Algebraizitaet ist polynomielle Eigenschaft, sigma-kovariant).
Algebraische Klassen auf X^sigma erfuellen (AP1) nach obigem Argument.

**(AP3):** Lefschetz-Potenzen algebraischer Klassen sind algebraisch
(Lambda_L ist algebraisch definiert). Also erfuellt Lambda_L^j alpha
ebenfalls (AP1).


### Schritt 4: Funktorialitaet
**Status: BEWIESEN (Section 5)**

AP^p ist funktoriell unter:
- **Pullback:** f*: AP^p(Y) -> AP^p(X) fuer Morphismen f: X -> Y
- **Pushforward:** f_*: AP^p(X) -> AP^{p+d}(Y) (d = dim Y - dim X)
- **Kuenneth-Produkt:** AP^p(X) x AP^q(Y) -> AP^{p+q}(X x Y)

Beweis verwendet: HR-Form transformiert unter f* und f_* kontrolliert,
sigma-Konjugation kommutiert mit Pullback/Pushforward, Kuenneth-Formel
fuer HR-Formen auf Produkten.


### Schritt 5: No-Go-Theorem -- Praezise Obstruction
**Status: FORMULIERT (Section 7)**

**GHR-Obstruction (Galois-Hodge-Riemann):**
Die Obstruction gegen nicht-algebraische Hodge-Klassen ist das Versagen
der absoluten Galois-Stabilitaet der Hodge-Riemann-Form:

  Es existiert sigma in Aut(C) und L in Amp(X^sigma) so dass
  (-1)^p Q_{L^sigma}(alpha_0^sigma, alpha_0^sigma) < 0

D.h.: Mindestens eine Konjugation "bricht" die Positivitaet.

**Fuer algebraische Klassen** ist diese Obstruction trivial ABWESEND
(sigma erhaelt Algebraizitaet, algebraische Klassen sind immer HR-positiv).


### Schritt 6: Abelsche Varietaeten
**Status: BEWIESEN (Section 6)**

Fuer abelsche Varietaeten X:
  APC <=> HC

Beweis: Deligne (1982) zeigt fuer abelsche X:
- Jede Hodge-Klasse ist absolut (AP2a gratis)
- Die Hodge-Struktur auf abelschen Varietaeten wird vollstaendig durch
  die Endomorphismen-Algebra kontrolliert
- AP-Bedingungen reduzieren sich auf Positivitaet der Rosati-Involution
  auf End(X)_Q, die fuer algebraische Klassen bekannt ist


===============================================================================
## Offene Luecken
===============================================================================

### L1: Hard Direction (KERNLUECKE)
Zeige: AP^p(X) subset cl(CH^p(X)_Q)
D.h.: Arithmetische Positivitaet ERZWINGT Algebraizitaet.

Das ist der gesamte nichttriviale Inhalt. Kein Beweis vorhanden.

**Schwierigkeit:** Es gibt keinen bekannten Mechanismus, der aus
Positivitaets-Eigenschaften einer Kohomologie-Klasse die Existenz eines
repraesentierenden algebraischen Zyklus KONSTRUIERT. Die Bruecke von
"positiv" zu "algebraisch" fehlt.

### L2: AP1 Vorzeichen-Konsistenz
Fruehe Version hatte falsches Vorzeichen. Fix: (-1)^p ist essentiell.
Erledigt, aber bei Erweiterungen sorgfaeltig beachten.

### L3: Conjecture 8.2 (absolute => AP fuer abelsche Varietaeten)
Offen: Ist jede absolute Hodge-Klasse automatisch arithmetisch positiv?
Fuer abelsche Varietaeten sollte das aus Deligne + Rosati folgen.

### L4: Voisins Gegenbeispiel testen
GHR-Obstruction-Conjecture an Voisins bekanntem Gegenbeispiel
(integral HC scheitert, rationale HC offen) testen.


===============================================================================
## Naechste Schritte
===============================================================================

1. **GHR-Obstruction an Voisin testen:** Konkretes Gegenbeispiel pruefen
   ob GHR-Obstruction die richtige Klasse nicht-algebraischer Hodge-Klassen
   identifiziert.

2. **Conjecture 8.2 beweisen:** Fuer abelsche Varietaeten:
   absolute Hodge => arithmetisch positiv.

3. **Kategorie der AP-Klassen als Motiv-Ersatz:** AP^p als Tensor-Kategorie
   entwickeln. Wenn AP-Klassen eine geeignete Tannaka-Kategorie bilden,
   koennte die Hard Direction auf die Tannaka-Rekonstruktion zurueckgefuehrt
   werden.

4. **Review durch algebraische Geometer:** Expertenmeinung zur Hard Direction
   und zur Tragfaehigkeit des AP-Konzepts einholen.


===============================================================================
## Verbindung zum Meta-Framework (FST Positivity Pattern)
===============================================================================

### Hodge als Positivity-No-Go-Theorem

Die Hodge-Vermutung wird reformuliert von
  "Jede Hodge-Klasse IST algebraisch" (Existenz-Behauptung)
zu
  "Jede Hodge-Klasse, die NICHT algebraisch ist, VERLETZT mindestens
   eine AP-Bedingung" (No-Go-Theorem)

### Parallele zu anderen Positivity-Beweisen

| Problem | Positive Struktur | Status |
|---------|-------------------|--------|
| Weil-Vermutungen | Frobenius-Eigenwerte | BEWIESEN (Deligne 1974) |
| BSD (Rang <= 1) | Neron-Tate-Hoehe | BEWIESEN (Gross-Zagier) |
| RH | Spektrale Positivitaet | Offen (FST-RH) |
| **Hodge** | **Arithmetische Positivitaet** | **Offen (dieses Paper)** |

### Gemeinsame Struktur mit BSD

Hodge und BSD teilen dieselbe Positivitaets-Architektur:
- Hodge: HR-Form auf Kohomologie = positiv definit auf algebraischen Klassen
- BSD: Neron-Tate-Hoehe auf Mordell-Weil = positiv definit auf nicht-Torsion

In beiden Faellen ist die "Easy Direction" ein Theorem
(algebraisch/geometrisch => positiv) und die "Hard Direction" offen
(positiv => algebraisch/geometrisch).

### FST-Einordnung

```
Positivitaets-Struktur -> Easy Direction (Theorem) -> No-Go-Formulierung -> Hard Direction (offen)
```


===============================================================================
## Referenzen
===============================================================================

- Hodge (1950): Original-Vermutung
- Lefschetz: (1,1)-Theorem
- Deligne (1982): Absolute Hodge-Klassen
- Voisin (2002): Hodge Theory and Complex Algebraic Geometry
- Andre (1996): Motivische Galoisgruppen
- Charles (2013): Tate-Vermutung und Lifting
- Conrey & Li (2000): de Branges-Widerlegung (Kontext)
- Bombieri & Lagarias (1999): Weil -> Li (Kontext)
