# Physics-Master — RFEP v1.8

**Rolle:** Ebene 2 der FST-Hierarchie — physikalische Instantiation von Functional Stability Theory via Renormalized Free Energy Principle.

**Status:** v1.8 ist seit 2026-05-30 auf Zenodo live: Record `20467623`, DOI `10.5281/zenodo.20467623`, Concept-DOI `10.5281/zenodo.19036190`. Der aktive lokale Dateisatz ist weiter die quellengeprüfte Manuskriptfassung v0.9, wurde aber nach dem Live-v1.8-Upload am 2026-06-04/12 erneut gebaut. EN/DE/Kombi weichen daher am 2026-06-18 per MD5 vom Live-Record ab und sind nur ein v1.9-Vorbereitungskandidat, nicht der Live-Dateisatz.

RFEP ist nicht nur ein Begleitpaper, sondern die Übersetzungsschicht des Physikstamms des Zeta-Zoo: Operatorgewinne aus Zookeeper, Prime-Hub und Zoo-Mapping werden hier in eine Sprache übertragen, die für Yang-Mills, Navier-Stokes, NS-LDI, Turbulenz und Cosmology nutzbar ist.

## Haupt-Repository

- GitHub: https://github.com/research-line/functional-stability-theory
- Altes Repo: https://github.com/research-line/rfep-framework (ausphasiert)

## Relevante tree-interne Dokumente

- `CORE/UNIFIED_TERMINOLOGY.md` — Namenskonvention FST/RFEP
- `CORE/CONSOLIDATION_PLAN.md` — Vier-Ebenen-Hierarchie
- `CORE/RFEP/RFEP_V1_8_UPDATE_PLAN.md` — Update-Protokoll v1.7 -> v1.8
- `CORE/RFEP/BEWEISNOTIZ.md` — Post-Zookeeper-Beweiskette für RFEP als Physik-Adapter
- `CORE/RFEP/notes/ZOOKEEPER_RFEP_TRANSFER.md` — technisches Transfer-Schema aus Zookeeper/Prime-Hub nach RFEP
- `CORE/zookeeper/proof_notes/` — Three-Lemma-Endgame, Mass Concentration, Coercive Complement
- `CORE/prime_hub/` — Bouquet-/Frontier-Prime-Adapter für physikalische Randdaten

## Lokaler Spiegel

Ein Spiegel der alten v1.7-Quellen liegt in `_archive/old_fst_framework_local_mirror/` zur Referenz. Aktive Paperdateien liegen lokal in `paper/`; der öffentliche GitHub-Pfad ist `https://github.com/research-line/functional-stability-theory/tree/main/masters/spectrum-duality`.

## Brücke zum Math-Master

- Abschnitt UCU <-> Pattern A wird im Math-Master (`CORE/zetazoo/paper/NE_B_BOUNDARY_v1_en.tex`) in Sektion 7.5 ausgearbeitet. Der Physics-Master verweist auf diese Brücke in seinem Intro.

## Post-Zookeeper-Gewinne für RFEP

**Guardrail 2026-05-17:** Der Zookeeper-Transfer wird hier als bedingte
Mikrocluster-Normalform gelesen. Die uniforme Kontrolle von `s_lambda`,
`p_lambda` und `g_*` bleibt ein explizites Gate; RFEP importiert daraus kein
automatisches Physik- oder RH-Theorem.

**Strict-Reviewer 2026-05-20:** FST-7 wird nur noch als bedingtes
Kompatibilitätsschema geführt, massendefinierte Clusterlücken sind
posteriori-Diagnostik, und Domain-Transfers benötigen unabhängige
Zielsubspace- oder Außengap-Zertifikate.

**Quellencheck 2026-05-26:** Die aktive lokale Manuskriptfassung wurde von
v0.8 auf v0.9 angehoben. Korrigiert wurden Connes 2026, CCM 2024,
RH-Landscape/Even-Dominance, Cornelissen--Marcolli 2010 sowie fehlende
DOI-Metadaten; der v0.8-Dateisatz liegt archiviert in
`_archive/quellencheck_2026-05-26_v0-8_pre-sourcecheck/`. Dieser v0.9-Dateisatz
ist seit dem Zenodo-v1.8-Upload live; die Live-Dateinamen enthalten das harmlose
Tool-Rename-Artefakt `FST_SPECTRUM_DUALITY-9_*`. Die späteren lokalen
post-v1.8-PDFs sind neuer als live und brauchen vor Veröffentlichung eine echte
v1.9-New-Version plus GitHub-Sync und finales lokales Upload-Gate.

Der Zookeeper-Gewinn wird in RFEP als wiederverwendbares Spektralpaket aufgenommen:

1. **Rank-one / finite-rank defect language:** Physikalische Obstruktionen werden als kleine Operator- oder Formdefekte formuliert.
2. **Mass-based cluster gap:** Nicht nur ein einzelner Eigenwert zählt, sondern die spektrale Masse in einem stabilen Cluster.
3. **Coercive complement:** Der relevante RFEP-Minimierer wird durch eine Lücke gegen das orthogonale Komplement geschützt.
4. **Galerkin faithfulness:** Endliche Trunkierungen werden nicht als bloße Numerik behandelt, sondern als Beweisverpflichtung mit explizitem Grenzschritt.
5. **Prime-Hub boundary transfer:** Prim-/Randdaten werden als Frontier- oder Bouquet-artige Randterme lesbar.

## Physik-Adapter

| Physikzweig | RFEP-Zookeeper-Nutzung | Status |
|---|---|---|
| Yang-Mills | Mass gap / Gibbs uniqueness als coercive-complement-Problem | Kandidat |
| Navier-Stokes / NS-LDI | Dissipative Selection und Log-Distance-Defekt als Spektralcluster-Kontrolle | Kandidat |
| Turbulenz K41 | Energie-Kaskade als stabiler RFEP-Minimierer mit Clusterlücke | Kandidat |
| Cosmology / CRM | UCU/strikte Konvexität bleibt Hauptweg; Zookeeper liefert Operator-Sprache für Perturbationen | unterstützend |
| Prime-Hub | Randdaten/Frontier-Primes als physikalisch übersetzbare Boundary-Operatoren | hoher Nutzen |
