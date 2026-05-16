# Quellencheck 2026-05-16 -- Dirichlet Atlas / Zoo-Mapping

## Gegenstand

- Projekt: `.LAB/.ZETA-ZOO/CORE/zoo-mapping/`
- Aktive neue Dateien:
  - `paper/DIRICHLET_CHARACTER_ATLAS_v2_en.tex`
  - `paper/DIRICHLET_CHARACTER_ATLAS_v2_de.tex`
  - `paper/DIRICHLET_CHARACTER_ATLAS_v2_en.pdf`
  - `paper/DIRICHLET_CHARACTER_ATLAS_v2_de.pdf`
  - `paper/DIRICHLET_CHARACTER_ATLAS_v2_kombi.pdf`
- Ausgangslage: Zenodo v1 ist live (`10.5281/zenodo.19960810`); der lokale v1-Stand enthielt bereits den Zitationscheck vom 2026-05-15 und wurde für diesen Quellencheck als v2-Kandidat fortgeführt.

## Externe Prüfung

Geprüft wurden die aktiven 16 Bibitems und die quellenbezogene LMFDB-Behauptung gegen Web- und Datenbankquellen:

- Zenodo Records/API für interne FST-/RH-/Zookeeper-/Zeta-Zoo-/Dirichlet-Atlas-DOIs:
  `19673226`, `19036190`, `19673126`, `19035640`, `19960810`, `19960809`, `19962588`, `19764771`.
- Crossref für Connes 1999 (`10.1007/s000290050042`).
- AMS Bookstore für Iwaniec--Kowalski 2004 und Paley--Wiener 1934.
- CERN/i-Scholar/Webkataloge für Weil 1952 und Selberg 1956.
- LMFDB für die Datenbankidentität und eine Dirichlet-L-Funktionsseite (`1-5-5.4-r0-0-0`), die berechnete erste Nullstellen auf der kritischen Geraden anzeigt.

## Befund

Die bibliographischen Kernangaben waren nach dem Zitationscheck weitgehend korrekt. Ein inhaltlicher Quellenfehler blieb:

- Die EN/DE-Fassungen formulierten, LMFDB berichte für alle getesteten Führer alle nicht-trivialen Nullstellen unterhalb `T ~ 10^4` als auf der kritischen Geraden verifiziert.
- Die geprüfte LMFDB-L-Funktionsseite belegt berechnete erste Nullstellen auf der kritischen Geraden, aber nicht den stärkeren Zertifikatsclaim für alle Nullstellen bis zu einer festen Höhe.

Zusätzlich wurden zwei Dokumentationsdrifts bestätigt:

- Das Projekt war live v1, aber lokal bereits bibliographisch weiter als die veröffentlichte Fassung.
- Einzelne Plan-/Credential-Labels sprachen noch von `RH Even Dominance v2.1`, obwohl die aktive DOI in der Dokumentation als Concept-DOI geführt wird und Zenodo latest v2.3 meldet.

## Korrekturen

- Neue v2-TeX-Dateien aus dem lokalen v1-Stand erstellt.
- EN/DE-LMFDB-Passage zurückgenommen: LMFDB wird jetzt als Quelle für berechnete niedrig liegende Nullstellen und endliche Konsistenzdaten beschrieben, nicht als Zertifikat für alle Nullstellen unterhalb einer bestimmten Höhe.
- LMFDB-Bibitem-Zugriffsdatum auf 2026-05-16 aktualisiert.
- v1-Dateisatz nach `paper/_archive/quellencheck_2026-05-16_v1/` verschoben.

## Verifikation

- `pdflatex` für EN und DE je zweimal erfolgreich.
- Kombi-PDF aus EN+DE mit `pypdf` neu gemergt.
- Vor Upload wurden die PDF-Metadaten bereinigt; die sichtbaren Textspuren bleiben unverändert und enthalten echte Umlaute.
- `_tools/check_refs.py`: EN 16/16 und DE 16/16, keine fehlenden oder unzitierten Keys. Die alphabetische Sortierwarnung ist bekannt und beabsichtigt, weil die Bibliographie thematisch gruppiert ist.
- Logscan: keine LaTeX-Fehler, keine undefinierten Referenzen/Zitate, keine Overfull-HBoxen; nur harmlose Underfulls, ein `h`→`ht`-Float-Hinweis in DE und MiKTeX-Updatehinweis.
- Deutsche PDF-Textspur enthält echte Umlaute, u. a. `Führer`, `Höhe`, `Einträge`, `für`.

## Artefakte

| Datei | Seiten | MD5 |
|---|---:|---|
| `DIRICHLET_CHARACTER_ATLAS_v2_en.pdf` | 13 | `E76FCCBBCA511E02A2094123753597A6` |
| `DIRICHLET_CHARACTER_ATLAS_v2_de.pdf` | 14 | `3E9704B96202430CCC3CBECEB42C423C` |
| `DIRICHLET_CHARACTER_ATLAS_v2_kombi.pdf` | 27 | `83A189DE0EE3D47C17E6286E11F2CAF3` |

## Quellenlinks

- Zenodo API Dirichlet Atlas live v1: https://zenodo.org/api/records/19960810
- Zenodo API Zeta Zoo Concept/latest: https://zenodo.org/api/records/19673226
- Zenodo API Zookeeper Concept/latest: https://zenodo.org/api/records/19673126
- Zenodo API RH Even-Dominance Concept/latest: https://zenodo.org/api/records/19035640
- Crossref Connes 1999: https://api.crossref.org/works/10.1007/s000290050042
- AMS Iwaniec--Kowalski: https://bookstore.ams.org/COLL/53
- AMS Paley--Wiener: https://bookstore.ams.org/coll-19
- CERN Weil 1952: https://cds.cern.ch/record/471308?ln=en
- Selberg 1956 journal page: https://www.i-scholar.in/index.php/JIMSIMS/article/view/146884
- LMFDB: https://www.lmfdb.org/
- LMFDB Dirichlet-L-function example page: https://www.lmfdb.org/L/Character/Dirichlet/5/4/

## Zenodo-Folge

Zenodo-v2-New-Version-Upload ist erforderlich, weil v1 bereits veröffentlicht ist und der lokale v2-Kandidat eine quellenbezogene Textkorrektur enthält.
