# Contributing to Functional Stability Theory / Mitwirken an der Funktionalen Stabilitätstheorie

[English](#english) | [Deutsch](#deutsch)

---

<a id="english"></a>
## English

Thank you for your interest in contributing to **Functional Stability Theory (FST)** (`research-line/functional-stability-theory`), a unified open-science mathematical research programme identifying *Functional Positivity under Gauge Constraint* (Pattern A) as the common foundation of open problems across number theory, mathematical physics, and cosmology.

### 1. Architectural Principles & 10 Governance Invariants

All contributions, numerical validation scripts, and theoretical supplements must strictly adhere to the 10 Governance and Research Invariants:

1. **100% Local-First & Zero-Egress (`INV-LOCAL-01`)**: All numerical validation scripts, test harnesses, and diagnostics execute completely offline. Zero network calls, telemetry, analytics, or remote cloud compute.
2. **Unprivileged Execution / RunAsInvoker (`INV-UNPRIV-02`)**: Execution operates strictly in unprivileged user-mode. Elevation to root or administrator privileges is never required or invoked.
3. **Deterministic Reproducibility (`INV-DETERM-03`)**: Mathematical algorithms and parameter sweeps yield bit-for-bit verifiable, reproducible results and machine-readable CSV/JSON ledgers across all supported platforms.
4. **Claim-Level Disambiguation (`INV-DISAMBIG-04`)**: Unconditional proven theorems (e.g., CCM microcluster closure) are strictly demarcated from conditional transfer hypotheses (e.g., Yang–Mills continuum limit) and open problems.
5. **Fail-Closed Ledger Gatekeeping (`INV-FAILCLOSED-05`)**: Verification scripts immediately fail-closed and reject circular certificates, degenerate scales, or ill-conditioned inputs rather than generating ambiguous approximations.
6. **Immutable Zenodo Anchors (`INV-ZENODO-06`)**: All milestone paper branches, preprints, and reproducible datasets are permanently bound to immutable Zenodo Concept-DOIs for verifiable scholarly provenance.
7. **Multi-OS Platform Parity (`INV-PLATFORM-07`)**: Core scripts and LaTeX document builds deliver identical mathematical logic across Windows, POSIX Linux, and macOS environments.
8. **Multi-Host & Cloud-Sync Hardening (`INV-SYNCHARD-08`)**: Rigorous `.gitignore` policies safeguard the working tree from cloud synchronization conflicts (`*-conflict-*`) and multi-agent lock contention (`LOCK.*`).
9. **Transparent Diagnostic Floor (`INV-FLOOR-09`)**: Unconditional positive controls are strictly segregated from conditional transfer hypotheses.
10. **48h Security Response & 5-Day Triage SLA (`INV-SLA-10`)**: Binding commitment to acknowledge security and integrity notifications within 48 hours and provide detailed triage within 5 business days via `security@open-bricks.org`, `security@ellmos.ai`, and `support@lukasgeiger.com`.

### 2. Plan D Local Development Workflow

In accordance with our cross-system architecture (Plan D), the local git clone at `C:\_Local_DEV\repos\functional-stability-theory` serves as the authoritative **Source of Truth**. Development, testing, and git operations take place exclusively within this local repository clone.

```bash
# Navigate to the canonical local clone
cd C:\_Local_DEV\repos\functional-stability-theory

# Verify git status and branch
git status
git branch --show-current

# Install testing dependencies
pip install pytest ruff

# Run the complete test suite
pytest

# Verify static linting and formatting
ruff check .
```

### 3. Version Freeze Discipline (`T-20260920-167562623`)

`functional-stability-theory` operates under strict version-freeze discipline (`T-20260920-167562623`). The version identifier (`1.0.6` in `pyproject.toml`) remains frozen. All ongoing refinements, hygiene enhancements, and discoverability updates are documented under `## [Unreleased]` in `CHANGELOG.md`.

### 4. Quality Gates

Before submitting contributions or pushing changes, verify all local quality gates:

1. `pytest`: 100% green test execution across all verification ledgers and contract tests.
2. `ruff check .`: Zero lint errors.
3. `python -m compileall -q .`: Zero bytecode compilation errors.
4. `git diff --check`: Zero whitespace anomalies or trailing whitespace errors.
5. `git diff -G"version = "`: Zero unauthorized version modifications.

### 5. Statutory Notice (§ 521 BGB) & Scholarly Licensing

This research repository is made available free of charge under the [Creative Commons Attribution 4.0 International (CC-BY-4.0)](LICENSE) license, with computational scripts licensed under permissive open-source terms. Under German statutory law (**§ 521 BGB** — *Haftung des Schenkers*), liability for gratuitously provided research assets is strictly limited to intentional misconduct (*Vorsatz*) and gross negligence (*grobe Fahrlässigkeit*).

---

<a id="deutsch"></a>
## Deutsch

Vielen Dank für Ihr Interesse an einer Mitarbeit an der **Funktionalen Stabilitätstheorie (FST)** (`research-line/functional-stability-theory`), einem vereinheitlichten mathematischen Open-Science-Forschungsprogramm, das die *Funktionale Positivität unter Eich-Constraint* (Muster A) als gemeinsames Substrat ungelöster Probleme in der Zahlentheorie, der mathematischen Physik und der Kosmologie identifiziert.

### 1. Architektur-Prinzipien & 10 Governance-Invarianten

Alle Beiträge, numerischen Validierungsskripte und theoretischen Ergänzungen müssen die 10 Governance- und Forschungsinvarianten strikt einhalten:

1. **100% Local-First & Zero-Egress (`INV-LOCAL-01`)**: Alle numerischen Validierungsskripte, Test-Suites und Diagnostiken laufen vollständig offline. Keine Netzwerkanfragen, keine Telemetrie, keine Analytik und keine externen Cloud-Ressourcen.
2. **Unprivilegierte Ausführung / RunAsInvoker (`INV-UNPRIV-02`)**: Die Ausführung erfolgt ausschließlich im unprivilegierten Benutzermodus. Root- oder Administratorrechte werden weder benötigt noch angefordert.
3. **Deterministische Reproduzierbarkeit (`INV-DETERM-03`)**: Mathematische Algorithmen und Parametersweeps liefern bitgenau überprüfbare, reproduzierbare Ergebnisse und maschinenlesbare CSV/JSON-Ledger über alle unterstützten Plattformen hinweg.
4. **Aussagen-Disambiguierung (`INV-DISAMBIG-04`)**: Vollständig bewiesene Theoreme (z. B. CCM-Mikrocluster-Abschluss) sind strikt von konditionalen Transfer-Hypothesen (z. B. Yang–Mills Kontinuumslimes) und offenen Problemen getrennt.
5. **Fail-Closed Ledger-Gatekeeping (`INV-FAILCLOSED-05`)**: Validierungsskripte terminieren bei zirkulären Zertifikaten, degenerierten Skalen oder unzulässigen Eingaben sofort fail-closed, anstatt ungenaue Näherungen zu erzeugen.
6. **Unveränderliche Zenodo-Anker (`INV-ZENODO-06`)**: Alle Paper-Zweige, Preprints und reproduzierbaren Datensätze sind dauerhaft an unveränderliche Zenodo-Konzept-DOIs gebunden, um wissenschaftliche Provenienz zu sichern.
7. **Plattform-Parität (`INV-PLATFORM-07`)**: Kernskripte und LaTeX-Kompilate liefern identische mathematische Logik unter Windows, POSIX Linux und macOS.
8. **Multi-Host- & Cloud-Sync-Härtung (`INV-SYNCHARD-08`)**: Eine strenge `.gitignore`-Richtlinie schützt den Arbeitsbaum vor Cloud-Synchronisationskonflikten (`*-conflict-*`) und Multi-Agenten Lock-Kollisionen (`LOCK.*`).
9. **Transparenter diagnostischer Mindeststandard (`INV-FLOOR-09`)**: Unbedingte Positivkontrollen sind strikt von konditionalen Transfer-Hypothesen isoliert.
10. **48h Sicherheits-SLA & 5-Tage-Triage (`INV-SLA-10`)**: Verbindliche Zusage zur Bestätigung von Sicherheits- und Integritätsmeldungen innerhalb von 48 Stunden und vollständige Triage innerhalb von 5 Werktagen via `security@open-bricks.org`, `security@ellmos.ai` und `support@lukasgeiger.com`.

### 2. Plan D Lokaler Entwicklungs-Workflow

Gemäß unserer systemweiten Architektur (Plan D) ist das lokale Git-Repository unter `C:\_Local_DEV\repos\functional-stability-theory` die verbindliche **Source of Truth**. Entwicklung, Tests und Commits finden ausschließlich im kanonischen lokalen Klon statt.

```bash
# In den kanonischen lokalen Klon wechseln
cd C:\_Local_DEV\repos\functional-stability-theory

# Git-Status und Branch prüfen
git status
git branch --show-current

# Test-Abhängigkeiten installieren
pip install pytest ruff

# Vollständige Testsuite ausführen
pytest

# Statische Code-Analyse ausführen
ruff check .
```

### 3. Version-Freeze-Disziplin (`T-20260920-167562623`)

`functional-stability-theory` unterliegt einer strikten Version-Freeze-Disziplin (`T-20260920-167562623`). Die Versionsnummer (`1.0.6` in `pyproject.toml`) bleibt unverändert eingefroren. Alle fortlaufenden Verfeinerungen, Hygieneanpassungen und Auffindbarkeitsverbesserungen werden unter `## [Unreleased]` in `CHANGELOG.md` dokumentiert.

### 4. Quality Gates

Vor dem Einreichen von Beiträgen müssen alle lokalen Quality Gates erfolgreich absolviert werden:

1. `pytest`: 100% grüne Testergebnisse über alle Validierungs-Ledger und Vertragstests.
2. `ruff check .`: 0 Linter-Fehler.
3. `python -m compileall -q .`: 0 Bytecode-Kompilierungsfehler.
4. `git diff --check`: 0 Whitespace- oder Zeilenumbruchfehler.
5. `git diff -G"version = "`: 0 unautorisierte Versionsmodifikationen.

### 5. Gesetzlicher Haftungshinweis (§ 521 BGB) & Wissenschaftliche Lizenzierung

Dieses Forschungsprojekt wird unentgeltlich unter der [Creative Commons Attribution 4.0 International (CC-BY-4.0)](LICENSE) Lizenz zur Verfügung gestellt; begleitende Rechenskripte unterliegen permissiven Open-Source-Lizenzen. Gemäß **§ 521 BGB** (*Haftung des Schenkers*) ist die Haftung für unentgeltlich überlassene wissenschaftliche Inhalte und Software auf Vorsatz und grobe Fahrlässigkeit beschränkt.
