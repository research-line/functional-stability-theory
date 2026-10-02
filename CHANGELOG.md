# Changelog

All notable changes to the Functional Stability Theory (FST) repository will be documented in this file.

## [Unreleased]

### Added
- Bilingual [`CONTRIBUTING.md`](CONTRIBUTING.md) developer guidelines formalizing all 10 Governance and Research Invariants (`INV-LOCAL-01` to `INV-SLA-10`), Plan D local clone development workflow (`C:\_Local_DEV\repos\functional-stability-theory`), unprivileged `RunAsInvoker` user-mode execution, strict version-freeze discipline (`T-20260920-167562623`), § 521 BGB statutory liability waiver, and 48-hour security response SLA.
- Level 1 SBOM plain-text companion [`THIRD_PARTY_LICENSES.txt`](THIRD_PARTY_LICENSES.txt) Stand 2026-10-03 with comprehensive runtime, scientific, and testing component inventory, zero-copyleft certification, and invariant compliance ledger.
- Standardized PEP 621 metadata URLs in `pyproject.toml` registering `Contributing`, `Third-Party Licenses (Text)`, `Level 1 SBOM`, `Level 1 SBOM (Text)`, and `Plain-Text License`, and added `THIRD_PARTY_LICENSES.txt` to `license-files` whitelist.
- Shields.io badges in `README.md` and `README_de.md` for Contributing Guidelines (`Contributing-Guidelines-blue.svg` / `Mitwirken-Leitfaden-blue.svg`), Level 1 SBOM Plain Text (`SBOM-Level_1_Text-informational.svg`), updated test suite count (`Tests-146+-Passed`), and refreshed audit timestamp (`2026-10-03`).
- Documented turnus Pfad B Discoverability, Level 1 SBOM text companion, and bilingual developer guidelines in `MARKETING-LOG.txt`.
- Expanded automated metadata contract test suite in `tests/test_metadata.py` asserting `CONTRIBUTING.md` presence and bilingual parity, `THIRD_PARTY_LICENSES.txt` structure, PEP 621 SBOM and Contributing URLs, badge parity, and audit currency (`2026-10-03`).
- Canonical [`NOTICE`](NOTICE) attribution file formalizing copyright Lukas Geiger, research-line parent organization, and open-bricks umbrella under CC-BY-4.0.
- Standardized PEP 621 metadata in `pyproject.toml`:
  - Saturated 20/20 `keywords` in exact parity with GitHub repository topics.
  - Included `NOTICE` in `license-files` whitelist (`["LICENSE", "NOTICE", "THIRD_PARTY_LICENSES.md"]`).
  - Added canonical `Notice` entry to `[project.urls]`.
  - Configured `addopts = "-ra -v --basetemp=.pytest_temp"` and added `.pytest_temp` and `.hypothesis` to `norecursedirs`.
- Automated contract test suite `tests/test_compute_turbulence.py` covering turbulence paper artifact presence, dual DFC1 scenario evaluations across canonical controls, phi profile bounds, and flux profile properties.
- Hardened `.gitignore` with `.pytest_temp/`, `.pytest_tmp*/`, `*-IDEAPAD*`, `.automation-lock`, and `Desktop.ini`.
- Appended Pfad A technical hygiene and metadata contract audit entry to `MARKETING-LOG.txt`.
- Expanded metadata contract tests in `tests/test_metadata.py` to assert NOTICE file presence, PEP 621 20/20 keywords parity, Notice URL, pytest options, and .gitignore hardening.
- Complete 18-point bilingual navigation parity across `README.md` and `README_de.md` with reciprocal HTML anchors (`<a id="..."></a>`).
- Elevated Theoretical Data Flow & Validation Sequence to Section 07 in both language surfaces.
- Integrated Section 16: Comparative Matrix vs. Alternatives (10 invariant dimensions across FST, Classical Number Theory, Connes Noncommutative Geometry, Traditional CFD Simulation, and Closed Math Suites).
- Integrated Section 18: Security Policy, Author & Statutory Notice (§ 521 BGB Gefälligkeitsrecht for gratuitous academic research software).
- Synchronized ASCII Proof Architecture overview (Section 10), Independent Foundations (Section 11), and Repository Structure directory tree (Section 13) in `README_de.md`.

### Changed
- Re-audited `THIRD_PARTY_LICENSES.md` to `2026-10-03` with cross-reference to plain-text Level 1 SBOM companion [`THIRD_PARTY_LICENSES.txt`](THIRD_PARTY_LICENSES.txt) and developer guidelines [`CONTRIBUTING.md`](CONTRIBUTING.md).
- Updated `llms.txt` verification timestamp to `2026-10-03` with 146+ passed tests baseline and links to `CONTRIBUTING.md` and `THIRD_PARTY_LICENSES.txt`.
- Synchronized Shields.io badges in `README.md` and `README_de.md` (`Contributing Guidelines` / `Mitwirken-Leitfaden`, `Level 1 SBOM Text`, `Tests 146+ Passed` / `146+ Bestanden`, `Last-checked 2026-10-03` / `Geprüft 2026-10-03`).
- Synchronized Turbulence Cascade companion paper artifacts in `fst-physics/turbulence/` (v1.8 maintenance release: isolated 3-page frontmatter architecture with table of contents on separate page, full bibliography Underfull/Badness-10000 elimination via `xurl` and `\doilink`, standardized bilingual AI disclosures `ai_disclosure_STANDARD_{en,de}.tex`, and verified 71-page combined bilingual PDF compilation; 35 pages EN, 36 pages DE, 71 pages kombi PDF; 0 replacement characters).
- Synchronized BSD Positivity domain supplement artifacts in `fst-mathematics/bsd/` (v1.5 maintenance release: 3-page frontmatter architecture, Chicago/APA title casing, author-pair en-dash typography, comprehensive RevTeX 4-2 table hardening across Tables 1–8, and 4-tier TikZ vector architecture schema `fig:bsd_architecture`; 29 pages EN, 31 pages DE, 60 pages kombi PDF; 0 replacement characters).
- Added automated contract unit test suite `tests/test_compute_bsd.py` covering BSD formula verification across 4 reference curves, regulator positivity (Axiom II), and Cremona database rank-2 regulator sampling.
- Re-audited `THIRD_PARTY_LICENSES.md` to `2026-09-26` with explicit cross-reference to [`NOTICE`](NOTICE) and Level 1 SBOM invariant verification (100% permissive open-source dependencies, 0% copyleft, 100% offline zero-egress).
- Synchronized Shields.io badges in `README.md` and `README_de.md` (`Attribution: NOTICE`, `Last-checked 2026-09-26` / `Geprüft 2026-09-26`).
- Updated `llms.txt` verification timestamp to `2026-09-26` with canonical NOTICE attribution entrypoint and references to the 18-point bilingual navigation table and comparative matrix.

## [1.0.6] - 2026-09-18

### Added
- Hardened GitHub Actions CI/CD workflows:
  - Top-level least-privilege `permissions: contents: read` and job `timeout-minutes: 15` in `.github/workflows/ci.yml`.
  - Concurrency control (`cancel-in-progress: true`) and job `timeout-minutes: 10` in `.github/workflows/stale.yml`.
  - Upgraded `actions/first-interaction` to `v3` with concurrency control and job `timeout-minutes: 5` in `.github/workflows/welcome.yml`.
- Extended multi-host cloud-sync defense patterns in `.gitignore` (`*conflicted copy*`, `* (Kopie)*`, `* (Copy)*`, `*-ASUS*`, `*-ASUS-GEI*`, `*-LAPTOP*`, `*-WORKSTATION*`, `*-WORKSTATION-LG*`, `*-Mac Studio*`, `*-MacBook*`).
- Extended canonical multi-agent lock patterns in `.gitignore` (`LOCK`, `LOCK.*`, `LOCK*.txt`, `LOCK.permissions.json`, `LOCK.user.*`, `LOCK.until.*`, `LOCK.condition.*`, `uv.lock`, with explicit negation `!package-lock.json`).
- Extended build, cache, and test defense patterns in `.gitignore` (`.hypothesis/`, `.turbo/`, `.nyc_output/`, `*.orig`, `*.rej`).
- Declared PEP 621 standard `license-files = ["LICENSE", "THIRD_PARTY_LICENSES.md"]` in `pyproject.toml`.
- Configured pytest execution guardrails in `pyproject.toml`: `minversion = "7.0"` and `norecursedirs` excluding all non-test domain, script, and build directories.
- Documented turnus Pfad A technical hygiene and CI hardening in `MARKETING-LOG.txt`.
- Expanded metadata contract tests in `tests/test_metadata.py` covering least-privilege CI permissions, workflow timeouts, action versions, extended `.gitignore` patterns, pytest guardrails, PEP 621 license files, and audit records.

### Changed
- Bumped project version to `1.0.6` across `pyproject.toml`, test suites, and documentation.
- Re-audited `THIRD_PARTY_LICENSES.md` to `2026-09-18` (v1.0.6, 100% permissive, 0% copyleft, 100% offline zero-egress).
- Synchronized Shields.io badges in `README.md` and `README_de.md` (`v1.0.6`, `Tests 128+ Passed` / `128+ Bestanden`, `Last-checked 2026-09-18` / `Geprüft 2026-09-18`).
- Updated `llms.txt` verification timestamp to `2026-09-18` with updated test suite count (128+ Passed) and release version `1.0.6`.

## [1.0.5] - 2026-09-12

### Added
- Audited software inventory in `THIRD_PARTY_LICENSES.md` certifying 100% permissive open-source dependencies (CC-BY-4.0, BSD-3-Clause, MIT, PSF-2.0), 0% copyleft risk (zero AGPL/GPL), 100% offline zero-egress compliance, and full alignment with 10 Governance Invariants (`INV-LOCAL-01` to `INV-SLA-10`).
- Target Personas & Value Propositions (4 personas: Field Theorists, Millennium Researchers, Open-Science Curators, AI Research Agents) in `README.md` and `README_de.md`.
- Expanded Quick Navigation to 16 points in `README.md` and 14 points in `README_de.md` with complete reciprocal anchor parity (`#target-personas--discoverability` / `#zielgruppen--auffindbarkeit`, `#third-party-licenses--transparency` / `#drittanbieter-lizenzen--transparenz`).
- Dedicated Third-Party Licenses & Transparency sections in both `README.md` and `README_de.md`.
- PEP 621 metadata URLs in `pyproject.toml` for `Third-Party Licenses`, `Marketing Log`, and `LLM Ready`.
- Comprehensive Pfad B audit entry in `MARKETING-LOG.txt` featuring a 10-dimension 5-way competitive matrix, bilingual high-intent keyword search queries, and historical release tracking.
- Expanded contract test suite in `tests/test_metadata.py` verifying `THIRD_PARTY_LICENSES.md` structure, 16-point quick navigation, target personas, third-party license disclosures, and PEP 621 URLs.

### Changed
- Bumped project version to `1.0.5` across `pyproject.toml`, test suites, and documentation.
- Synchronized Shields.io badges in `README.md` and `README_de.md` (`v1.0.5`, `Tests 118+ Passed` / `118+ Bestanden`, `Third-Party Licenses: Audited`, `Marketing Log: Active`, `Last-checked 2026-09-12`).
- Updated `llms.txt` verification timestamp to `2026-09-12` with updated test suite count (118+ Passed) and `THIRD_PARTY_LICENSES.md` reference.

## [1.0.4] - 2026-09-10

### Added
- Standardized pytest execution flags (`addopts = "-ra -v"`) in `pyproject.toml` and GitHub Actions CI runner (`.github/workflows/ci.yml`).
- Extended `.gitignore` multi-host protection patterns (`*-ASUS-GEI.*`, `*-WORKSTATION-LG.*`, `*-WORKSTATION.*`, `* (kopie)*`, `* (copy)*`), multi-agent lock files (`LOCK.permissions.json`, `uv.lock`), and packaging caches (`wheelhouse/`, `.wheel-smoke/`).
- Automated contract test expansion in `tests/test_metadata.py` verifying standard pytest CLI flags, CI workflow test step configuration, extended gitignore rules, and Pfad A changelog entry.

### Changed
- Bumped project version to `1.0.4` across `pyproject.toml`, test suites, and documentation.
- Synchronized Shields.io badges in `README.md` and `README_de.md` (`v1.0.4`, test suite status, `Last-checked 2026-09-10`).
- Updated `llms.txt` verification timestamp to `2026-09-10` with updated test suite count and release version.

## [1.0.3] - 2026-09-09

### Added
- 14-point quick navigation index in both `README.md` and `README_de.md` with full bidirectional anchor parity.
- Auto-numbered Mermaid sequence diagram (`sequenceDiagram`) detailing the reproducible Zero-Egress Numerical Validation & Contract Audit Lifecycle.
- 10-point Governance & Runtime Invariants table in both `README.md` and `README_de.md` covering local-first execution, non-elevation, claim-level disambiguation, and fail-closed gatekeeping.
- Supported versions table (`1.0.x`), 48-hour response SLA, and binding 5 business days triage commitment in `SECURITY.md` (bilingual English & German).
- Umbrella security contact (`security@open-bricks.org`) alongside existing disclosure routes.
- GitHub Actions CI concurrency control (`cancel-in-progress: true`) and bytecode compilation gate (`python -m compileall -q .`).
- `.gitignore` hardening against multi-host conflict copies (`*-conflict-*`, `*.sync-temp-*`, etc.) and multi-agent locks (`LOCK.*`, `*.lock`, `LOCK*.txt`).
- Local `MARKETING-LOG.txt` documenting Pfad B discoverability, SEO badges, and contract testing metrics.
- Comprehensive contract test suite expansion in `tests/test_metadata.py` verifying all metadata invariants and security commitments.

### Changed
- Bumped project version to `1.0.3` in `pyproject.toml` and updated project URLs to include `Parent Organization` and `Umbrella Ecosystem`.
- Synchronized Shields.io badges in `README.md` and `README_de.md` (v1.0.3, 110+ passed tests, 48h Security SLA, RunAsInvoker, Last-checked 2026-09-09).
- Updated `llms.txt` verification timestamp to 2026-09-09 with updated test suite count and invariants summary.

## [1.0.2] - 2026-08-23

### Added
- Multi-OS (`ubuntu-latest`, `windows-latest`, `macos-latest`) and multi-Python (`3.10`, `3.11`, `3.12`, `3.13`) GitHub Actions CI workflow in `.github/workflows/ci.yml`.
- PEP 621 Operating System (`OS Independent`, `Microsoft Windows`, `POSIX Linux`, `MacOS`) and Python 3.13 classifiers in `pyproject.toml`.
- Standard project URLs (`Issues`, `Changelog`, `Security`) in `pyproject.toml`.
- CI workflow integrity and PEP 621 contract validation test cases in `tests/test_metadata.py` expanding test suite to 12/12 passing tests.
- CI status badges in `README.md` and `README_de.md`.

### Changed
- Updated `llms.txt` verification timestamp to 2026-08-23 and updated test suite / CI matrix specifications.
- Updated test suite badges in `README.md` and `README_de.md` to reflect 12/12 passing tests across Python 3.10–3.13.

## [1.0.1] - 2026-08-21

### Added
- Bilingual `SECURITY.md` defining Open-Science Research Integrity, 100% offline and zero-egress guarantees for numerical scripts, and coordinated vulnerability disclosure channels (`security@ellmos.ai`).
- Interactive Mermaid sequence diagram `Theoretical Data Flow & Validation Sequence` illustrating the progression from axioms (RFEP / Pattern A / DS1-DS3) to Master proofs, domain instantiations, and local numerical diagnostics in `README.md` and `README_de.md`.
- Expanded sibling research and toolchain ecosystem matrix across `research-line`, `biotec-line`, `doc-bricks`, `dev-bricks`, `ellmos-ai`, and `open-bricks`.
- Security policy badges and updated 10/10 test suite metrics in `README.md` and `README_de.md`.
- Automated security policy validation and Ruff lint compliance tests in `tests/test_metadata.py`.
- 5 new Yang-Mills verification and transfer scripts in `scripts/yang-mills/`: `compute_rp_os_transfer_ledger.py`, `compute_rp_os_rfep_transfer_ledger.py`, `compute_u1_2d_strong_coupling_positive_control.py`, `compute_ym_waisen_transfer_ledger.py`, and `u1_strong_coupling_positive_control.py`.

### Changed
- Synchronized Yang-Mills domain preprint (`fst-physics/yang-mills/`) to latest LaTeX design/layout standards (frontmatter isolation, abstract roman p.1, TOC p.2, arabic main body p.1, 2-pass clean pdflatex build, EN 48 S., DE 49 S., Kombi 97 S.).
- Synchronized P vs NP domain preprint (`fst-mathematics/p-vs-np/`) with English style and US typography normalization.
- Updated `fst-physics/yang-mills/README.md`, `README.md`, and `README_de.md` numerical validation tables and reproduction runbooks.
- Configured `[tool.ruff]` in `pyproject.toml` with clean package exclusions for standalone numerical research scripts, achieving 100% clean repository-wide linting (`ruff check .` 0 errors).
- Updated `llms.txt` verification timestamp to 2026-08-21.

## [1.0.0] - 2026-08-16

### Added
- Standardized `pyproject.toml` configuration with PEP 621 metadata, pytest configuration targeting `tests/`, and Ruff lint rules.
- Automated metadata, manifest, and parity test suite in `tests/test_metadata.py` (8/8 passed).
- Test suite and Python version badges in both `README.md` and `README_de.md`.
- Comprehensive bilingual research sibling and ecosystem matrix linking `research-line`, `ellmos-ai`, `dev-bricks`, and `open-bricks`.
- Synchronized domain supplement version tables and numerical validation script tables between `README.md` and `README_de.md`.

### Changed
- Updated `llms.txt` verification timestamp to 2026-08-16 and added test suite context.

## [0.1.2] - 2026-08-04

### Added
- German language documentation `README_de.md` for bilingual accessibility.
- Organization (`research-line`) and umbrella (`open-bricks`) ecosystem badges in `README.md` & `README_de.md`.
- Language switcher between English (`README.md`) and German (`README_de.md`).

### Changed
- Updated `llms.txt` verification timestamp to 2026-08-04.

## [0.1.1] - 2026-07-27

### Added
- Interactive Mermaid architecture diagram for FST Master Foundations and Domain Supplements in `README.md`.

### Changed
- Updated `llms.txt` verification timestamp to 2026-07-27.

## [0.1.0] - 2026-07-25

### Added
- Shields.io metadata badges (License CC-BY 4.0, ORCID, Zenodo Concept-DOI, LLM context) in `README.md`.
- Machine-readable AI / LLM integration callout (`> [!NOTE]`) in `README.md`.

### Changed
- Updated `llms.txt` index verification date to 2026-07-25.
