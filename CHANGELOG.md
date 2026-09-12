# Changelog

All notable changes to the Functional Stability Theory (FST) repository will be documented in this file.

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
