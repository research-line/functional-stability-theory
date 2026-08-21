# Changelog

All notable changes to the Functional Stability Theory (FST) repository will be documented in this file.

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
