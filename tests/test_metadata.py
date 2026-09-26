"""Metadata and discoverability parity tests for Functional Stability Theory (FST)."""

import pathlib
try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 compatibility
    import tomli as tomllib


REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent


def test_pyproject_metadata():
    """Validate pyproject.toml PEP 621 metadata."""
    pyproject_path = REPO_ROOT / "pyproject.toml"
    assert pyproject_path.is_file(), "pyproject.toml must exist"

    with open(pyproject_path, "rb") as f:
        data = tomllib.load(f)

    project = data.get("project", {})
    assert project.get("name") == "functional-stability-theory"
    assert project.get("version") == "1.0.6"
    assert project.get("requires-python") == ">=3.10"
    assert project.get("license", {}).get("text") == "CC-BY-4.0"
    assert len(project.get("authors", [])) >= 1
    assert project["authors"][0]["name"] == "Lukas Geiger"


def test_pyproject_pep621_classifiers_and_urls():
    """Validate PEP 621 classifiers and project URLs in pyproject.toml."""
    pyproject_path = REPO_ROOT / "pyproject.toml"
    with open(pyproject_path, "rb") as f:
        data = tomllib.load(f)

    project = data.get("project", {})
    classifiers = project.get("classifiers", [])
    assert "Operating System :: OS Independent" in classifiers
    assert "Operating System :: Microsoft :: Windows" in classifiers
    assert "Operating System :: POSIX :: Linux" in classifiers
    assert "Operating System :: MacOS" in classifiers
    assert "Programming Language :: Python :: 3.13" in classifiers

    urls = project.get("urls", {})
    assert "Homepage" in urls
    assert "Documentation" in urls
    assert "Repository" in urls
    assert "Issues" in urls
    assert "Changelog" in urls
    assert "Security" in urls
    assert "Parent Organization" in urls
    assert "Umbrella Ecosystem" in urls
    assert "Third-Party Licenses" in urls
    assert "Notice" in urls
    assert "Marketing Log" in urls
    assert "LLM Ready" in urls


def test_ci_workflow_integrity():
    """Validate GitHub Actions CI matrix workflow configuration."""
    ci_path = REPO_ROOT / ".github" / "workflows" / "ci.yml"
    assert ci_path.is_file(), ".github/workflows/ci.yml must exist"

    content = ci_path.read_text(encoding="utf-8")
    assert "ubuntu-latest" in content
    assert "windows-latest" in content
    assert "macos-latest" in content
    assert "3.10" in content
    assert "3.11" in content
    assert "3.12" in content
    assert "3.13" in content
    assert "ruff check ." in content
    assert "pytest" in content
    assert "concurrency:" in content
    assert "cancel-in-progress: true" in content
    assert "python -m compileall -q ." in content


def test_readme_and_readme_de_badges():
    """Validate badge parity in README.md and README_de.md."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    # Required badge signatures in English README
    assert "License-CC_BY_4.0" in readme_en
    assert "version-1.0.6" in readme_en
    assert "actions/workflows/ci.yml" in readme_en
    assert "ORCID-0009--0005--7296--1534" in readme_en
    assert "Zenodo-10.5281" in readme_en
    assert "Security-Research%20Integrity" in readme_en
    assert "Security-RunAsInvoker" in readme_en
    assert "Security%20SLA" in readme_en
    assert "Ecosystem-research--line" in readme_en
    assert "Umbrella-open--bricks" in readme_en
    assert "Attribution-NOTICE" in readme_en
    assert "Third--Party%20Licenses-Audited-brightgreen" in readme_en
    assert "Marketing%20Log-Active-blue" in readme_en
    assert "llms.txt" in readme_en
    assert "Tests-128%2B%20Passed" in readme_en
    assert "Python-3.10--3.13" in readme_en or "Python-3.10" in readme_en
    assert "Platform-Windows" in readme_en
    assert ("Audit-Last--checked%202026--09--26" in readme_en or "Audit-Last--checked%202026--09--20" in readme_en)

    # Required badge signatures in German README
    assert "Lizenz-CC_BY_4.0" in readme_de or "License-CC_BY_4.0" in readme_de
    assert "Attribution-NOTICE" in readme_de
    assert "Version-1.0.6" in readme_de
    assert "actions/workflows/ci.yml" in readme_de
    assert "ORCID-0009--0005--7296--1534" in readme_de
    assert "Zenodo-10.5281" in readme_de
    assert "Sicherheit-Open%20Science%20Integrit%C3%A4t" in readme_de
    assert "Sicherheit-RunAsInvoker" in readme_de
    assert "Sicherheits--SLA" in readme_de
    assert "research--line" in readme_de
    assert "open--bricks" in readme_de
    assert "Drittanbieter--Lizenzen-Gepr%C3%BCft-brightgreen" in readme_de
    assert "Marketing--Log-Aktiv-blue" in readme_de
    assert "llms.txt" in readme_de
    assert "Tests-128%2B%20Bestanden" in readme_de
    assert "Python-3.10--3.13" in readme_de or "Python-3.10" in readme_de
    assert "Plattform-Windows" in readme_de
    assert ("Audit-Gepr%C3%BCft%202026--09--26" in readme_de or "Audit-Gepr%C3%BCft%202026--09--20" in readme_de)


def test_quick_navigation_parity():
    """Validate quick navigation index parity in English and German (18-point) READMEs."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    assert "## Quick Navigation" in readme_en
    assert "## Schnellnavigation" in readme_de

    en_points = [
        "01. Start Here",
        "02. Discovery Context",
        "03. The Five Masters",
        "04. Domain Supplements",
        "05. Glossary",
        "06. Proof Architecture",
        "07. Theoretical Data Flow",
        "08. Numerical Lifecycle",
        "09. Governance Invariants",
        "10. ASCII Architecture",
        "11. Independent Foundations",
        "12. Validation Scripts",
        "13. Repository Structure",
        "14. Sibling Ecosystem",
        "15. Target Personas",
        "16. Comparative Matrix",
        "17. Third-Party Licenses",
        "18. Security Policy & Statutory Notice",
    ]
    for point in en_points:
        assert point in readme_en, f"Quick nav point {point} missing in README.md"

    de_points = [
        "01. Einstieg",
        "02. Auffindbarkeit",
        "03. Die Fünf Master-Arbeiten",
        "04. Domänen-Ergänzungen",
        "05. Glossar",
        "06. Beweisarchitektur",
        "07. Theoretischer Datenfluss",
        "08. Numerischer Lebenszyklus",
        "09. Governance-Invarianten",
        "10. ASCII-Übersicht",
        "11. Unabhängige Fundamente",
        "12. Validierungsskripte",
        "13. Repository-Struktur",
        "14. Ökosystem & Repositories",
        "15. Zielgruppen & Auffindbarkeit",
        "16. Vergleichsmatrix",
        "17. Drittanbieter-Lizenzen",
        "18. Sicherheitsrichtlinie & Haftung",
    ]
    for point in de_points:
        assert point in readme_de, f"Schnellnavigation point {point} missing in README_de.md"


def test_governance_invariants_table_parity():
    """Validate 10-point governance and runtime invariants table in both READMEs."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    assert "## Governance & Runtime Invariants" in readme_en
    assert "## Governance- & Laufzeit-Invarianten" in readme_de or "### Governance- & Laufzeit-Invarianten" in readme_de

    # 10 numbered rows
    for i in range(1, 11):
        num_str = f"**{i:02d}**"
        assert num_str in readme_en, f"Invariant {num_str} missing in README.md"
        assert num_str in readme_de, f"Invariant {num_str} missing in README_de.md"

    assert "100% Local-First & Zero-Egress" in readme_en
    assert "Unprivileged Execution (RunAsInvoker)" in readme_en
    assert "Deterministic Reproducibility" in readme_en
    assert "48h Security & 5-Day Triage SLA" in readme_en

    assert "100% Local-First & Zero-Egress" in readme_de
    assert "Rechtefreie Ausführung (RunAsInvoker)" in readme_de
    assert "Deterministische Reproduzierbarkeit" in readme_de
    assert "48h Sicherheits- & 5-Tage-Triage-SLA" in readme_de


def test_five_masters_structure():
    """Validate that all 5 Core Masters exist in the repository."""
    masters_dir = REPO_ROOT / "masters"
    assert masters_dir.is_dir(), "masters/ directory must exist"

    expected_masters = [
        "zookeeper",
        "zeta-zoo",
        "spectrum-duality",
        "atlas",
        "selberg",
    ]
    for master in expected_masters:
        master_path = masters_dir / master
        assert master_path.is_dir(), f"Master {master} must exist in masters/"
        assert (master_path / "README.md").is_file(), f"Master {master} must have a README.md"


def test_domain_supplements_structure():
    """Validate that domain supplement directories and script hubs exist."""
    expected_dirs = [
        "fst-mathematics",
        "fst-physics",
        "fst-cosmology",
        "fst-biology",
        "applications",
        "scripts",
    ]
    for domain in expected_dirs:
        domain_path = REPO_ROOT / domain
        assert domain_path.is_dir(), f"Domain supplement or script folder {domain} must exist"


def test_concept_dois_consistency():
    """Validate Concept-DOIs consistency across README.md, README_de.md, and llms.txt."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")
    llms_txt = (REPO_ROOT / "llms.txt").read_text(encoding="utf-8")

    key_dois = [
        "10.5281/zenodo.19673126",  # Zookeeper
        "10.5281/zenodo.19673226",  # Zeta Zoo
        "10.5281/zenodo.19036190",  # Spectrum Duality / RFEP
        "10.5281/zenodo.19960809",  # Atlas
        "10.5281/zenodo.19962588",  # Selberg
    ]

    for doi in key_dois:
        assert doi in readme_en, f"DOI {doi} missing from README.md"
        assert doi in readme_de, f"DOI {doi} missing from README_de.md"
        assert doi in llms_txt, f"DOI {doi} missing from llms.txt"


def test_llms_txt_structure_and_timestamp():
    """Validate llms.txt structure, canonical metadata, and verification date."""
    llms_txt = (REPO_ROOT / "llms.txt").read_text(encoding="utf-8")

    assert ("## Last-checked: 2026-09-26" in llms_txt or "## Last-checked: 2026-09-20" in llms_txt)
    assert "https://github.com/research-line/functional-stability-theory" in llms_txt
    assert "research-line" in llms_txt
    assert "Renormalized Free-Energy Principle" in llms_txt
    assert "Pattern A" in llms_txt
    assert "Version: 1.0.6" in llms_txt
    assert "NOTICE" in llms_txt
    assert "SECURITY.md" in llms_txt
    assert "## Search Phrases" in llms_txt
    assert "128+ Passed" in llms_txt
    assert "Governance & Runtime Invariants" in llms_txt
    assert "THIRD_PARTY_LICENSES.md" in llms_txt
    assert "MARKETING-LOG.txt" in llms_txt


def test_changelog_entry():
    """Validate that CHANGELOG.md documents the latest release and hygiene audits."""
    changelog = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert "2026-09-18" in changelog
    assert "1.0.6" in changelog
    assert "THIRD_PARTY_LICENSES.md" in changelog
    assert "SECURITY.md" in changelog
    assert "MARKETING-LOG.txt" in changelog


def test_mermaid_diagrams_parity():
    """Validate that both English and German READMEs contain all Mermaid diagrams."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    for doc, name in [(readme_en, "README.md"), (readme_de, "README_de.md")]:
        assert "```mermaid" in doc, f"{name} missing mermaid codeblock"
        assert "flowchart TD" in doc, f"{name} missing flowchart TD"
        assert "flowchart LR" in doc, f"{name} missing flowchart LR"
        assert "sequenceDiagram" in doc, f"{name} missing sequenceDiagram"
        assert "autonumber" in doc, f"{name} missing autonumber"
        assert "subgraph HYP" in doc, f"{name} missing subgraph HYP"


def test_security_policy_bilingual_parity():
    """Validate SECURITY.md exists and contains bilingual invariants and contact info."""
    security_file = REPO_ROOT / "SECURITY.md"
    assert security_file.is_file(), "SECURITY.md must exist"

    content = security_file.read_text(encoding="utf-8")
    assert "## English" in content
    assert "## Deutsche Fassung" in content
    assert "Zero-Egress" in content or "0% Netzwerk-Egress" in content
    assert "security@open-bricks.org" in content
    assert "security@ellmos.ai" in content
    assert "support@lukasgeiger.com" in content
    assert "48" in content
    assert "5" in content


def test_sibling_matrix_parity():
    """Validate that sibling and ecosystem repositories are linked in both READMEs."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    expected_siblings = [
        "research-line/fst-nash",
        "research-line/rh-even-dominance",
        "research-line/crm-cosmology",
        "research-line/prompt-archaeology-casestudy2",
        "research-line/ai-elite-swr",
        "research-line/economic-sanctions-coercive-diplomacy",
        "biotec-line/VFDistiller",
        "doc-bricks/MediaBrain",
        "dev-bricks/CodeBox",
        "dev-bricks/DevCenter",
        "dev-bricks/githubbot",
        "ellmos-ai/skills",
        "ellmos-ai/sqlite-transit-sync",
        "ellmos-ai/ellmos-development-system",
        "open-bricks/governance",
        "open-bricks",
    ]

    for sibling in expected_siblings:
        assert sibling in readme_en, f"Sibling {sibling} missing in README.md"
        assert sibling in readme_de, f"Sibling {sibling} missing in README_de.md"


def test_gitignore_hygiene():
    """Validate .gitignore hardening against conflict copies, lock files, and cache."""
    gitignore = (REPO_ROOT / ".gitignore").read_text(encoding="utf-8")
    assert "*-conflict-*" in gitignore
    assert "*.sync-temp-*" in gitignore
    assert "LOCK.*" in gitignore
    assert "*.lock" in gitignore
    assert ".pytest_cache/" in gitignore
    assert ".ruff_cache/" in gitignore


def test_marketing_log_exists():
    """Validate that local MARKETING-LOG.txt exists and contains Pfad B audit metadata."""
    log_file = REPO_ROOT / "MARKETING-LOG.txt"
    assert log_file.is_file(), "MARKETING-LOG.txt must exist"

    content = log_file.read_text(encoding="utf-8")
    assert "Pfad B" in content
    assert "1.0.3" in content
    assert "research-line/functional-stability-theory" in content


def test_gitignore_hygiene_patterns():
    """Validate comprehensive .gitignore patterns for multi-host, lock, cache, and editor hygiene."""
    gitignore = (REPO_ROOT / ".gitignore").read_text(encoding="utf-8")

    # Multi-host sync conflict patterns
    assert "*-conflict-*" in gitignore
    assert "*.sync-temp-*" in gitignore
    assert "*.sync-conflict-*" in gitignore
    assert "*.conflict" in gitignore
    assert "*-CONFLIT-*" in gitignore
    assert "*-ASUS-GEI.*" in gitignore
    assert "*-WORKSTATION-LG.*" in gitignore
    assert "*-WORKSTATION.*" in gitignore
    assert "* (kopie)*" in gitignore
    assert "* (copy)*" in gitignore

    # Multi-agent lock patterns
    assert "LOCK.*" in gitignore
    assert "*.lock" in gitignore
    assert "LOCK*.txt" in gitignore
    assert "LOCK" in gitignore
    assert "LOCK.permissions.json" in gitignore
    assert "uv.lock" in gitignore

    # Packaging, build, and test cache patterns
    assert ".pytest_cache/" in gitignore
    assert ".ruff_cache/" in gitignore
    assert ".coverage" in gitignore
    assert "coverage/" in gitignore
    assert "htmlcov/" in gitignore
    assert "wheelhouse/" in gitignore
    assert ".wheel-smoke/" in gitignore

    # Editor and temp patterns
    assert "*.tmp" in gitignore
    assert "*.bak" in gitignore
    assert "*.swp" in gitignore
    assert "*~" in gitignore
    assert "*.log" in gitignore


def test_pytest_configuration_and_flags():
    """Validate standard pytest configuration and addopts flags in pyproject.toml."""
    pyproject_path = REPO_ROOT / "pyproject.toml"
    with open(pyproject_path, "rb") as f:
        data = tomllib.load(f)

    pytest_opts = data.get("tool", {}).get("pytest", {}).get("ini_options", {})
    assert pytest_opts.get("testpaths") == ["tests"]
    assert pytest_opts.get("python_files") == ["test_*.py"]
    assert pytest_opts.get("python_functions") == ["test_*"]
    assert "-ra -v" in pytest_opts.get("addopts", "")


def test_ci_workflow_pytest_flags():
    """Validate that CI workflow executes pytest with standard -ra -v flags."""
    ci_path = REPO_ROOT / ".github" / "workflows" / "ci.yml"
    content = ci_path.read_text(encoding="utf-8")
    assert "pytest -ra -v" in content
    assert "python -m compileall -q ." in content
    assert "cancel-in-progress: true" in content


def test_changelog_recent_pfad_a_entry():
    """Validate that CHANGELOG.md contains the latest Pfad A release entry."""
    changelog = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert "## [1.0.4] - 2026-09-10" in changelog
    assert "addopts = \"-ra -v\"" in changelog
    assert ".gitignore" in changelog


def test_third_party_licenses_audit():
    """Validate that THIRD_PARTY_LICENSES.md exists and contains permissive inventory and invariants."""
    licenses_file = REPO_ROOT / "THIRD_PARTY_LICENSES.md"
    assert licenses_file.is_file(), "THIRD_PARTY_LICENSES.md must exist"

    content = licenses_file.read_text(encoding="utf-8")
    assert "CC-BY-4.0" in content
    assert "PSF-2.0" in content
    assert "BSD-3-Clause" in content
    assert "MIT" in content
    assert "Zero AGPL" in content or "0 AGPL" in content
    assert "Zero-Egress" in content
    assert "INV-LOCAL-01" in content
    assert "INV-SLA-10" in content
    assert "mpmath" in content
    assert "sympy" in content
    assert "pytest" in content
    assert "ruff" in content


def test_target_personas_sections():
    """Validate that target personas are documented in English and German READMEs."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    assert "## Target Personas & Discoverability" in readme_en
    assert "Theoretical Physicists & Field Theorists" in readme_en
    assert "Analytic Number Theorists & Millennium Researchers" in readme_en
    assert "Open-Science Curators & Formal Verification Reviewers" in readme_en
    assert "AI Research Agents & Literature Synthesizers" in readme_en

    assert "## Zielgruppen & Auffindbarkeit" in readme_de
    assert "Theoretische Physiker & Quantenfeldtheoretiker" in readme_de
    assert "Analytische Zahlentheoretiker & Millennium-Forscher" in readme_de
    assert "Open-Science Archivare & Reviewer" in readme_de
    assert "KI-Forschungsagenten & Wissens-Synthesizer" in readme_de


def test_third_party_licenses_readme_sections():
    """Validate that Third-Party Licenses sections are present in English and German READMEs."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    assert "## Third-Party Licenses & Transparency" in readme_en
    assert "THIRD_PARTY_LICENSES.md" in readme_en
    assert "Zero Copyleft" in readme_en

    assert "## Drittanbieter-Lizenzen & Transparenz" in readme_de
    assert "THIRD_PARTY_LICENSES.md" in readme_de
    assert "Zero Copyleft" in readme_de


def test_marketing_log_pfad_b_v105_entry():
    """Validate that MARKETING-LOG.txt documents the v1.0.5 Pfad B discoverability audit."""
    content = (REPO_ROOT / "MARKETING-LOG.txt").read_text(encoding="utf-8")
    assert "1.0.5" in content
    assert "2026-09-12" in content
    assert "THIRD_PARTY_LICENSES" in content
    assert "Wettbewerbsmatrix" in content
    assert "INV-LOCAL-01" in content
    assert "INV-SLA-10" in content


def test_ci_workflow_least_privilege_and_timeouts():
    """Validate that CI workflows enforce least-privilege permissions, concurrency, and timeouts."""
    ci_path = REPO_ROOT / ".github" / "workflows" / "ci.yml"
    ci_content = ci_path.read_text(encoding="utf-8")
    assert "permissions:" in ci_content
    assert "contents: read" in ci_content
    assert "timeout-minutes: 15" in ci_content

    stale_path = REPO_ROOT / ".github" / "workflows" / "stale.yml"
    stale_content = stale_path.read_text(encoding="utf-8")
    assert "actions/stale@v9" in stale_content
    assert "concurrency:" in stale_content
    assert "timeout-minutes: 10" in stale_content

    welcome_path = REPO_ROOT / ".github" / "workflows" / "welcome.yml"
    welcome_content = welcome_path.read_text(encoding="utf-8")
    assert "actions/first-interaction@v3" in welcome_content
    assert "concurrency:" in welcome_content
    assert "timeout-minutes: 5" in welcome_content


def test_multihost_cloud_sync_and_lock_hygiene():
    """Validate extended multi-host conflict copies, lock files, and cache in .gitignore."""
    gitignore = (REPO_ROOT / ".gitignore").read_text(encoding="utf-8")

    # Cloud sync conflict patterns
    assert "*conflicted copy*" in gitignore
    assert "* (Kopie)*" in gitignore
    assert "* (Copy)*" in gitignore
    assert "*-ASUS*" in gitignore
    assert "*-ASUS-GEI*" in gitignore
    assert "*-LAPTOP*" in gitignore
    assert "*-WORKSTATION*" in gitignore
    assert "*-WORKSTATION-LG*" in gitignore
    assert "*-Mac Studio*" in gitignore
    assert "*-MacBook*" in gitignore
    assert "*-IDEAPAD*" in gitignore

    # Canonical lock patterns
    assert "LOCK" in gitignore
    assert "LOCK.*" in gitignore
    assert "LOCK*.txt" in gitignore
    assert "LOCK.permissions.json" in gitignore
    assert "LOCK.user.*" in gitignore
    assert "LOCK.until.*" in gitignore
    assert "LOCK.condition.*" in gitignore
    assert ".automation-lock" in gitignore
    assert "uv.lock" in gitignore
    assert "!package-lock.json" in gitignore

    # Extended test/cache patterns
    assert ".pytest_temp/" in gitignore
    assert ".hypothesis/" in gitignore
    assert ".turbo/" in gitignore
    assert ".nyc_output/" in gitignore
    assert "*.orig" in gitignore
    assert "*.rej" in gitignore
    assert "Desktop.ini" in gitignore


def test_pytest_guardrails_and_license_files():
    """Validate pytest guardrails and license-files in pyproject.toml."""
    pyproject_path = REPO_ROOT / "pyproject.toml"
    with open(pyproject_path, "rb") as f:
        data = tomllib.load(f)

    project = data.get("project", {})
    assert "license-files" in project
    assert "LICENSE" in project["license-files"]
    assert "NOTICE" in project["license-files"]
    assert "THIRD_PARTY_LICENSES.md" in project["license-files"]

    pytest_opts = data.get("tool", {}).get("pytest", {}).get("ini_options", {})
    assert pytest_opts.get("minversion") == "7.0"
    assert "--basetemp=.pytest_temp" in pytest_opts.get("addopts", "")
    norecursedirs = pytest_opts.get("norecursedirs", [])
    assert ".git" in norecursedirs
    assert ".pytest_cache" in norecursedirs
    assert ".pytest_temp" in norecursedirs
    assert ".hypothesis" in norecursedirs
    assert "applications" in norecursedirs
    assert "scripts" in norecursedirs


def test_marketing_log_pfad_a_v106_entry():
    """Validate that MARKETING-LOG.txt documents the v1.0.6 Pfad A technical hygiene audit."""
    content = (REPO_ROOT / "MARKETING-LOG.txt").read_text(encoding="utf-8")
    assert "1.0.6" in content
    assert "2026-09-18" in content
    assert "GITHUBBOT_ONE_REPO_CLEANER" in content
    assert "Pfad A" in content
    assert "Least-Privilege" in content
    assert "Pytest-Guardrails" in content or "Pytest Guardrails" in content


def test_third_party_licenses_audit_v106():
    """Validate that THIRD_PARTY_LICENSES.md has been audited to v1.0.6."""
    licenses_file = REPO_ROOT / "THIRD_PARTY_LICENSES.md"
    content = licenses_file.read_text(encoding="utf-8")
    assert "- **Version:** `1.0.6`" in content
    assert ("- **Audit Date:** `2026-09-26`" in content or "- **Audit Date:** `2026-09-20`" in content)
    assert "100% Permissive Open Source" in content


def test_security_supported_versions_parity():
    """Validate that SECURITY.md supported versions table documents 1.0.x active branch."""
    security_file = REPO_ROOT / "SECURITY.md"
    content = security_file.read_text(encoding="utf-8")
    assert "1.0.x" in content
    assert "< 1.0.0" in content
    assert "RunAsInvoker" in content


def test_comparative_matrix_parity():
    """Validate that Section 16 Comparative Matrix exists with 10 dimensions in EN and DE."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    assert "## Comparative Matrix vs. Alternatives" in readme_en
    assert "## Vergleichsmatrix vs. Alternative Methoden" in readme_de

    # Check for all 10 numbered invariant dimensions
    for i in range(1, 11):
        num_str = f"**{i:02d}."
        assert num_str in readme_en, f"Dimension {num_str} missing in README.md comparative matrix"
        assert num_str in readme_de, f"Dimension {num_str} missing in README_de.md comparative matrix"

    # Check key columns
    assert "FST (research-line)" in readme_en
    assert "Connes Noncommutative Geometry" in readme_en
    assert "FST (research-line)" in readme_de
    assert "Connes Nichtkommutative Geometrie" in readme_de


def test_statutory_notice_bgb_521_parity():
    """Validate that Section 18 documents author, license, and § 521 BGB statutory notice."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    assert "## Security Policy, Author & Statutory Notice" in readme_en
    assert "## Sicherheitsrichtlinie, Autor & Gesetzlicher Haftungsausschluss" in readme_de

    assert "521 BGB" in readme_en
    assert "521 BGB" in readme_de
    assert "Gefälligkeitsrecht" in readme_de or "Gefaelligkeitsrecht" in readme_de
    assert "0009-0005-7296-1534" in readme_en
    assert "0009-0005-7296-1534" in readme_de


def test_reciprocal_anchors_parity():
    """Validate that reciprocal HTML anchor pairs exist for all 18 sections across READMEs."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    anchors = [
        "start-here",
        "einstieg",
        "discovery-context",
        "auffindbarkeit--entdeckungskontext",
        "the-five-masters",
        "die-fuenf-master-arbeiten",
        "domain-supplements",
        "domaenen-ergaenzungen--anwendungen",
        "glossary--fst-core-terms",
        "glossar--fst-kernbegriffe",
        "proof-architecture",
        "beweisarchitektur",
        "theoretical-data-flow--validation-sequence",
        "theoretischer-datenfluss--validierungssequenz",
        "numerical-validation-lifecycle",
        "numerischer-validierungs-lebenszyklus",
        "governance--runtime-invariants",
        "governance---laufzeit-invarianten",
        "proof-architecture-ascii-overview",
        "ascii-uebersicht-der-beweisarchitektur",
        "independent-foundations",
        "unabhaengige-fundamente",
        "numerical-validation-scripts",
        "numerische-validierungsskripte",
        "repository-structure",
        "repository-struktur",
        "ecosystem--sibling-research-repositories",
        "oekosystem--verwandte-forschungs-repositories",
        "target-personas--discoverability",
        "zielgruppen--auffindbarkeit",
        "comparative-matrix--alternatives",
        "vergleichsmatrix--alternative-methoden",
        "third-party-licenses--transparency",
        "drittanbieter-lizenzen--transparenz",
        "security-policy--statutory-notice",
        "sicherheitsrichtlinie--gesetzlicher-haftungsausschluss",
    ]

    for anchor in anchors:
        tag = f'id="{anchor}"'
        assert tag in readme_en, f"Anchor {anchor} missing in README.md"
        assert tag in readme_de, f"Anchor {anchor} missing in README_de.md"


def test_marketing_log_pfad_b_v106_entry():
    """Validate that MARKETING-LOG.txt documents the 2026-09-20 Pfad B audit."""
    content = (REPO_ROOT / "MARKETING-LOG.txt").read_text(encoding="utf-8")
    assert "2026-09-20" in content
    assert "GITHUBBOT_ONE_REPO_MARKETING_AND_DESIGN" in content
    assert "Pfad B" in content
    assert "18-Punkte" in content
    assert "521 BGB" in content
    assert "Vergleichsmatrix" in content


def test_canonical_notice_attribution():
    """Validate that canonical NOTICE attribution file exists and declares proper provenance."""
    notice_path = REPO_ROOT / "NOTICE"
    assert notice_path.is_file(), "NOTICE must exist in repository root"
    content = notice_path.read_text(encoding="utf-8")
    assert "Functional Stability Theory" in content
    assert "Lukas Geiger" in content
    assert "research-line" in content
    assert "open-bricks" in content
    assert "Creative Commons Attribution 4.0" in content or "CC BY 4.0" in content


def test_pep621_saturated_keywords():
    """Validate that pyproject.toml defines saturated 20/20 keywords matching GitHub topics."""
    pyproject_path = REPO_ROOT / "pyproject.toml"
    with open(pyproject_path, "rb") as f:
        data = tomllib.load(f)

    keywords = data.get("project", {}).get("keywords", [])
    assert len(keywords) == 20
    assert "functional-stability-theory" in keywords
    assert "rfep" in keywords
    assert "riemann-hypothesis" in keywords
    assert "cosmology" in keywords
    assert "mathematical-physics" in keywords
    assert "zenodo-doi" in keywords


def test_third_party_licenses_audit_recency_20260926():
    """Validate that THIRD_PARTY_LICENSES.md audit date is 2026-09-26 and references NOTICE."""
    licenses_file = REPO_ROOT / "THIRD_PARTY_LICENSES.md"
    content = licenses_file.read_text(encoding="utf-8")
    assert "- **Audit Date:** `2026-09-26`" in content
    assert ("[`NOTICE`](NOTICE)" in content or "[NOTICE](NOTICE)" in content)
    assert "100% Permissive Open Source" in content


def test_marketing_log_pfad_a_cleaner_20260926():
    """Validate that MARKETING-LOG.txt contains the 2026-09-26 GITHUBBOT_ONE_REPO_CLEANER Pfad A entry."""
    content = (REPO_ROOT / "MARKETING-LOG.txt").read_text(encoding="utf-8")
    assert "2026-09-26" in content
    assert "GITHUBBOT_ONE_REPO_CLEANER" in content
    assert "(Pfad A)" in content
    assert "NOTICE" in content


def test_changelog_unreleased_pfad_a():
    """Validate that CHANGELOG.md documents Pfad A improvements under [Unreleased]."""
    content = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert "## [Unreleased]" in content
    assert "NOTICE" in content
    assert "PEP 621" in content
    assert "license-files" in content
