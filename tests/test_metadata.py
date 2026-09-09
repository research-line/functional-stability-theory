"""Metadata and discoverability parity tests for Functional Stability Theory (FST)."""

import pathlib
import tomllib


REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent


def test_pyproject_metadata():
    """Validate pyproject.toml PEP 621 metadata."""
    pyproject_path = REPO_ROOT / "pyproject.toml"
    assert pyproject_path.is_file(), "pyproject.toml must exist"

    with open(pyproject_path, "rb") as f:
        data = tomllib.load(f)

    project = data.get("project", {})
    assert project.get("name") == "functional-stability-theory"
    assert project.get("version") == "1.0.3"
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
    assert "version-1.0.3" in readme_en
    assert "actions/workflows/ci.yml" in readme_en
    assert "ORCID-0009--0005--7296--1534" in readme_en
    assert "Zenodo-10.5281" in readme_en
    assert "Security-Research%20Integrity" in readme_en
    assert "Security-RunAsInvoker" in readme_en
    assert "Security%20SLA" in readme_en
    assert "Ecosystem-research--line" in readme_en
    assert "Umbrella-open--bricks" in readme_en
    assert "llms.txt" in readme_en
    assert "Tests-110%2B%20Passed" in readme_en
    assert "Python-3.10--3.13" in readme_en or "Python-3.10" in readme_en
    assert "Platform-Windows" in readme_en
    assert "Audit-Last--checked%202026--09--09" in readme_en

    # Required badge signatures in German README
    assert "Lizenz-CC_BY_4.0" in readme_de or "License-CC_BY_4.0" in readme_de
    assert "Version-1.0.3" in readme_de
    assert "actions/workflows/ci.yml" in readme_de
    assert "ORCID-0009--0005--7296--1534" in readme_de
    assert "Zenodo-10.5281" in readme_de
    assert "Sicherheit-Open%20Science%20Integrit%C3%A4t" in readme_de
    assert "Sicherheit-RunAsInvoker" in readme_de
    assert "Sicherheits--SLA" in readme_de
    assert "research--line" in readme_de
    assert "open--bricks" in readme_de
    assert "llms.txt" in readme_de
    assert "Tests-110%2B%20Bestanden" in readme_de
    assert "Python-3.10--3.13" in readme_de or "Python-3.10" in readme_de
    assert "Plattform-Windows" in readme_de
    assert "Audit-Gepr%C3%BCft%202026--09--09" in readme_de


def test_quick_navigation_parity():
    """Validate 14-point quick navigation index parity in English and German READMEs."""
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
        "10. Numerische Validierungsskripte",
        "11. Ökosystem & Repositories",
        "12. Autor & Lizenz",
    ]
    for point in de_points:
        assert point in readme_de, f"Schnellnavigation point {point} missing in README_de.md"


def test_governance_invariants_table_parity():
    """Validate 10-point governance and runtime invariants table in both READMEs."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    assert "## Governance & Runtime Invariants" in readme_en
    assert "### Governance- & Laufzeit-Invarianten" in readme_de

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

    assert "## Last-checked: 2026-09-09" in llms_txt
    assert "https://github.com/research-line/functional-stability-theory" in llms_txt
    assert "research-line" in llms_txt
    assert "Renormalized Free-Energy Principle" in llms_txt
    assert "Pattern A" in llms_txt
    assert "SECURITY.md" in llms_txt
    assert "## Search Phrases" in llms_txt
    assert "110+ Passed" in llms_txt
    assert "Governance & Runtime Invariants" in llms_txt
    assert "MARKETING-LOG.txt" in llms_txt


def test_changelog_entry():
    """Validate that CHANGELOG.md documents the latest release and hygiene audits."""
    changelog = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert "2026-09-09" in changelog
    assert "1.0.3" in changelog
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
