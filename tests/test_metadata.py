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
    assert project.get("version") == "1.0.0"
    assert project.get("requires-python") == ">=3.10"
    assert project.get("license", {}).get("text") == "CC-BY-4.0"
    assert len(project.get("authors", [])) >= 1
    assert project["authors"][0]["name"] == "Lukas Geiger"


def test_readme_and_readme_de_badges():
    """Validate badge parity in README.md and README_de.md."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    # Required badge signatures in English README
    assert "License-CC_BY_4.0" in readme_en
    assert "ORCID-0009--0005--7296--1534" in readme_en
    assert "Zenodo-10.5281" in readme_en
    assert "Security-Research%20Integrity" in readme_en
    assert "Ecosystem-research--line" in readme_en
    assert "Umbrella-open--bricks" in readme_en
    assert "llms.txt" in readme_en
    assert "Tests-10%2F10%20Passed" in readme_en
    assert "Python-3.10%2B" in readme_en

    # Required badge signatures in German README
    assert "Lizenz-CC_BY_4.0" in readme_de or "License-CC_BY_4.0" in readme_de
    assert "ORCID-0009--0005--7296--1534" in readme_de
    assert "Zenodo-10.5281" in readme_de
    assert "Sicherheit-Open%20Science%20Integrit%C3%A4t" in readme_de
    assert "research--line" in readme_de
    assert "open--bricks" in readme_de
    assert "llms.txt" in readme_de
    assert "Tests-10%2F10%20Bestanden" in readme_de
    assert "Python-3.10%2B" in readme_de


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

    assert "## Last-checked: 2026-08-21" in llms_txt
    assert "https://github.com/research-line/functional-stability-theory" in llms_txt
    assert "research-line" in llms_txt
    assert "Renormalized Free-Energy Principle" in llms_txt
    assert "Pattern A" in llms_txt
    assert "SECURITY.md" in llms_txt
    assert "## Search Phrases" in llms_txt


def test_changelog_entry():
    """Validate that CHANGELOG.md documents the 2026-08-21 discoverability audit."""
    changelog = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert "2026-08-21" in changelog
    assert "1.0.1" in changelog
    assert "SECURITY.md" in changelog


def test_mermaid_diagrams_parity():
    """Validate that both English and German READMEs contain all Mermaid diagrams."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    for doc, name in [(readme_en, "README.md"), (readme_de, "README_de.md")]:
        assert "```mermaid" in doc, f"{name} missing mermaid codeblock"
        assert "flowchart TD" in doc, f"{name} missing flowchart TD"
        assert "flowchart LR" in doc, f"{name} missing flowchart LR"
        assert "subgraph MASTERS" in doc, f"{name} missing subgraph MASTERS"
        assert "subgraph HYP" in doc, f"{name} missing subgraph HYP"


def test_security_policy_bilingual_parity():
    """Validate SECURITY.md exists and contains bilingual invariants and contact info."""
    security_file = REPO_ROOT / "SECURITY.md"
    assert security_file.is_file(), "SECURITY.md must exist"

    content = security_file.read_text(encoding="utf-8")
    assert "## English" in content
    assert "## Deutsche Fassung" in content
    assert "Zero-Egress" in content or "0% Netzwerk-Egress" in content
    assert "security@ellmos.ai" in content
    assert "support@lukasgeiger.com" in content


def test_sibling_matrix_parity():
    """Validate that sibling and ecosystem repositories are linked in both READMEs."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    expected_siblings = [
        "research-line/fst-nash",
        "research-line/rh-even-dominance",
        "research-line/crm-cosmology",
        "biotec-line/VFDistiller",
        "doc-bricks/MediaBrain",
        "dev-bricks/CodeBox",
        "dev-bricks/DevCenter",
        "ellmos-ai/skills",
        "ellmos-ai/sqlite-transit-sync",
        "open-bricks",
    ]

    for sibling in expected_siblings:
        assert sibling in readme_en, f"Sibling {sibling} missing in README.md"
        assert sibling in readme_de, f"Sibling {sibling} missing in README_de.md"

