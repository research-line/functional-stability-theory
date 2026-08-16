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
    assert "Ecosystem-research--line" in readme_en
    assert "Umbrella-open--bricks" in readme_en
    assert "llms.txt" in readme_en
    assert "Tests-8%2F8%20Passed" in readme_en
    assert "Python-3.10%2B" in readme_en

    # Required badge signatures in German README
    assert "Lizenz-CC_BY_4.0" in readme_de or "License-CC_BY_4.0" in readme_de
    assert "ORCID-0009--0005--7296--1534" in readme_de
    assert "Zenodo-10.5281" in readme_de
    assert "research--line" in readme_de
    assert "open--bricks" in readme_de
    assert "llms.txt" in readme_de
    assert "Tests-8%2F8%20Bestanden" in readme_de or "Tests-8%2F8%20Passed" in readme_de
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

    assert "## Last-checked: 2026-08-16" in llms_txt
    assert "https://github.com/research-line/functional-stability-theory" in llms_txt
    assert "research-line" in llms_txt
    assert "Renormalized Free-Energy Principle" in llms_txt
    assert "Pattern A" in llms_txt
    assert "## Search Phrases" in llms_txt


def test_changelog_entry():
    """Validate that CHANGELOG.md documents the 2026-08-16 discoverability overhaul."""
    changelog = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert "2026-08-16" in changelog
    assert "Discoverability" in changelog or "discoverability" in changelog or "metadata" in changelog


def test_mermaid_diagrams_parity():
    """Validate that both English and German READMEs contain Mermaid diagrams."""
    readme_en = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    readme_de = (REPO_ROOT / "README_de.md").read_text(encoding="utf-8")

    assert "```mermaid" in readme_en
    assert "flowchart TD" in readme_en
    assert "subgraph MASTERS" in readme_en

    assert "```mermaid" in readme_de
    assert "flowchart TD" in readme_de
    assert "subgraph MASTERS" in readme_de
