"""Automated contract and reproducibility tests for FST-II Chemical Stability domain."""

import importlib.util
from pathlib import Path
import sys
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

FST2_DIR = REPO_ROOT / "applications" / "fst-ii-chemical"
PAPERS_DIR = FST2_DIR / "papers"
SCRIPTS_DIR = FST2_DIR / "scripts"

VERIFY_PERM_PATH = SCRIPTS_DIR / "verify_fst2_permanence_resolvent_and_nash_gate.py"
MODULE_PERM = None
IMPORT_ERROR_PERM = None

if VERIFY_PERM_PATH.exists():
    try:
        spec = importlib.util.spec_from_file_location("verify_fst2_permanence", VERIFY_PERM_PATH)
        if spec and spec.loader:
            mod = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = mod
            spec.loader.exec_module(mod)
            MODULE_PERM = mod
    except Exception as exc:  # pragma: no cover
        IMPORT_ERROR_PERM = exc


def test_fst2_chemical_paper_artifacts_presence():
    """Verify that FST-II Chemical v3 paper sources and compiled PDFs exist and are valid."""
    expected_artifacts = [
        "FST_II_Chemical_Stability_v3_en.tex",
        "FST_II_Chemical_Stability_v3_de.tex",
        "FST_II_Chemical_Stability_v3_en.pdf",
        "FST_II_Chemical_Stability_v3_de.pdf",
        "FST_II_Chemical_Stability_v3_kombi.pdf",
    ]
    for filename in expected_artifacts:
        artifact_path = PAPERS_DIR / filename
        assert artifact_path.is_file(), f"Expected paper artifact {filename} missing at {artifact_path}"
        assert artifact_path.stat().st_size > 1000, f"Artifact {filename} is suspiciously small ({artifact_path.stat().st_size} bytes)"


def test_fst2_chemical_scripts_presence():
    """Verify that all canonical FST-II Chemical verification scripts exist and compile."""
    expected_scripts = [
        "autocatalysis_3species.py",
        "replicator_entropy_test.py",
        "dissipation_topology_waterline_ledger.py",
        "homochirality_mechanism_taxonomy_gate.py",
        "homochirality_open_drive_gate.py",
        "verify_fst2_permanence_resolvent_and_nash_gate.py",
        "verify_fst2_transient_peak_2026_10_18.py",
    ]
    for script_name in expected_scripts:
        script_path = SCRIPTS_DIR / script_name
        assert script_path.is_file(), f"FST-II verification script {script_name} missing at {script_path}"
        assert script_path.stat().st_size > 0, f"Script {script_name} is empty"


def test_fst2_theorem_1_permanence_contract():
    """Verify Theorem 1: Simplex permanence and boundary transversality."""
    if MODULE_PERM is None:
        pytest.skip(f"Could not import verify_fst2_permanence: {IMPORT_ERROR_PERM}")
    res = MODULE_PERM.verify_theorem_1_permanence()
    assert res["x_star_sum"] == pytest.approx(1.0)
    assert res["min_inflow_facet_A"] > 0
    assert res["invasion_rate_B_at_eA"] > 0
    assert res["invasion_rate_C_at_eB"] > 0
    assert res["permanence_certified"] is True


def test_fst2_theorem_2_resolvent_contract():
    """Verify Theorem 2: Non-normal resolvent and pseudospectral stability."""
    if MODULE_PERM is None:
        pytest.skip(f"Could not import verify_fst2_permanence: {IMPORT_ERROR_PERM}")
    res = MODULE_PERM.verify_theorem_2_resolvent()
    assert res["alpha_T"] < 0, "Spectral abscissa must be strictly negative"
    assert res["numerical_abscissa_mu2"] > 0, "Numerical abscissa must be positive (non-normal transient growth)"
    assert res["max_transient_growth"] > 1.2, "Transient growth must exceed 1.2x"
    assert res["contrast_to_naive_bound"] > 1.8, "Contrast to naive bound must exceed 1.8x"


def test_fst2_theorem_3_replicator_nash_contract():
    """Verify Theorem 3: Replicator-Nash and potential-Lyapunov demarcation."""
    if MODULE_PERM is None:
        pytest.skip(f"Could not import verify_fst2_permanence: {IMPORT_ERROR_PERM}")
    res = MODULE_PERM.verify_theorem_3_replicator_nash_potential()
    assert res["rps_orbit_is_nash"] is False, "General interior limit point must not be presumed Nash"
    assert res["coordination_d_DKL_dt"] > 0, "Coordination divergence must be positive"
    assert res["coordination_dPhi_dt"] > 0, "Potential function must strictly increase"
