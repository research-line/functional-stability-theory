"""Automated contract and reproducibility tests for BSD Positivity domain."""

import importlib.util
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent

VERIFY_SCRIPT_PATH = REPO_ROOT / "scripts" / "bsd" / "compute_bsd_verification.py"
SPEC_VERIFY = importlib.util.spec_from_file_location("compute_bsd_verification", VERIFY_SCRIPT_PATH)
MODULE_VERIFY = importlib.util.module_from_spec(SPEC_VERIFY)
sys.modules[SPEC_VERIFY.name] = MODULE_VERIFY
SPEC_VERIFY.loader.exec_module(MODULE_VERIFY)

RANK2_SCRIPT_PATH = REPO_ROOT / "scripts" / "bsd" / "compute_rank2_lmfdb.py"
SPEC_RANK2 = importlib.util.spec_from_file_location("compute_rank2_lmfdb", RANK2_SCRIPT_PATH)
MODULE_RANK2 = importlib.util.module_from_spec(SPEC_RANK2)
sys.modules[SPEC_RANK2.name] = MODULE_RANK2
SPEC_RANK2.loader.exec_module(MODULE_RANK2)


def test_bsd_curves_database_integrity():
    """Verify that BSD reference curve dataset has required keys and valid data."""
    curves = MODULE_VERIFY.BSD_CURVES
    assert len(curves) >= 4
    for label in ["11a1", "37a1", "389a1", "5077a1"]:
        assert label in curves
        data = curves[label]
        assert data["omega"] > 0
        assert data["regulator"] > 0
        assert data["torsion_order"] >= 1
        assert data["sha"] >= 1
        assert data["L_value"] > 0


def test_bsd_formula_verification():
    """Verify that BSD formula numerical ratio is within 0.1% for all curves."""
    results = MODULE_VERIFY.verify_bsd_formula()
    assert len(results) == 4
    for label, rank, match, ratio in results:
        assert match, f"Curve {label} (rank {rank}) did not match BSD formula: ratio={ratio}"
        assert abs(ratio - 1.0) < 0.001


def test_bsd_regulator_positivity():
    """Verify that regulators and height matrices are strictly positive (Axiom II)."""
    curves = MODULE_VERIFY.BSD_CURVES
    for label, data in curves.items():
        assert data["regulator"] > 0
        if data["rank"] > 0:
            h_mat = data["height_matrix"]
            eigvals = MODULE_VERIFY.np.linalg.eigvalsh(h_mat)
            assert (eigvals > 0).all()


def test_bsd_rank2_cremona_regulators_strictly_positive():
    """Verify that all rank-2 sample curves have strictly positive regulators."""
    curves = MODULE_RANK2.RANK2_CURVES
    assert len(curves) == 16
    for _label, n, r, reg, _sha_an, omega, c_prod in curves:
        assert reg > 0
        assert n > 0
        assert r == 2
        assert omega > 0
        assert c_prod >= 1
        l_star = 2 * omega * reg * c_prod
        assert l_star > 0


def test_bsd_supplement_artifacts_exist():
    """Verify that paper artifacts and reproducibility scripts exist."""
    bsd_dir = REPO_ROOT / "fst-mathematics" / "bsd"
    assert (bsd_dir / "BSD_Positivity_EN.tex").is_file()
    assert (bsd_dir / "BSD_Positivity_DE.tex").is_file()
    assert (bsd_dir / "BSD_Positivity_EN.pdf").is_file()
    assert (bsd_dir / "BSD_Positivity_DE.pdf").is_file()
    assert (bsd_dir / "BSD_Positivity_kombi.pdf").is_file()
    assert (bsd_dir / "README.md").is_file()
    assert (REPO_ROOT / "scripts" / "bsd" / "compute_bsd_verification.py").is_file()
    assert (REPO_ROOT / "scripts" / "bsd" / "compute_height_saturation.py").is_file()
    assert (REPO_ROOT / "scripts" / "bsd" / "compute_rank2_lmfdb.py").is_file()
