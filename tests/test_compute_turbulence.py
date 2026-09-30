"""Automated contract and reproducibility tests for Turbulence domain."""

import importlib.util
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent

TURB_DIR = REPO_ROOT / "fst-physics" / "turbulence"
SCRIPTS_DIR = REPO_ROOT / "scripts" / "turbulence"

# Dynamically import compute_dual_dfc
DUAL_DFC_PATH = SCRIPTS_DIR / "compute_dual_dfc.py"
SPEC_DUAL = importlib.util.spec_from_file_location("compute_dual_dfc", DUAL_DFC_PATH)
MODULE_DUAL = importlib.util.module_from_spec(SPEC_DUAL)
sys.modules[SPEC_DUAL.name] = MODULE_DUAL
SPEC_DUAL.loader.exec_module(MODULE_DUAL)


def test_turbulence_paper_artifacts_presence():
    """Verify that paper sources and compiled PDFs exist and are non-empty."""
    en_tex = TURB_DIR / "FST-TU_Turbulence_Skeleton_v1_en.tex"
    de_tex = TURB_DIR / "FST-TU_Turbulence_Skeleton_v1_de.tex"
    en_pdf = TURB_DIR / "FST-TU_Turbulence_Skeleton_v1_en.pdf"
    de_pdf = TURB_DIR / "FST-TU_Turbulence_Skeleton_v1_de.pdf"
    kombi_pdf = TURB_DIR / "FST-TU_Turbulence_Skeleton_v1_kombi.pdf"

    for path in [en_tex, de_tex, en_pdf, de_pdf, kombi_pdf]:
        assert path.exists(), f"Missing paper artifact: {path.name}"
        assert path.stat().st_size > 1000, f"Artifact too small: {path.name}"


def test_compute_dual_dfc_scenarios():
    """Verify dual DFC1 evaluation across canonical controls."""
    # Positive control: smooth forward
    res_smooth = MODULE_DUAL.evaluate(
        "smooth_forward", "positive control", "smooth", "smooth_forward"
    )
    assert res_smooth.dual_dfc1_status == "pass"
    assert res_smooth.pointwise_dfc1_status == "pass"
    assert res_smooth.verdict == "pass"
    assert res_smooth.weighted_flux_ratio > 1.0

    # Negative control: bad corridor
    res_bad = MODULE_DUAL.evaluate(
        "bad_corridor", "backscatter in high weights", "front_loaded", "bad_corridor"
    )
    assert res_bad.dual_dfc1_status == "fail"
    assert "fail" in res_bad.verdict.lower()
    assert res_bad.weighted_flux_ratio < 1.0

    # Vacuous reference control: k41
    res_k41 = MODULE_DUAL.evaluate(
        "k41_reference", "reference profile", "k41", "k41"
    )
    assert res_k41.dual_dfc1_status == "pass"
    assert res_k41.weighted_flux_ratio == 1.0


def test_turbulence_phi_profile_bounds():
    """Verify that phi profiles maintain non-negative values and proper lengths."""
    for kind in ["k41", "smooth", "front_loaded", "slush"]:
        prof = MODULE_DUAL.phi_profile(kind, 10)
        assert len(prof) == 10
        assert all(val >= 0.0 for val in prof)


def test_turbulence_flux_profile_properties():
    """Verify that flux profiles generate expected length and forward properties."""
    for kind in ["k41", "smooth_forward", "alternating_tolerated", "bad_corridor", "same_mean_shuffle"]:
        fl = MODULE_DUAL.flux_profile(kind, 10)
        assert len(fl) == 9
        if kind == "smooth_forward":
            assert all(v > 0 for v in fl)
