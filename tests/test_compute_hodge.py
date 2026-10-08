"""Automated contract and reproducibility tests for Hodge Positivity domain."""

import importlib.util
from pathlib import Path
import sys
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

HODGE_DIR = REPO_ROOT / "fst-mathematics" / "hodge"
SCRIPTS_DIR = REPO_ROOT / "scripts" / "hodge"

VOISIN_PATH = SCRIPTS_DIR / "compute_voisin_test.py"
MODULE_VOISIN = None
IMPORT_ERROR = None

if VOISIN_PATH.exists():
    try:
        SPEC_VOISIN = importlib.util.spec_from_file_location("compute_voisin_test", VOISIN_PATH)
        if SPEC_VOISIN and SPEC_VOISIN.loader:
            mod = importlib.util.module_from_spec(SPEC_VOISIN)
            sys.modules[SPEC_VOISIN.name] = mod
            SPEC_VOISIN.loader.exec_module(mod)
            MODULE_VOISIN = mod
    except Exception as exc:  # pragma: no cover
        IMPORT_ERROR = exc


def test_hodge_paper_artifacts_presence():
    """Verify that paper sources and compiled PDFs exist and are non-empty."""
    en_tex = HODGE_DIR / "Hodge_Positivity_EN.tex"
    de_tex = HODGE_DIR / "Hodge_Positivity_DE.tex"
    en_pdf = HODGE_DIR / "Hodge_Positivity_EN.pdf"
    de_pdf = HODGE_DIR / "Hodge_Positivity_DE.pdf"
    kombi_pdf = HODGE_DIR / "Hodge_Positivity_kombi.pdf"

    for path in [en_tex, de_tex, en_pdf, de_pdf, kombi_pdf]:
        assert path.exists(), f"Missing paper artifact: {path.name}"
        assert path.stat().st_size > 1000, f"Artifact too small: {path.name}"


def test_hodge_riemann_form_computation():
    """Verify Hodge-Riemann form computation and primitive component extraction."""
    if MODULE_VOISIN is None:
        pytest.skip(f"compute_voisin_test could not be imported: {IMPORT_ERROR}")

    import numpy as np

    H = np.eye(4)
    Q = MODULE_VOISIN.hodge_riemann_form(H, None, 2, 4)
    assert np.allclose(Q, H)

    Q_prim, basis = MODULE_VOISIN.primitive_component(Q, None, 4)
    assert np.allclose(Q_prim, Q)
    assert basis.shape == (4, 4)


def test_hodge_test_cases_numerical_reproducibility():
    """Verify that Hodge obstruction test case functions evaluate consistently."""
    if MODULE_VOISIN is None:
        pytest.skip(f"compute_voisin_test could not be imported: {IMPORT_ERROR}")

    ap1_1, spec_1 = MODULE_VOISIN.test_generic_torus()
    assert bool(ap1_1) is True
    assert len(spec_1) == 800

    ap1_3, spec_3 = MODULE_VOISIN.test_k3_product()
    assert bool(ap1_3) is False
    assert len(spec_3) == 900

