"""Automated contract and reproducibility tests for K41 Variational Minimiser domain."""

import importlib.util
from pathlib import Path
import sys
from unittest.mock import MagicMock
import numpy as np

try:
    import fitz  # PyMuPDF
except ImportError:
    fitz = None

REPO_ROOT = Path(__file__).resolve().parent.parent

K41_DIR = REPO_ROOT / "fst-physics" / "k41-variational-minimiser"
SCRIPTS_DIR = REPO_ROOT / "scripts" / "k41"

# Dynamically import compute_F_spectrum with headless matplotlib mock isolation
_mocked_matplotlib = False
if "matplotlib" not in sys.modules:
    mock_plt = MagicMock()
    mock_plt.subplots.return_value = (MagicMock(), (MagicMock(), MagicMock()))
    mock_mpl = MagicMock()
    mock_mpl.pyplot = mock_plt
    sys.modules["matplotlib"] = mock_mpl
    sys.modules["matplotlib.pyplot"] = mock_plt
    _mocked_matplotlib = True

F_SPEC_PATH = SCRIPTS_DIR / "compute_F_spectrum.py"
SPEC_F = importlib.util.spec_from_file_location("compute_F_spectrum", F_SPEC_PATH)
MODULE_F = importlib.util.module_from_spec(SPEC_F)
sys.modules[SPEC_F.name] = MODULE_F
SPEC_F.loader.exec_module(MODULE_F)

if _mocked_matplotlib:
    sys.modules.pop("matplotlib.pyplot", None)
    sys.modules.pop("matplotlib", None)


def test_k41_paper_artifacts_presence():
    """Verify that K41 paper sources and compiled PDFs exist and are valid."""
    en_tex = K41_DIR / "K41_Variational_Minimiser_v1_en.tex"
    de_tex = K41_DIR / "K41_Variational_Minimiser_v1_de.tex"
    en_pdf = K41_DIR / "K41_Variational_Minimiser_v1_en.pdf"
    de_pdf = K41_DIR / "K41_Variational_Minimiser_v1_de.pdf"
    kombi_pdf = K41_DIR / "K41_Variational_Minimiser_v1_kombi.pdf"

    for path in [en_tex, de_tex, en_pdf, de_pdf, kombi_pdf]:
        assert path.exists(), f"Missing paper artifact: {path.name}"
        assert path.stat().st_size > 10000, f"Artifact unexpectedly small: {path.name}"

    if fitz is not None:
        doc_en = fitz.open(en_pdf)
        assert len(doc_en) >= 20, f"EN PDF should have at least 20 pages, got {len(doc_en)}"
        txt_en = "".join(page.get_text() for page in doc_en)
        assert txt_en.count("\ufffd") == 0, "EN PDF textlayer contains replacement characters"

        doc_de = fitz.open(de_pdf)
        assert len(doc_de) >= 20, f"DE PDF should have at least 20 pages, got {len(doc_de)}"
        txt_de = "".join(page.get_text() for page in doc_de)
        assert txt_de.count("\ufffd") == 0, "DE PDF textlayer contains replacement characters"

        doc_kombi = fitz.open(kombi_pdf)
        assert len(doc_kombi) >= 40, f"Kombi PDF should have at least 40 pages, got {len(doc_kombi)}"
        txt_kombi = "".join(page.get_text() for page in doc_kombi)
        assert txt_kombi.count("\ufffd") == 0, "Kombi PDF textlayer contains replacement characters"


def test_k41_free_energy_minimizer_uniqueness():
    """Verify that F[E] achieves unique minimum at reference spectrum E_star."""
    k_centers = MODULE_F.k_centers
    delta_k = MODULE_F.delta_k
    E_ref = MODULE_F.E_star(k_centers, delta_k)

    # Free energy at minimum must vanish identically: F[E*] = 0
    F_min = MODULE_F.compute_F(E_ref, E_ref)
    assert np.isclose(F_min, 0.0, atol=1e-12), f"F[E*] must be 0, got {F_min}"

    # Strict positivity for any non-trivial scaled perturbation c != 1
    for c in [0.25, 0.5, 0.8, 1.25, 1.75, 2.5]:
        E_pert = c * E_ref
        F_pert = MODULE_F.compute_F(E_pert, E_ref)
        assert F_pert > 0.0, f"F[c*E*] must be strictly positive for c={c}, got {F_pert}"

    # Strict positivity for power-law tilt perturbations
    for alpha in [-0.5, -0.2, 0.2, 0.5]:
        tilt = (k_centers / k_centers[0]) ** alpha
        E_tilt = E_ref * tilt
        F_tilt = MODULE_F.compute_F(E_tilt, E_ref)
        assert F_tilt > 0.0, f"F[tilt] must be strictly positive for alpha={alpha}"


def test_k41_kolmogorov_scaling():
    """Verify that E_star follows the exact Kolmogorov k^(-5/3) power law."""
    k_centers = MODULE_F.k_centers
    delta_k = MODULE_F.delta_k
    E_ref = MODULE_F.E_star(k_centers, delta_k)

    # Spectral density e(k) = E_ref / delta_k
    e_density = E_ref / delta_k

    # Log-log slope must match -5/3 across the inertial range
    log_k = np.log(k_centers)
    log_e = np.log(e_density)
    slope, _ = np.polyfit(log_k, log_e, 1)

    expected_slope = -5.0 / 3.0
    assert np.isclose(slope, expected_slope, atol=1e-8), (
        f"Kolmogorov slope must be {expected_slope}, got {slope}"
    )


def test_k41_dfc_spectral_steepness():
    """Verify DFC1 and DFC2 spectral steepness criteria for Kolmogorov scaling."""
    k_centers = MODULE_F.k_centers
    delta_k = MODULE_F.delta_k
    E_ref = MODULE_F.E_star(k_centers, delta_k)

    dfc1_ok = MODULE_F.check_DFC1(E_ref, k_centers, delta_k)
    assert dfc1_ok, "K41 reference spectrum must satisfy DFC1 forward cascade"

    dfc2_weak, dfc2_strong, mean_slope, max_slope, phi_nonincreasing = MODULE_F.check_DFC2(
        E_ref, k_centers, delta_k
    )
    assert dfc2_weak, "K41 must satisfy weak DFC2 threshold"
    assert dfc2_strong, "K41 must satisfy strong K41 DFC2 threshold"
    assert phi_nonincreasing, "K41 must satisfy non-increasing phi criterion"
