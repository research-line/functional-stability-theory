"""Automated contract and reproducibility tests for Navier-Stokes domain."""

import importlib.util
from pathlib import Path
import sys
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent

NS_DIR = REPO_ROOT / "fst-physics" / "navier-stokes"
SCRIPTS_DIR = REPO_ROOT / "scripts" / "navier-stokes"

# Dynamically import compute_ds3_lorenz
DS3_PATH = SCRIPTS_DIR / "compute_ds3_lorenz.py"
SPEC_DS3 = importlib.util.spec_from_file_location("compute_ds3_lorenz", DS3_PATH)
MODULE_DS3 = importlib.util.module_from_spec(SPEC_DS3)
sys.modules[SPEC_DS3.name] = MODULE_DS3
SPEC_DS3.loader.exec_module(MODULE_DS3)


def test_navier_stokes_paper_artifacts_presence():
    """Verify that paper sources and compiled PDFs exist and are non-empty."""
    expected_artifacts = [
        "FST-NS_NavierStokes_Skeleton_v1_en.tex",
        "FST-NS_NavierStokes_Skeleton_v1_de.tex",
        "FST-NS_NavierStokes_Skeleton_v1_en.pdf",
        "FST-NS_NavierStokes_Skeleton_v1_de.pdf",
        "FST-NS_NavierStokes_Skeleton_v1_kombi.pdf",
        "FST-NS_LogDistanceIntegrability_v1_en.tex",
        "FST-NS_LogDistanceIntegrability_v1_de.tex",
        "FST-NS_LogDistanceIntegrability_v1_en.pdf",
        "FST-NS_LogDistanceIntegrability_v1_de.pdf",
        "FST-NS_LogDistanceIntegrability_v1_kombi.pdf",
    ]

    for name in expected_artifacts:
        path = NS_DIR / name
        assert path.exists(), f"Missing paper artifact: {name}"
        assert path.stat().st_size > 1000, f"Artifact too small: {name}"


def test_navier_stokes_lorenz_vector_field():
    """Verify Lorenz vector field fixed points and velocity evaluations."""
    rho, beta = MODULE_DS3.RHO, MODULE_DS3.BETA

    # Origin is a fixed point
    f_origin = MODULE_DS3.lorenz(0.0, [0.0, 0.0, 0.0])
    assert np.allclose(f_origin, [0.0, 0.0, 0.0])

    # Fixed points C+ and C-
    x_c = np.sqrt(beta * (rho - 1.0))
    z_c = rho - 1.0
    f_cp = MODULE_DS3.lorenz(0.0, [x_c, x_c, z_c])
    assert np.allclose(f_cp, [0.0, 0.0, 0.0], atol=1e-12)

    f_cm = MODULE_DS3.lorenz(0.0, [-x_c, -x_c, z_c])
    assert np.allclose(f_cm, [0.0, 0.0, 0.0], atol=1e-12)


def test_navier_stokes_total_variation_operator():
    """Verify total variation calculation on discrete test trajectories."""
    # Stationary trajectory has TV = 0
    stat_traj = np.ones((10, 3))
    assert np.isclose(MODULE_DS3.total_variation(stat_traj), 0.0)

    # Unit step in x direction at each step
    unit_step_traj = np.array([[i, 0.0, 0.0] for i in range(5)])
    assert np.isclose(MODULE_DS3.total_variation(unit_step_traj), 4.0)
