#!/usr/bin/env python3
"""Verify the published FST-I 2D log-Hessian pilot.

This is intentionally a small, dependency-free rebuild wrapper.  The published
JSON is the frozen input: the wrapper does not pretend to rerun the stellar
model, but independently recomputes the symmetric 2x2 matrix eigenvalues and
checks the recorded finite-difference step convergence and pilot warning.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
APP_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT = APP_DIR / "results" / "stellar" / "log_hessian_alpha_me_2026-05-18.json"
EIGENVALUE_TOLERANCE = 1e-10
SYMMETRY_TOLERANCE = 1e-15
REQUIRED_STEPS = ("0.005", "0.01", "0.02", "0.05")


class ValidationError(ValueError):
    """Raised when the published pilot does not satisfy its declared schema."""


def _finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValidationError(f"{label} must be a number")
    number = float(value)
    if not math.isfinite(number):
        raise ValidationError(f"{label} must be finite")
    return number


def _vector(value: Any, label: str) -> list[float]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValidationError(f"{label} must contain exactly two values")
    return [_finite_number(item, f"{label}[{index}]") for index, item in enumerate(value)]


def _eigenvalues(matrix: list[list[float]]) -> list[float]:
    """Return ascending eigenvalues of a real symmetric 2x2 matrix."""

    a, b = matrix[0]
    _, d = matrix[1]
    centre = (a + d) / 2.0
    radius = math.hypot((a - d) / 2.0, b)
    return [centre - radius, centre + radius]


def _display_path(path: Path) -> str:
    repo_root = APP_DIR.parents[1]
    try:
        return path.resolve().relative_to(repo_root).as_posix()
    except ValueError:
        return path.name


def _validate(payload: dict[str, Any], input_path: Path) -> dict[str, Any]:
    required = {
        "created",
        "project",
        "source_script",
        "objective",
        "main_step",
        "H_log_main",
        "eigvalsh_main",
        "step_convergence",
        "interpretation",
    }
    missing = sorted(required.difference(payload))
    if missing:
        raise ValidationError(f"missing required fields: {', '.join(missing)}")

    matrix = payload["H_log_main"]
    if (
        not isinstance(matrix, list)
        or len(matrix) != 2
        or any(not isinstance(row, list) or len(row) != 2 for row in matrix)
    ):
        raise ValidationError("H_log_main must be a 2x2 matrix")
    matrix = [[_finite_number(item, f"H_log_main[{row}][{column}]") for column, item in enumerate(values)] for row, values in enumerate(matrix)]
    if abs(matrix[0][1] - matrix[1][0]) > SYMMETRY_TOLERANCE:
        raise ValidationError("H_log_main must be symmetric")

    reported = _vector(payload["eigvalsh_main"], "eigvalsh_main")
    recomputed = _eigenvalues(matrix)
    errors = [abs(left - right) for left, right in zip(recomputed, reported)]
    if max(errors) > EIGENVALUE_TOLERANCE:
        raise ValidationError(
            "reported eigenvalues disagree with H_log_main "
            f"(max error {max(errors):.3e})"
        )
    if reported != sorted(reported):
        raise ValidationError("eigvalsh_main must be ascending")

    main_step = _finite_number(payload["main_step"], "main_step")
    convergence = payload["step_convergence"]
    if not isinstance(convergence, dict):
        raise ValidationError("step_convergence must be an object")
    if set(convergence) != set(REQUIRED_STEPS):
        raise ValidationError("step_convergence must contain exactly 0.005, 0.01, 0.02 and 0.05")
    main_key = format(main_step, "g")
    if main_key not in convergence:
        raise ValidationError(f"main_step {main_step:g} is absent from step_convergence")

    convergence_report: dict[str, dict[str, Any]] = {}
    for step in REQUIRED_STEPS:
        row = convergence[step]
        if not isinstance(row, dict):
            raise ValidationError(f"step_convergence[{step}] must be an object")
        values = _vector(row.get("eigvalsh"), f"step_convergence[{step}].eigvalsh")
        if values != sorted(values):
            raise ValidationError(f"step_convergence[{step}].eigvalsh must be ascending")
        convergence_report[step] = {
            "eigvalsh_reported": values,
            "max_abs_delta_to_main": max(abs(value - reference) for value, reference in zip(values, reported)),
        }

    interpretation = payload["interpretation"]
    if not isinstance(interpretation, str) or "pilot" not in interpretation.lower() or "proof" not in interpretation.lower():
        raise ValidationError("interpretation must retain the pilot/non-proof warning")

    return {
        "schema": "fst-i.log-hessian-pilot-verification.v1",
        "status": "verified",
        "source_artifact": _display_path(input_path),
        "source_created": payload["created"],
        "project": payload["project"],
        "source_script": payload["source_script"],
        "objective": payload["objective"],
        "main_step": main_step,
        "H_log_main": matrix,
        "eigvalsh_main_reported": reported,
        "eigvalsh_main_recomputed": recomputed,
        "max_abs_eigenvalue_error": max(errors),
        "step_convergence": convergence_report,
        "checks": {
            "matrix_shape": True,
            "matrix_symmetric": True,
            "main_eigenvalues_match": True,
            "step_convergence_complete": True,
            "pilot_warning_present": True,
        },
        "interpretation": interpretation,
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="published pilot JSON (default: the repository artifact)")
    parser.add_argument("--output", type=Path, help="write the deterministic verification report here; stdout when omitted")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    try:
        with args.input.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict):
            raise ValidationError("input JSON must be an object")
        report = _validate(payload, args.input)
        rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
        else:
            sys.stdout.write(rendered)
        return 0
    except (OSError, json.JSONDecodeError, ValidationError) as error:
        print(f"verification failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
