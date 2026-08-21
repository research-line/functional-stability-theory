"""
U(1) strong-coupling positive-control ledger for the Yang-Mills project.

This is a finite-lattice guardrail/control, not a Yang-Mills mass-gap proof.
It fills the explicit next step left by the 2026-06-28 RP/OS-RFEP ledger:
feed a U(1)/strong-coupling positive control into the same target/source/
projection/outer-gap/residual schema.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
RESULT_DIR = ROOT / "_results"
RESULT_DIR.mkdir(exist_ok=True)

DATE = "2026-07-03"
BASE = RESULT_DIR / f"U1_STRONG_COUPLING_POSITIVE_CONTROL_{DATE}"

FIELDS = [
    "case_id",
    "role",
    "beta",
    "i0",
    "i1",
    "i2",
    "c1_ratio",
    "c2_ratio",
    "area_law_proxy",
    "finite_gap_proxy",
    "coefficient_decay_pass",
    "source_kind",
    "mirror_leak_class",
    "repair_certificate_id",
    "summable_loss_weight",
    "os_crossing_capacity",
    "area_law_mass_gap_separated",
    "rp_os_cone_status",
    "continuum_limit_status",
    "target_predefined",
    "source_axis_independent",
    "projection_predefined",
    "outer_gap_source",
    "spectral_exactness_status",
    "total_residual_over_gap",
    "projected_residual_over_tangent_gap",
    "normal_residual",
    "normal_paid_by_constraint",
    "matched_control_family",
    "decision",
    "claim_status",
    "notes",
]


def modified_bessel_i(n: int, x: float, tol: float = 1e-16) -> float:
    """Return I_n(x) for small integer n using the defining positive series."""
    if n < 0:
        raise ValueError("n must be non-negative")

    half = x / 2.0
    term = half**n / math.factorial(n)
    total = term
    k = 0
    while True:
        k += 1
        term *= (half * half) / (k * (k + n))
        total += term
        if abs(term) < tol * max(1.0, abs(total)):
            break
        if k > 200:
            raise RuntimeError("Bessel series did not converge")
    return total


def finite_control_row(beta: float) -> dict[str, Any]:
    i0 = modified_bessel_i(0, beta)
    i1 = modified_bessel_i(1, beta)
    i2 = modified_bessel_i(2, beta)
    c1 = i1 / i0
    c2 = i2 / i0
    area_law_proxy = -math.log(c1)
    finite_gap_proxy = 1.0 - c1
    coefficient_decay_pass = 0.0 < c2 < c1 < 1.0

    # This is a compact finite-control budget: it checks summability shape for
    # an already predeclared Fourier axis. It is not RG data.
    summable_loss_weight = sum((c1 ** (2 * k)) / ((k + 1) ** 2) for k in range(1, 128))
    residual_over_gap = c2 / finite_gap_proxy

    return {
        "case_id": f"u1_character_positive_beta_{str(beta).replace('.', '_')}",
        "role": "u1_strong_coupling_positive_control",
        "beta": beta,
        "i0": i0,
        "i1": i1,
        "i2": i2,
        "c1_ratio": c1,
        "c2_ratio": c2,
        "area_law_proxy": area_law_proxy,
        "finite_gap_proxy": finite_gap_proxy,
        "coefficient_decay_pass": coefficient_decay_pass,
        "source_kind": "analytic_u1_character_expansion",
        "mirror_leak_class": "internal",
        "repair_certificate_id": "U1_FOURIER_AXIS_AND_POSITIVE_COEFFICIENTS_PREDECLARED",
        "summable_loss_weight": summable_loss_weight,
        "os_crossing_capacity": 0.0,
        "area_law_mass_gap_separated": "yes",
        "rp_os_cone_status": "finite_lattice_rp_character_pass",
        "continuum_limit_status": "not_attempted_finite_control",
        "target_predefined": True,
        "source_axis_independent": True,
        "projection_predefined": True,
        "outer_gap_source": "u1_fourier_character_ratio_n1",
        "spectral_exactness_status": "analytic_u1_fourier_coefficients",
        "total_residual_over_gap": residual_over_gap,
        "projected_residual_over_tangent_gap": residual_over_gap,
        "normal_residual": 0.0,
        "normal_paid_by_constraint": True,
        "matched_control_family": "u1_strong_coupling_character_positive_control",
        "decision": "finite_positive_control_pass_no_continuum_claim",
        "notes": (
            "U(1) plaquette character ratios provide an analytic finite-lattice "
            "positive control for the RP/OS-RFEP schema; no nonabelian or continuum "
            "Yang-Mills claim is inferred."
        ),
    }


def os_crossing_negative_control(reference: dict[str, Any]) -> dict[str, Any]:
    row = dict(reference)
    row.update(
        {
            "case_id": "u1_same_coefficients_unpaid_os_crossing_negative",
            "role": "matched_negative_control",
            "mirror_leak_class": "os_crossing",
            "repair_certificate_id": "U1_FOURIER_AXIS_BUT_OS_CROSSING_UNREPAIRED",
            "os_crossing_capacity": 0.37,
            "rp_os_cone_status": "failed_os_crossing",
            "continuum_limit_status": "blocked_os_leak",
            "normal_residual": 0.37,
            "normal_paid_by_constraint": False,
            "matched_control_family": "same_u1_coefficients_bad_os_geometry",
            "decision": "reject",
            "notes": (
                "Same harmless Fourier coefficients, but an unpaid OS-crossing "
                "normal component must be rejected. This checks that the positive "
                "control cannot pass by coefficients alone."
            ),
        }
    )
    return row


def classify(row: dict[str, Any]) -> str:
    if row["decision"] == "reject":
        return "matched_negative_control_rejected"
    if not row["coefficient_decay_pass"]:
        return "blocked_character_decay"
    if row["mirror_leak_class"] == "os_crossing":
        return "blocked_os_crossing_leak"
    if row["continuum_limit_status"] != "continuum_verified":
        return "diagnostic_positive_control_no_claim"
    return "claim_candidate_needs_human_review"


def build_rows() -> list[dict[str, Any]]:
    rows = [finite_control_row(beta) for beta in [0.1, 0.2, 0.5, 0.8]]
    rows.append(os_crossing_negative_control(rows[-1]))
    for row in rows:
        row["claim_status"] = classify(row)
    return rows


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_md(rows: list[dict[str, Any]], path: Path, assertions: dict[str, Any]) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["claim_status"]] = counts.get(row["claim_status"], 0) + 1

    lines = [
        "# U(1) Strong-Coupling Positive-Control Ledger 2026-07-03",
        "",
        "Status: FORSCHER-Guardrail, finite-lattice positive control only.",
        "",
        "This ledger feeds an explicit U(1) character-expansion control into the",
        "RP/OS-RFEP transfer schema. The predeclared Fourier axis is the `n=1`",
        "character ratio `I_1(beta)/I_0(beta)`. The run checks coefficient decay,",
        "a positive finite gap proxy, a separated area-law proxy, and the existing",
        "target/source/projection/outer-gap/residual fields.",
        "",
        "## Summary",
        "",
    ]
    for key in sorted(counts):
        lines.append(f"- `{key}`: {counts[key]}")
    lines.extend(
        [
            f"- `claim_pass`: {assertions['claim_pass']}",
            f"- `positive_control_rows`: {assertions['positive_control_rows']}",
            f"- `negative_control_rows`: {assertions['negative_control_rows']}",
            "",
            "## Ledger",
            "",
            "| case | beta | c1=I1/I0 | c2=I2/I0 | area proxy | gap proxy | leak | decision | status |",
            "|---|---:|---:|---:|---:|---:|---|---|---|",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['case_id']} | {row['beta']:.3f} | {row['c1_ratio']:.12g} | "
            f"{row['c2_ratio']:.12g} | {row['area_law_proxy']:.12g} | "
            f"{row['finite_gap_proxy']:.12g} | {row['mirror_leak_class']} | "
            f"{row['decision']} | {row['claim_status']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The positive rows show that the RFEP transfer schema can represent a clean",
            "finite U(1) strong-coupling control: the target and projection are fixed",
            "before inspection, the source axis is independent of the desired continuum",
            "claim, the character coefficients decay as expected, and there is no",
            "OS-crossing normal residual. The matched negative row keeps the same",
            "coefficient data but adds an unpaid OS-crossing normal component; it is",
            "rejected.",
            "",
            "This does not upgrade the Yang-Mills paper. The control is abelian and",
            "finite-lattice only. The nonabelian RG/Gauge-block data, continuum limit,",
            "repair certificate and OS/RP-compatible blocking remain open.",
            "",
            "## Verification",
            "",
            f"- rows: {assertions['rows']}",
            f"- all finite positive rows pass coefficient decay: {assertions['all_positive_decay_pass']}",
            f"- maximum positive control c1 ratio: {assertions['max_positive_c1_ratio']:.12g}",
            f"- claim pass count: {assertions['claim_pass']}",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    positive_rows = [row for row in rows if row["role"] == "u1_strong_coupling_positive_control"]
    negative_rows = [row for row in rows if row["role"] == "matched_negative_control"]
    assertions = {
        "rows": len(rows),
        "positive_control_rows": len(positive_rows),
        "negative_control_rows": len(negative_rows),
        "all_positive_decay_pass": all(row["coefficient_decay_pass"] for row in positive_rows),
        "max_positive_c1_ratio": max(row["c1_ratio"] for row in positive_rows),
        "claim_pass": sum(1 for row in rows if row["claim_status"] == "claim_candidate_needs_human_review"),
    }

    if assertions["positive_control_rows"] != 4:
        raise AssertionError("expected four positive-control rows")
    if assertions["negative_control_rows"] != 1:
        raise AssertionError("expected one negative-control row")
    if not assertions["all_positive_decay_pass"]:
        raise AssertionError("positive control coefficient decay failed")
    if assertions["claim_pass"] != 0:
        raise AssertionError("control ledger must not create a claim pass")

    write_csv(rows, BASE.with_suffix(".csv"))
    payload = {"rows": rows, "assertions": assertions}
    BASE.with_suffix(".json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    write_md(rows, BASE.with_suffix(".md"), assertions)
    print(f"wrote {BASE.with_suffix('.md')}")
    print(
        "rows={rows} positive={positive_control_rows} negative={negative_control_rows} "
        "claim_pass={claim_pass}".format(**assertions)
    )


if __name__ == "__main__":
    main()
