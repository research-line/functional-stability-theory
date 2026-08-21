"""
compute_u1_2d_strong_coupling_positive_control.py
==================================================
Yang-Mills / Station-6-Review Route YM-R1: minimal U(1) strong-coupling
positive control, drafted 2026-07-03. Read the header before wiring
anything from this script into compute_rp_os_rfep_transfer_ledger.py.

WHAT THIS DOES: 2D compact U(1) lattice gauge theory (Wilson action
S = -beta * sum_plaquette cos(theta_p)) via Metropolis Monte Carlo. Measures
Wilson loops W(R,R) for several R and fits the area-law string tension
sigma_MC = -ln<W(R,R)> / R^2. Compares against the EXACT analytic result for
2D compact U(1) (confines at every coupling, exactly solvable via
independent per-plaquette integration):

    sigma_exact(beta) = -ln( I_1(beta) / I_0(beta) )

I_0, I_1 = modified Bessel functions of the first kind. Source: WebSearch
2026-07-03, "Lattice Gauge Theory and Wilson-Loop Confinement: A
Statistical-Mechanical Survey" (arXiv:2605.02156) and standard lattice
gauge theory references (Kogut 1979-style strong-coupling expansion,
exact in 2D because each plaquette factorizes). NOT independently
re-derived here -- treat as a literature-sourced ground truth, re-verify
before using in a publication-grade claim.

WHAT THIS DOES NOT DO (explicit scope limit, 2026-07-03):
- Does NOT compute os_crossing_capacity, rp_os_cone_status,
  projected_residual_over_tangent_gap, or any other
  compute_rp_os_rfep_transfer_ledger.py CASES[] field. Those are
  RFEP-Tangential-Waterline-specific quantities with an operational
  definition this script's author (Claude, Station-6 follow-up session)
  did not verify from the project's own RFEP core documents
  (CORE/RFEP/). Wiring a real CASES[] row therefore requires someone
  who can define those fields for this toy model -- NOT attempted here,
  to avoid fabricating ledger fields just to fill the schema.
- Does NOT touch SU(2) or any non-abelian gauge group.
- 2D only. No claim about the 4D continuum Yang-Mills mass gap.

Status: diagnostic self-check of a Monte-Carlo implementation against a
literature-sourced exact result. Not a Yang-Mills proof step. claim_upgrade
is not applicable here because this script does not write into any
project claim register -- it only produces its own status JSON.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.special import iv  # modified Bessel function of the first kind


def exact_string_tension(beta: float) -> float:
    i0 = iv(0, beta)
    i1 = iv(1, beta)
    return -math.log(i1 / i0)


class U1Lattice2D:
    """Periodic L x L 2D lattice, one link angle per (site, direction)."""

    def __init__(self, length: int, beta: float, rng: np.random.Generator) -> None:
        self.length = length
        self.beta = beta
        self.rng = rng
        # theta[dir, x, y], dir=0 -> x-links, dir=1 -> y-links
        self.theta = rng.uniform(-math.pi, math.pi, size=(2, length, length))

    def plaquette(self, x: int, y: int) -> float:
        length = self.length
        xp, yp = (x + 1) % length, (y + 1) % length
        return (
            self.theta[0, x, y]
            + self.theta[1, xp, y]
            - self.theta[0, x, yp]
            - self.theta[1, x, y]
        )

    def total_action(self) -> float:
        length = self.length
        s = 0.0
        for x in range(length):
            for y in range(length):
                s += math.cos(self.plaquette(x, y))
        return -self.beta * s

    def metropolis_sweep(self) -> None:
        length = self.length
        for direction in (0, 1):
            for x in range(length):
                for y in range(length):
                    old = self.theta[direction, x, y]
                    delta_action = self._local_action_delta(direction, x, y, old)
                    proposal = old + self.rng.uniform(-math.pi, math.pi)
                    new_delta_action = self._local_action_delta(direction, x, y, proposal)
                    d_s = new_delta_action - delta_action
                    if d_s <= 0 or self.rng.random() < math.exp(-d_s):
                        self.theta[direction, x, y] = proposal

    def _local_action_delta(self, direction: int, x: int, y: int, value: float) -> float:
        # Sum of cos(plaquette) over the (at most 2) plaquettes touching this link,
        # holding the link at `value`, times -beta. Only the two plaquettes that
        # actually contain link (direction, x, y) are affected by a change here.
        length = self.length
        old_value = self.theta[direction, x, y]
        self.theta[direction, x, y] = value
        total = 0.0
        if direction == 0:
            total += math.cos(self.plaquette(x, y))
            total += math.cos(self.plaquette(x, (y - 1) % length))
        else:
            total += math.cos(self.plaquette(x, y))
            total += math.cos(self.plaquette((x - 1) % length, y))
        self.theta[direction, x, y] = old_value
        return -self.beta * total

    def wilson_loop(self, size: int) -> float:
        """R x R square Wilson loop (convenience wrapper)."""
        return self.wilson_loop_rect(size, size)

    def wilson_loop_rect(self, width: int, height: int) -> float:
        """width x height rectangular Wilson loop starting at a fixed anchor site."""
        length = self.length
        assert max(width, height) < length, "loop must be smaller than the lattice to avoid wraparound bias"
        x0, y0 = 0, 0
        angle_sum = 0.0
        for step in range(width):
            angle_sum += self.theta[0, (x0 + step) % length, y0]
        for step in range(height):
            angle_sum += self.theta[1, (x0 + width) % length, (y0 + step) % length]
        for step in range(width):
            angle_sum -= self.theta[0, (x0 + width - 1 - step) % length, (y0 + height) % length]
        for step in range(height):
            angle_sum -= self.theta[1, x0 % length, (y0 + height - 1 - step) % length]
        return math.cos(angle_sum)


def run(
    length: int,
    beta: float,
    n_thermalize: int,
    n_measure: int,
    measure_every: int,
    loop_sizes: tuple[int, ...],
    seed: int,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    lattice = U1Lattice2D(length, beta, rng)

    for _ in range(n_thermalize):
        lattice.metropolis_sweep()

    # Rectangular pairs needed for Creutz ratios chi(R) =
    # -ln[ W(R,R)*W(R-1,R-1) / W(R,R-1)^2 ], which cancel the leading
    # perimeter-law term and converge to sigma_exact as R grows -- standard
    # technique (Creutz 1980), not an invented method.
    rect_pairs: set[tuple[int, int]] = set()
    for size in loop_sizes:
        rect_pairs.add((size, size))
        if size - 1 >= 1:
            rect_pairs.add((size - 1, size - 1))
            rect_pairs.add((size, size - 1))

    rect_samples: dict[tuple[int, int], list[float]] = {pair: [] for pair in rect_pairs}
    for step in range(n_measure * measure_every):
        lattice.metropolis_sweep()
        if step % measure_every == 0:
            for width, height in rect_pairs:
                rect_samples[(width, height)].append(lattice.wilson_loop_rect(width, height))

    rect_means = {pair: float(np.mean(values)) for pair, values in rect_samples.items()}
    rect_stderr = {
        pair: float(np.std(values, ddof=1) / math.sqrt(len(values)))
        for pair, values in rect_samples.items()
    }
    loop_means = {size: rect_means[(size, size)] for size in loop_sizes}
    loop_stderr = {size: rect_stderr[(size, size)] for size in loop_sizes}

    # Naive area-law fit (kept for comparison, known to be perimeter-biased
    # at small R -- see notes.perimeter_term_not_separated).
    sizes = np.array(sorted(loop_sizes), dtype=float)
    neg_log_w = np.array([-math.log(max(loop_means[int(r)], 1e-12)) for r in sizes])
    areas = sizes**2
    a_matrix = np.vstack([areas, np.ones_like(areas)]).T
    slope, intercept = np.linalg.lstsq(a_matrix, neg_log_w, rcond=None)[0]
    sigma_mc_naive = float(slope)

    creutz_ratios: dict[int, float] = {}
    for size in sorted(loop_sizes):
        if size - 1 < 1:
            continue
        w_rr = max(rect_means[(size, size)], 1e-12)
        w_r1r1 = max(rect_means[(size - 1, size - 1)], 1e-12)
        w_rr1 = max(rect_means[(size, size - 1)], 1e-12)
        creutz_ratios[size] = -math.log((w_rr * w_r1r1) / (w_rr1 * w_rr1))

    # Best Creutz estimate: largest R available (closest to the R->infinity
    # asymptotic sigma), which is the standard convergence direction.
    sigma_mc_creutz = creutz_ratios[max(creutz_ratios)] if creutz_ratios else float("nan")

    sigma_exact = exact_string_tension(beta)
    delta_naive = abs(sigma_mc_naive - sigma_exact)
    delta_creutz = abs(sigma_mc_creutz - sigma_exact) if creutz_ratios else float("nan")
    relative_delta_naive = delta_naive / sigma_exact if sigma_exact > 0 else float("nan")
    relative_delta_creutz = (
        delta_creutz / sigma_exact if sigma_exact > 0 and creutz_ratios else float("nan")
    )

    return {
        "run_date": "2026-07-03",
        "job": "compute_u1_2d_strong_coupling_positive_control",
        "lattice": {"length": length, "beta": beta, "seed": seed},
        "mc_params": {
            "n_thermalize_sweeps": n_thermalize,
            "n_measurements": n_measure,
            "measure_every_sweeps": measure_every,
        },
        "wilson_loops": {
            str(size): {"mean": loop_means[size], "stderr": loop_stderr[size]}
            for size in loop_sizes
        },
        "area_law_fit_naive": {
            "sigma_mc": sigma_mc_naive,
            "intercept": float(intercept),
            "method": "least_squares_-ln<W(R,R)>_vs_R^2_no_perimeter_term",
            "known_bias": "perimeter-law contamination at small R, see notes",
        },
        "creutz_ratios": {
            str(size): value for size, value in creutz_ratios.items()
        },
        "creutz_estimate": {
            "sigma_mc_creutz": sigma_mc_creutz,
            "at_size_R": max(creutz_ratios) if creutz_ratios else None,
            "method": "chi(R) = -ln[W(R,R)*W(R-1,R-1) / W(R,R-1)^2] (Creutz 1980), perimeter-term cancels",
        },
        "exact_reference": {
            "sigma_exact": sigma_exact,
            "formula": "-ln(I1(beta)/I0(beta))",
            "source": (
                "WebSearch 2026-07-03: arXiv:2605.02156 'Lattice Gauge Theory and "
                "Wilson-Loop Confinement: A Statistical-Mechanical Survey'; "
                "standard result for 2D compact U(1) (exact, not just leading "
                "strong-coupling order). NOT independently re-derived by this script."
            ),
        },
        "comparison": {
            "delta_naive": delta_naive,
            "relative_delta_naive": relative_delta_naive,
            "delta_creutz": delta_creutz,
            "relative_delta_creutz": relative_delta_creutz,
            "matches_within_10_percent_creutz": (
                bool(relative_delta_creutz < 0.10) if not math.isnan(relative_delta_creutz) else False
            ),
        },
        "notes": {
            "perimeter_term_not_separated_in_naive_fit": (
                "area_law_fit_naive ignores a possible perimeter-law contribution "
                "(-ln<W> ~ sigma*Area + c*Perimeter); creutz_ratios/creutz_estimate "
                "is the standard correction and should be treated as the primary result."
            ),
            "not_wired_into_rp_os_rfep_ledger": (
                "See module docstring: os_crossing_capacity / rp_os_cone_status / "
                "projected_residual_over_tangent_gap etc. are NOT computed here."
            ),
        },
        "status": "done",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=16)
    parser.add_argument("--beta", type=float, default=3.0)
    parser.add_argument("--n-thermalize", type=int, default=300)
    parser.add_argument("--n-measure", type=int, default=1500)
    parser.add_argument("--measure-every", type=int, default=5)
    parser.add_argument("--loop-sizes", type=str, default="1,2,3,4")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", type=str, required=True)
    args = parser.parse_args()

    loop_sizes = tuple(int(s) for s in args.loop_sizes.split(","))
    payload = run(
        length=args.length,
        beta=args.beta,
        n_thermalize=args.n_thermalize,
        n_measure=args.n_measure,
        measure_every=args.measure_every,
        loop_sizes=loop_sizes,
        seed=args.seed,
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {out_path}")
    print(f"sigma_mc_naive={payload['area_law_fit_naive']['sigma_mc']:.6f} "
          f"sigma_mc_creutz={payload['creutz_estimate']['sigma_mc_creutz']:.6f} "
          f"sigma_exact={payload['exact_reference']['sigma_exact']:.6f} "
          f"relative_delta_creutz={payload['comparison']['relative_delta_creutz']:.4f}")


if __name__ == "__main__":
    main()
