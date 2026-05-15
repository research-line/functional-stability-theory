"""Finite toy check for the Mobius-graded Prime-Hub identity.

This is not a Prime-Hub construction. It verifies, on finite prime
sets, the algebraic core discussed in the Zeta Zoo v2 paper and writes
the public results in _results/MOBIUS_GRADED_PRIMEHUB_TOY.{md,json}.

For a finite prime set P_N it compares:
  1. Product determinant prod_p (1 - p^-s)
  2. Exterior-algebra supertrace sum_S (-1)^|S| prod_{p in S} p^-s
  3. Squarefree Mobius sum sum_n mu(n)n^-s over n supported on P_N
"""

from __future__ import annotations

import cmath
import itertools
import json
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "_results"
PRIMES = [2, 3, 5, 7, 11]
S_VALUES = [
    ("s=2", 2 + 0j),
    ("s=3", 3 + 0j),
    ("s=0.75+5i", 0.75 + 5j),
]


def cformat(z: complex, digits: int = 12) -> str:
    if abs(z.imag) < 10 ** (-digits + 2):
        return f"{z.real:.{digits}g}"
    sign = "+" if z.imag >= 0 else "-"
    return f"{z.real:.{digits}g}{sign}{abs(z.imag):.{digits}g}i"


def cdict(z: complex) -> dict[str, float]:
    return {"real": z.real, "imag": z.imag}


def product_determinant(primes: Iterable[int], s: complex) -> complex:
    out = 1 + 0j
    for p in primes:
        out *= 1 - cmath.exp(-s * cmath.log(p))
    return out


def exterior_supertrace(primes: list[int], s: complex) -> complex:
    out = 0 + 0j
    for k in range(len(primes) + 1):
        sign = -1 if k % 2 else 1
        for subset in itertools.combinations(primes, k):
            n = 1
            for p in subset:
                n *= p
            out += sign * cmath.exp(-s * cmath.log(n))
    return out


def squarefree_terms(primes: list[int]) -> list[dict[str, object]]:
    terms = []
    for k in range(len(primes) + 1):
        sign = -1 if k % 2 else 1
        for subset in itertools.combinations(primes, k):
            n = 1
            for p in subset:
                n *= p
            terms.append({"n": n, "mu": sign, "support": list(subset)})
    return sorted(terms, key=lambda item: int(item["n"]))


def squarefree_mobius_sum(terms: list[dict[str, object]], s: complex) -> complex:
    out = 0 + 0j
    for term in terms:
        out += int(term["mu"]) * cmath.exp(-s * cmath.log(int(term["n"])))
    return out


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    terms = squarefree_terms(PRIMES)

    rows = []
    for label, s in S_VALUES:
        det = product_determinant(PRIMES, s)
        supertrace = exterior_supertrace(PRIMES, s)
        mobius_sum = squarefree_mobius_sum(terms, s)
        positive_graph_zeta = 1 / det
        max_abs_error = max(
            abs(det - supertrace),
            abs(det - mobius_sum),
            abs(supertrace - mobius_sum),
        )
        rows.append(
            {
                "label": label,
                "s": cdict(s),
                "product_determinant": cdict(det),
                "exterior_supertrace": cdict(supertrace),
                "squarefree_mobius_sum": cdict(mobius_sum),
                "positive_graph_zeta_inverse_orientation": cdict(positive_graph_zeta),
                "max_abs_error": max_abs_error,
            }
        )

    data = {
        "prime_set": PRIMES,
        "n_squarefree_terms": len(terms),
        "terms": terms,
        "checks": rows,
        "interpretation": (
            "Finite product determinant, exterior supertrace, and squarefree "
            "Mobius sum agree up to floating-point error. The reciprocal is "
            "the positive graph-zeta orientation and has the opposite pole/zero behavior."
        ),
    }

    json_path = RESULTS / "MOBIUS_GRADED_PRIMEHUB_TOY.json"
    json_path.write_text(json.dumps(data, indent=2), encoding="utf-8")

    lines = [
        "# MOBIUS_GRADED_PRIMEHUB_TOY",
        "",
        "**Purpose.** Finite check of the identity",
        "",
        "`det(I-A_s) = Str(Lambda A_s) = sum mu(n)n^-s = prod_p(1-p^-s)`.",
        "",
        f"**Prime set:** `{PRIMES}`",
        f"**Squarefree terms:** `{len(terms)}`",
        "",
        "| s | product determinant | exterior supertrace | squarefree Mobius sum | max abs error | reciprocal orientation |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        det = complex(row["product_determinant"]["real"], row["product_determinant"]["imag"])
        supertrace = complex(row["exterior_supertrace"]["real"], row["exterior_supertrace"]["imag"])
        mobius_sum = complex(row["squarefree_mobius_sum"]["real"], row["squarefree_mobius_sum"]["imag"])
        reciprocal = complex(
            row["positive_graph_zeta_inverse_orientation"]["real"],
            row["positive_graph_zeta_inverse_orientation"]["imag"],
        )
        lines.append(
            "| {label} | `{det}` | `{supertrace}` | `{mobius_sum}` | `{err:.3e}` | `{reciprocal}` |".format(
                label=row["label"],
                det=cformat(det),
                supertrace=cformat(supertrace),
                mobius_sum=cformat(mobius_sum),
                err=row["max_abs_error"],
                reciprocal=cformat(reciprocal),
            )
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The three left columns agree to floating-point precision. This is the finite",
            "exterior-algebra core of the Möbius-graded Prime-Hub route.",
            "",
            "The final column is the reciprocal orientation. It corresponds to the",
            "ordinary positive Euler-product zeta direction. This is why the OP5",
            "orientation must be audited before looking for eigenvalue-one events",
            "at zeros of `zeta`: `1/zeta` has poles where `zeta` has zeros.",
            "",
            "## First squarefree terms",
            "",
            "| n | mu | support |",
            "|---:|---:|---|",
        ]
    )
    for term in terms[:16]:
        lines.append(f"| {term['n']} | {term['mu']} | `{term['support']}` |")

    md_path = RESULTS / "MOBIUS_GRADED_PRIMEHUB_TOY.md"
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
