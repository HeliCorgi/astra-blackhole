"""Independent finite-fraction check; never substitutes for Lean compilation.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
THEOREMS = [
    "p4_low_moments", "p4_first_visible",
    "p6_low_moments", "p6_first_visible",
    "p8_low_moments", "p8_first_visible",
    "equal_observations_equal_predictions",
    "no_exact_forecast_for_distinct_futures",
]


def check_sources() -> None:
    source = (ROOT / "lean/AstraBlackhole/AngularCertificates.lean").read_text()
    expected = {4: Fraction(8, 315), 6: Fraction(16, 3003), 8: Fraction(128, 109395)}
    for ell, first in expected.items():
        match = re.search(rf"def p{ell} : List Rat := \[(.*?)\]", source)
        if match is None:
            raise ValueError(f"Missing coefficient definition P{ell}")
        coeff = [Fraction(v.strip()) for v in match.group(1).split(",")]
        values = [sum((c / (2 * (s + k) + 1) for k, c in enumerate(coeff)), Fraction(0))
                  for s in range(ell // 2 + 1)]
        if any(values[:-1]) or values[-1] != first:
            raise ValueError(f"Invalid exact finite moments for P{ell}: {values}")
        print(f"P{ell}: finite rational moments {list(map(str, values))}")
    for path in (ROOT / "lean").rglob("*.lean"):
        if ".lake" in path.parts:
            continue
        text = re.sub(r"/-.*?-/", "", path.read_text(), flags=re.S)
        text = re.sub(r"--[^\n]*", "", text)
        if re.search(r"\b(sorry|admit|axiom|native_decide)\b", text):
            raise ValueError(f"Disallowed proof shortcut in {path}")
    print("Python fraction/source checks passed; this is NOT a Lean build result.")


def check_axioms(path: Path) -> None:
    lines = path.read_text().splitlines()
    for name in THEOREMS:
        full = "AstraBlackhole." + name
        matching = [line.strip() for line in lines if line.startswith("'" + full + "'")]
        if len(matching) != 1:
            raise ValueError(f"Missing or duplicate axiom report for {full}")
        if matching[0] not in [f"'{full}' does not depend on any axioms",
                               f"'{full}' depends on axioms: []"]:
            raise ValueError(f"Unexpected axiom dependencies: {matching[0]}")
    print(f"Axiom reports verified: {len(THEOREMS)} declarations, no axiom dependencies.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--axioms", type=Path)
    args = parser.parse_args()
    check_sources()
    if args.axioms:
        check_axioms(args.axioms)


if __name__ == "__main__":
    main()
