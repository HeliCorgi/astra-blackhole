"""Independent finite-fraction checks and explicit Lean axiom-dependency policy.

Python checks never substitute for Lean compilation. The finite certificates
use only Lean's standard propext axiom; the observation lemmas use no axioms.
Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from pathlib import Path
import re
import tempfile

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_AXIOMS = {
    "p4_low_moments": {"propext"}, "p4_first_visible": {"propext"},
    "p6_low_moments": {"propext"}, "p6_first_visible": {"propext"},
    "p8_low_moments": {"propext"}, "p8_first_visible": {"propext"},
    "equal_observations_equal_predictions": set(),
    "no_exact_forecast_for_distinct_futures": set(),
    "time_dependent_factorization_residual": set(),
    "common_operator_preserves_recombination": set(),
    "commuting_operator_preserves_constraint_kernel": set(),
}


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


def check_axioms(path: Path, *, verbose: bool = True) -> None:
    lines = [line.strip() for line in path.read_text().splitlines()]
    for name, expected in EXPECTED_AXIOMS.items():
        full = "AstraBlackhole." + name
        matching = [line for line in lines if line.startswith("'" + full + "'")]
        if len(matching) != 1:
            raise ValueError(f"Missing or duplicate axiom report for {full}")
        line = matching[0]
        if line == f"'{full}' does not depend on any axioms":
            actual = set()
        else:
            match = re.fullmatch(re.escape(f"'{full}' depends on axioms: ") + r"\[([^\]]*)\]", line)
            if match is None:
                raise ValueError(f"Unrecognized axiom report: {line}")
            actual = {item.strip() for item in match.group(1).split(",") if item.strip()}
        if actual != expected:
            raise ValueError(f"Unexpected axioms for {full}: {actual}; expected {expected}")
    if verbose:
        print("Eleven reports verified: six use only standard propext; five use no axioms.")


def self_test() -> None:
    """Test the parser with synthetic reports, not as proof evidence."""
    lines = [f"'AstraBlackhole.{name}' depends on axioms: [{', '.join(sorted(axioms))}]"
             for name, axioms in EXPECTED_AXIOMS.items()]
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "axioms.txt"
        path.write_text("\n".join(lines))
        check_axioms(path, verbose=False)
        mutations = [lines[:-1], lines + [lines[0]],
                     [line.replace("[propext]", "[]") for line in lines],
                     [line.replace("[]", "[propext]") for line in lines]]
        for bad in ("sorryAx", "Lean.ofReduceBool", "Classical.choice", "Quot.sound", "FakeAxiom"):
            mutations.append([lines[0].replace("[propext]", f"[propext, {bad}]")] + lines[1:])
        for report in mutations:
            path.write_text("\n".join(report))
            try:
                check_axioms(path, verbose=False)
            except ValueError:
                continue
            raise AssertionError("Invalid synthetic report passed the axiom policy")
    print(f"Axiom parser self-test passed: one valid and {len(mutations)} rejected fixtures.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--axioms", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    check_sources()
    if args.self_test:
        self_test()
    if args.axioms:
        check_axioms(args.axioms)


if __name__ == "__main__":
    main()
