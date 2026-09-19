"""Regression checker for the KS continuum threshold/domain audit.

Scientific PENDING states are expected outputs, not CI failures.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--artifact",type=Path,required=True)
    p.add_argument("--baseline",type=Path,required=True)
    a=p.parse_args()
    got=json.loads(a.artifact.read_text())
    base=json.loads(a.baseline.read_text())

    assert got["status"]=="PARTIAL"
    assert got["exact_continuum_spectral_representation"]["status"]=="PASS"
    assert got["positive_frequency_KG_completion"]["status"]=="PASS"
    assert got["continuum_velocity"]["status"]=="PASS"
    assert got["mass_continuum_form"]["status"]=="PENDING"
    assert got["Y_clock_implication"]["status"]=="PENDING"

    ref=got["reference_packet_threshold"]
    if not ref["B0_nonzero"] or ref["relative_abs_B0_spread"]>1e-8:
        raise AssertionError(ref)
    tests={row["s"]: row for row in ref["H_minus_s_threshold_tests"]}
    if not tests[.5]["finite_at_zero_for_generic_nonzero_B0"]:
        raise AssertionError(tests)
    if not tests[1.0]["finite_at_zero_for_generic_nonzero_B0"]:
        raise AssertionError(tests)
    if tests[1.5]["finite_at_zero_for_generic_nonzero_B0"]:
        raise AssertionError(tests)

    rows=got["continuum_velocity"]["finite_box_controls"]
    if not all(r["sylvester_relative_residual"]<5e-11 for r in rows):
        raise AssertionError(rows)
    gaps=[r["min_I_plus_velocity_eigenvalue"] for r in rows]
    if not all(g>0 for g in gaps) or not all(b<a for a,b in zip(gaps,gaps[1:])):
        raise AssertionError(gaps)

    b0=ref["B0_controls"][1]["abs_B0"]
    b0_ref=base["reference_packet_threshold"]["B0_controls"][1]["abs_B0"]
    if abs(b0-b0_ref)>1e-10:
        raise AssertionError((b0,b0_ref))
    for r,br in zip(rows,base["continuum_velocity"]["finite_box_controls"]):
        if abs(r["min_I_plus_velocity_eigenvalue"]-br["min_I_plus_velocity_eigenvalue"])>2e-9:
            raise AssertionError((r,br))

    print("KS continuum threshold/domain regression passed; KG PASS and mass/Y-domain PENDING findings preserved.")


if __name__=="__main__":
    main()
