"""Regression checker for the KS affine-mass obstruction audit.

Scientific PENDING/selection-failure states are expected outputs, not CI
failures.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--artifact",type=Path,required=True)
    p.add_argument("--baseline",type=Path,required=True)
    a=p.parse_args()
    got=json.loads(a.artifact.read_text())
    base=json.loads(a.baseline.read_text())

    assert got["status"]=="PARTIAL"
    assert got["classical_affine_structure"]["status"]=="PASS"
    assert got["null_coordinate_bridge"]["status"]=="PASS"
    assert got["quantum_affine_covariance_no_go"]["status"]=="PASS"
    assert got["mass_continuum_form"]["status"]=="PENDING"
    assert got["ordering_selection"]["status"]=="FAIL_FOR_EXACT_AFFINE_COVARIANCE_CRITERION"

    rows=got["declared_mass_candidate_negative_control"]["finite_box_controls"]
    brows=base["declared_mass_candidate_negative_control"]["finite_box_controls"]
    if not all(r["relative_residual_for_[H,M0]=-ihM0"]>0.1 for r in rows):
        raise AssertionError(rows)
    if not all(r["commutator_antihermiticity_residual"]<2e-12 for r in rows):
        raise AssertionError(rows)
    for r,b in zip(rows,brows):
        if abs(r["relative_residual_for_[H,M0]=-ihM0"]-b["relative_residual_for_[H,M0]=-ihM0"])>2e-9:
            raise AssertionError((r,b))

    print("KS affine-mass obstruction regression passed; exact-affine selection no-go PASS and mass closability PENDING preserved.")


if __name__=="__main__":
    main()
