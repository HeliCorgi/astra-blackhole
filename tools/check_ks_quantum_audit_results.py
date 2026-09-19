"""Regression checker for the KS black-hole quantum audit artifacts.

Negative/PARTIAL scientific findings are expected outputs, not CI failures.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--artifact-dir",type=Path,required=True)
    p.add_argument("--baseline",type=Path,required=True)
    a=p.parse_args()
    base=json.loads(a.baseline.read_text())
    got={name:json.loads((a.artifact_dir/f"{name}.json").read_text())
         for name in ("clock","inner_product","factor_ordering","observable_domain")}

    assert got["clock"]["status"]=="PARTIAL"
    assert got["clock"]["classical_T_clock"]=="PASS"
    assert got["clock"]["classical_areal_clock_Y"]=="PASS"
    assert got["clock"]["quantum_multi_clock_equivalence"]=="PENDING"
    assert "inverse p_X" in got["clock"]["alternative_clock"]["quantum_domain_warning"]

    inner=got["inner_product"]
    assert inner["status"]=="PARTIAL"
    assert inner["finite_box_positive_frequency_same_state_map"]=="PASS"
    assert inner["continuum_H_minus_half_domain"]=="PENDING"
    maxmap=max(row["max_map_back_l2_error"] for row in inner["cases"])
    maxkg=max(row["max_kg_norm_error"] for row in inner["cases"])
    maxx=max(row["max_mapped_x_expectation_error"] for row in inner["cases"])
    maxr=max(row["max_mapped_r_expectation_error"] for row in inner["cases"])
    if max(maxmap,maxkg,maxx,maxr)>2e-11:
        raise AssertionError((maxmap,maxkg,maxx,maxr))
    mins=[row["min_energy"] for row in inner["threshold_diagnostic"]]
    if not all(b<a for a,b in zip(mins,mins[1:])):
        raise AssertionError("finite-box threshold energy did not decrease as right wall expanded")

    ordering=got["factor_ordering"]
    assert ordering["status"]=="FAIL" and ordering["finding"]=="sensitive"
    if ordering["max_mean_x_ordering_spread"] <= 100*ordering["existing_numerical_detectability_threshold"]:
        raise AssertionError("ordering effect no longer resolved over numerical tolerance")
    np.testing.assert_allclose(
        ordering["max_mean_x_ordering_spread"],
        base["factor_ordering"]["max_mean_x_ordering_spread"],
        rtol=2e-10,atol=2e-12)

    obs=got["observable_domain"]
    assert obs["status"]=="PARTIAL"
    assert obs["constraint_kernel_ordering_family"]["residual"]=="0"
    assert obs["flat_kinematical_L2_symmetry"]["status"]=="NO_SOLUTION"
    assert not obs["decision"]["quantum_mass_operator_selected"]
    assert not obs["decision"]["quantum_kretschmann_operator_selected"]
    fine=[row for row in obs["finite_box_ordering_scan"] if row["dx"]==.02][0]
    wd=fine["rows"]["weyl6"]["relative_expectation_drift"]
    sd=fine["rows"]["sym"]["relative_expectation_drift"]
    if not (0.005<wd<0.05 and 0.02<sd<0.08):
        raise AssertionError((wd,sd))
    np.testing.assert_allclose(
        wd,base["mass_curvature_observable"]["finite_box_dx02"]["weyl6_mass_expectation_drift"],
        rtol=2e-9,atol=2e-11)

    print("KS quantum audit regression passed; negative and partial gates preserved.")

if __name__=="__main__":
    main()
