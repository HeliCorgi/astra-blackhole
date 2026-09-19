"""Regression checker for the KS physical mass/Y-clock/Kretschmann audit.

Scientific FAIL/PARTIAL states are expected outputs, not CI failures.

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
    mass=json.loads((a.artifact_dir/"mass_dirac.json").read_text())
    y=json.loads((a.artifact_dir/"y_clock.json").read_text())
    k=json.loads((a.artifact_dir/"kretschmann.json").read_text())

    assert mass["status"]=="PARTIAL"
    assert mass["finite_box_positive_mass_dirac_candidate"]=="PASS"
    assert mass["finite_box_physical_KG_adjointness"]=="PASS"
    assert mass["continuum_quadratic_form_closure"]=="PENDING"
    assert mass["semiclassical_error_decreases_with_h"]
    errs=mass["semiclassical_same_wigner_relative_error_sequence"]
    if not all(b<a for a,b in zip(errs,errs[1:])) or errs[-1]>.05:
        raise AssertionError(errs)
    max_drift=max(c["max_relative_mass_expectation_drift"] for c in mass["cases"])
    max_kg=max(c["max_KG_adjoint_bilinear_residual"] for c in mass["cases"])
    if max_drift>2e-11 or max_kg>1e-10:
        raise AssertionError((max_drift,max_kg))
    box=[r["mass"] for r in mass["box_diagnostic"]]
    if (max(box)-min(box))/abs(box[1]) > 5e-6:
        raise AssertionError(box)
    if errs[-1]>.02 or errs[1]>.08:
        raise AssertionError(("same-Wigner semiclassical sequence too large",errs))

    assert y["status"]=="PARTIAL"
    assert y["same_Dirac_state_null_flux"]=="PASS"
    assert y["same_mass_observable_null_flux"]=="PASS"
    assert y["Y_clock_deparametrized_inverse_pX_domain"]=="PENDING"
    flux=y["same_state_flux"]
    if flux["max_null_norm_difference"]>5e-4 or flux["max_null_mass_difference"]>1e-5:
        raise AssertionError(flux)
    dom=y["inverse_pX_threshold_diagnostic"]
    mins=[r["min_pX_eigenvalue"] for r in dom]
    inv=[r["max_inverse_pX_eigenvalue"] for r in dom]
    if not all(b<a for a,b in zip(mins,mins[1:])):
        raise AssertionError(mins)
    if not all(b>a for a,b in zip(inv,inv[1:])) or inv[-1]<1e6:
        raise AssertionError(inv)

    assert k["status"]=="FAIL"
    assert k["regulator_decision"]=="FAIL"
    assert k["finding"]=="selected_mass_candidate_does_not_cure_explicit_radius_domain"
    assert k["mass_applied_threshold_overlap_status"]=="NONZERO_NUMERICALLY_RESOLVED"
    for key,row in k["evolved_T2_growth"].items():
        assert row["monotone_increasing"]
        if row["total_log10_growth"]<50:
            raise AssertionError((key,row))
    ovs=k["mass_applied_threshold_overlap"]
    if ovs[-1]["abs"]<.1:
        raise AssertionError(ovs)
    rel=abs(ovs[-1]["abs"]/ovs[-2]["abs"]-1)
    if rel>.01:
        raise AssertionError(("threshold overlap refinement",rel))

    print("KS physical mass/Y-clock/Kretschmann regression passed; scoped PARTIAL/FAIL findings preserved.")

if __name__=="__main__":
    main()
