"""Same-state map between the repo L2 positive-frequency branch and KG form.

For the time-independent positive operator H=sqrt(A), the positive-frequency
KG norm is (2/h)<psi,H phi>_L2.  The map
    psi_KG = sqrt(h/2) H^{-1/2} chi_L2
is an isometry on its domain.  This script verifies that statement on the
finite-box discretization and keeps the continuum H^{-1/2} domain separate.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))
from wdw_model import Packet,WDWModel

def mapped_observable_expectation(chi,psi,model,diag_values):
    # S=sqrt(h/2) H^-1/2; S^{-1}=sqrt(2/h) H^1/2.
    U=model.U
    E=model.energy
    h=model.h
    def S(v):
        return math.sqrt(h/2)*(U@((U.T@v)/np.sqrt(E)))
    # Field amplitude representing O chi.
    phi=S(diag_values*chi)
    # KG inner product on positive-frequency initial data.
    return (2/h)*np.vdot(psi,U@(E*(U.T@phi)))

def run_case(h):
    sigma=math.sqrt(h)
    model=WDWModel(h,left=-4.,right=16.,dx=.04)
    packet=Packet(h,sigma)
    chi0=packet.initial(model.x)
    U=model.U; E=model.energy
    def S(v):
        return math.sqrt(h/2)*(U@((U.T@v)/np.sqrt(E)))
    def Sinv(v):
        return math.sqrt(2/h)*(U@(np.sqrt(E)*(U.T@v)))
    psi0=S(chi0)
    kg0=(2/h)*np.vdot(psi0,U@(E*(U.T@psi0)))
    times=[0.,2.,6.]
    rows=[]
    max_back=0.; max_kg=0.; max_x=0.; max_r=0.; max_naive_x=0.
    coeff=U.T@chi0
    for T in times:
        phase=np.exp(-1j*E*T/h)
        chi=U@(coeff*phase)
        psi=S(chi)
        back=Sinv(psi)
        kg=(2/h)*np.vdot(psi,U@(E*(U.T@psi)))
        x_l2=np.vdot(chi,model.x*chi)
        x_kg=mapped_observable_expectation(chi,psi,model,model.x)
        rdiag=.25*np.exp(-model.x-T)
        r_l2=np.vdot(chi,rdiag*chi)
        r_kg=mapped_observable_expectation(chi,psi,model,rdiag)
        naive_prob=np.abs(psi)**2
        naive_prob/=naive_prob.sum()
        x_naive=float(model.x@naive_prob)
        max_back=max(max_back,float(np.linalg.norm(back-chi)))
        max_kg=max(max_kg,float(abs(kg-1)))
        max_x=max(max_x,float(abs(x_kg-x_l2)))
        max_r=max(max_r,float(abs(r_kg-r_l2)))
        max_naive_x=max(max_naive_x,float(abs(x_naive-x_l2.real)))
        rows.append({
          "T":T,
          "kg_norm_real":float(kg.real),
          "map_back_l2_error":float(np.linalg.norm(back-chi)),
          "x_l2":float(x_l2.real),
          "x_kg_mapped":float(x_kg.real),
          "r_l2":float(r_l2.real),
          "r_kg_mapped":float(r_kg.real),
          "naive_same_field_x":x_naive
        })
    return {
      "h":h,
      "min_finite_box_energy":float(E.min()),
      "initial_kg_norm_real":float(kg0.real),
      "max_map_back_l2_error":max_back,
      "max_kg_norm_error":max_kg,
      "max_mapped_x_expectation_error":max_x,
      "max_mapped_r_expectation_error":max_r,
      "max_naive_same_field_x_difference":max_naive_x,
      "rows":rows
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()

    cases=[run_case(h) for h in (.4,.2,.1)]
    if max(x["max_map_back_l2_error"] for x in cases)>2e-11:
        raise RuntimeError("L2/KG map-back failed")
    if max(x["max_kg_norm_error"] for x in cases)>2e-11:
        raise RuntimeError("KG norm not preserved by same-state map")
    if max(x["max_mapped_x_expectation_error"] for x in cases)>2e-11:
        raise RuntimeError("mapped x expectation mismatch")
    if max(x["max_mapped_r_expectation_error"] for x in cases)>2e-11:
        raise RuntimeError("mapped r expectation mismatch")

    # Finite-box evidence that H^{-1/2} becomes increasingly singular as the
    # right wall is moved outward. This is diagnostic, not a continuum proof.
    threshold=[]
    for right in (16.,24.,32.):
        m=WDWModel(.2,left=-4.,right=right,dx=.04)
        threshold.append({
          "right":right,
          "min_energy":float(m.energy.min()),
          "inverse_sqrt_min_energy":float(1/math.sqrt(m.energy.min()))
        })

    chiba=json.loads((ROOT/"research/chiba_results/summary.json").read_text())
    plateau=chiba["plateau"]
    direct_chiba_l2_status="FAIL"
    if plateau["squared"]<=0:
        raise RuntimeError("expected nonzero Chiba pullback plateau was not present")

    result={
      "schema":1,
      "status":"PARTIAL",
      "finite_box_positive_frequency_same_state_map":"PASS",
      "continuum_H_minus_half_domain":"PENDING",
      "chiba_scalar_pullback_direct_L2_identification":direct_chiba_l2_status,
      "map":{
        "L2_to_KG":"psi=sqrt(h/2) H^(-1/2) chi",
        "KG_to_L2":"chi=sqrt(2/h) H^(1/2) psi",
        "positive_frequency_KG_inner_product":"(psi,phi)_KG=(2/h)<psi,H phi>_L2",
        "observable_map":"O_KG=S O_L2 S^(-1), not the same multiplication operator in general"
      },
      "cases":cases,
      "threshold_diagnostic":threshold,
      "chiba_negative_control":{
        "pullback_constant_tail_squared":plateau["squared"],
        "source_conclusion":plateau["conclusion"],
        "interpretation":"The published KG scalar pullback is not directly the same state as the repo L2 Gaussian. Same-state comparison requires a sector choice and the nonlocal H^(±1/2) map."
      },
      "scope_note":"The finite-box positive-frequency sectors are isometric under the explicit map. Continuum H^{-1/2} is unbounded at the zero-energy threshold, so its domain/completion must be specified before claiming a global L2<->KG equivalence."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
