"""Kretschmann ordering/domain audit using the declared mass Dirac candidate.

Classically K=48 mu^2/r^6.  With the positive-frequency transported mass
candidate M_D(T), this script tests three positive quadratic-form orderings:

  K_MR = 48 || r^-3 M_D psi ||^2
  K_RM = 48 || M_D r^-3 psi ||^2
  K_S  = 48 || r^-3/2 M_D r^-3/2 psi ||^2

All have the same classical principal symbol 48 mu^2/r^6, but different
operator domains.  Finite boxes are regulators, not physical boundaries.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import k0

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))

from wdw_model import Packet,WDWModel
from ks_mass_dirac_observable_audit import evolve,mass0_apply,massD_apply


def log10_norm_sq(v):
    a=np.abs(v)
    m=float(a.max())
    if m==0:
        return float("-inf")
    return 2*math.log10(m)+math.log10(float(np.sum((a/m)**2)))


def forms_for_box(right,T,h=.2,dx=.04):
    model=WDWModel(h,left=-4.,right=right,dx=dx)
    packet=Packet(h,math.sqrt(h))
    chi0=packet.initial(model.x)
    chi=evolve(model,chi0,T)
    mchi=massD_apply(model,T,chi)

    r_m3=64*np.exp(3*(model.x+T))
    r_m32=8*np.exp(1.5*(model.x+T))

    a=r_m3*mchi
    b=massD_apply(model,T,r_m3*chi)
    c=r_m32*massD_apply(model,T,r_m32*chi)

    # log10(48 ||...||^2)
    add=math.log10(48)
    return {
      "right":right,"T":T,"dx":dx,"n":len(model.x),
      "log10_K_mass_then_radius":add+log10_norm_sq(a),
      "log10_K_radius_then_mass":add+log10_norm_sq(b),
      "log10_K_balanced":add+log10_norm_sq(c),
      "log10_bare_r_minus6":math.log10(48)+log10_norm_sq(r_m3*chi),
      "last_unit_probability":float(np.sum(abs(chi[model.x>right-1])**2)),
      "last_unit_mass_state_probability":float(
          np.sum(abs(mchi[model.x>right-1])**2)/max(np.sum(abs(mchi)**2),1e-300)
      )
    }


def mass_threshold_overlap(right,dx,h=.2):
    model=WDWModel(h,left=-4.,right=right,dx=dx)
    packet=Packet(h,math.sqrt(h))
    chi=packet.initial(model.x)
    mchi=mass0_apply(model,chi)
    kernel=k0(np.exp(-model.x)/h)
    # Euclidean amplitude v=sqrt(dx) psi => integral ~sqrt(dx) sum kernel*v.
    B=math.sqrt(dx)*np.sum(kernel*mchi)
    return {
      "right":right,"dx":dx,
      "real":float(B.real),"imag":float(B.imag),"abs":float(abs(B))
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()

    rows=[]
    for T in (0.,2.):
        for right in (12.,16.,24.,32.):
            rows.append(forms_for_box(right,T))

    overlap=[
      mass_threshold_overlap(16.,.04),
      mass_threshold_overlap(24.,.04),
      mass_threshold_overlap(24.,.02),
    ]
    ov=[x["abs"] for x in overlap]
    if min(ov) <= 1e-10:
        overlap_status="NOT_RESOLVED"
    else:
        overlap_status="NONZERO_NUMERICALLY_RESOLVED"

    evolved=[x for x in rows if x["T"]==2.]
    keys=[
      "log10_K_mass_then_radius",
      "log10_K_radius_then_mass",
      "log10_K_balanced",
    ]
    growth={}
    for key in keys:
        vals=[x[key] for x in evolved]
        growth[key]={
          "values":vals,
          "right":[x["right"] for x in evolved],
          "total_log10_growth":float(vals[-1]-vals[0]),
          "monotone_increasing":bool(all(b>a for a,b in zip(vals,vals[1:])))
        }

    # Existing continuum asymptotic says a nonzero threshold overlap produces
    # a polynomial x^-3 amplitude tail at fixed T>0. Any explicit positive
    # exponential radius factor then leaves the corresponding local quadratic
    # form outside L2. We only apply that logic to the declared forms here.
    all_grow=all(v["monotone_increasing"] for v in growth.values())
    status="FAIL" if overlap_status=="NONZERO_NUMERICALLY_RESOLVED" and all_grow else "PARTIAL"

    result={
      "schema":1,
      "status":status,
      "finding":"selected_mass_candidate_does_not_cure_explicit_radius_domain" if status=="FAIL" else "inconclusive",
      "mass_candidate":"M_D(T)=U(T) M0 U(T)^dagger from ks_mass_dirac_observable_audit.py",
      "classical_target":"K=48 mu^2/r^6",
      "positive_orderings":{
        "mass_then_radius":"48 ||r^-3 M_D psi||^2",
        "radius_then_mass":"48 ||M_D r^-3 psi||^2",
        "balanced":"48 ||r^-3/2 M_D r^-3/2 psi||^2"
      },
      "box_scan":rows,
      "evolved_T2_growth":growth,
      "mass_applied_threshold_overlap":overlap,
      "mass_applied_threshold_overlap_status":overlap_status,
      "domain_decision":{
        "mass_then_radius":"FAIL for the audited packet if the nonzero mass-applied threshold overlap persists in the continuum: r^-3 M_D psi is not L2 at fixed T>0.",
        "radius_then_mass":"FAIL already when r^-3 psi is outside the mass-operator domain; existing tail asymptotic has nonzero original threshold overlap.",
        "balanced":"FAIL already when r^-3/2 psi is outside the intermediate form domain for the same positive-exponential-tail reason.",
        "nonlocal_escape_hatch":"A different genuinely nonlocal relational curvature observable could have cancellations not represented by these explicit-radius positive forms; this audit does not rule that out."
      },
      "regulator_decision":"FAIL" if status=="FAIL" else "PARTIAL",
      "scope_note":"The failure is for the declared positive quadratic-form orderings built from the selected mass candidate. It is not a theorem that every possible quantum Kretschmann observable diverges."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
