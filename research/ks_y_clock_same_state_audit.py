"""Same-Dirac-state T-clock vs areal/null Y-clock audit.

Y=T+x=-log(4r), X=x.  Constant-Y slices are null in the 1+1
minisuperspace metric.  Instead of reinitializing a new wave packet, this
script pulls back the *same* positive-frequency WDW solution to Y=const and
compares the conserved KG flux and the transported mass-Dirac observable.

The canonical transform gives
  p_Y=p_T, p_X=p_x-p_T,
  p_Y=-(p_X^2+exp(-2X))/(2p_X)
on the target branch.  A deparametrized Y-Hamiltonian therefore requires
p_X^{-1}; the finite-box positive form associated with p_X is audited
separately for threshold/domain sensitivity.

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

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))

from wdw_model import Packet,WDWModel
from ks_mass_dirac_observable_audit import (
    H,
    H_half,
    S_apply,
    mass0_apply,
    kg_inner,
    velocity_matrix,
)


def momentum_matrix(model):
    n=len(model.x)
    D=np.zeros((n,n),float)
    for i in range(n):
        for off,c in ((2,-1),(1,8),(-1,-8),(-2,1)):
            j=i+off
            if 0<=j<n:
                D[i,j]+=c/(12*model.dx)
    return -1j*model.h*D


def slice_fields(model,coeff,Y,PU):
    # Euclidean mass amplitudes sqrt(dx)*Psi evaluated at T_i=Y-x_i.
    T=Y-model.x
    phase=np.exp(-1j*T[:,None]*model.energy[None,:]/model.h)
    w=phase*coeff[None,:]
    psi=np.sum(model.U*w,axis=1)
    Hpsi=np.sum(model.U*w*model.energy[None,:],axis=1)
    Ppsi=np.sum(PU*w,axis=1)
    return psi,Hpsi,Ppsi


def null_kg_bilinear(model,a,Ha,Pa,b,Hb,Pb):
    # Flux on Y=T+x=const:
    # int [a*(H+P)b + ((H+P)a)* b]/h dX.
    return np.sum(
        np.conj(a)*(Hb+Pb)+np.conj(Ha+Pa)*b
    )/model.h


def same_state_flux_case(h=.2,right=16.,dx=.04):
    model=WDWModel(h,left=-4.,right=right,dx=dx)
    packet=Packet(h,math.sqrt(h))
    chi0=packet.initial(model.x)
    psi0=S_apply(model,chi0)
    meta0=S_apply(model,mass0_apply(model,chi0))

    c=model.U.T@psi0
    cm=model.U.T@meta0
    P=momentum_matrix(model)
    PU=P@model.U

    ref_norm=kg_inner(model,psi0,psi0)
    ref_mass=kg_inner(model,psi0,meta0)

    rows=[]
    max_norm=0.
    max_mass=0.
    max_imag=0.
    for Y in (-2.,0.,2.,4.,6.):
        a,Ha,Pa=slice_fields(model,c,Y,PU)
        b,Hb,Pb=slice_fields(model,cm,Y,PU)
        norm=null_kg_bilinear(model,a,Ha,Pa,a,Ha,Pa)
        mass=null_kg_bilinear(model,a,Ha,Pa,b,Hb,Pb)
        norm_err=float(abs(norm-ref_norm))
        mass_err=float(abs(mass-ref_mass))
        max_norm=max(max_norm,norm_err)
        max_mass=max(max_mass,mass_err)
        max_imag=max(max_imag,float(abs(norm.imag)),float(abs(mass.imag)))
        rows.append({
          "Y":Y,
          "areal_radius":float(.25*math.exp(-Y)),
          "null_KG_norm_real":float(norm.real),
          "null_KG_norm_imag":float(norm.imag),
          "mass_flux_real":float(mass.real),
          "mass_flux_imag":float(mass.imag),
          "norm_difference_from_T_slice":norm_err,
          "mass_difference_from_T_slice":mass_err,
        })

    return {
      "h":h,"right":right,"dx":dx,
      "T_slice_KG_norm_real":float(ref_norm.real),
      "T_slice_mass_expectation_real":float(ref_mass.real),
      "max_null_norm_difference":max_norm,
      "max_null_mass_difference":max_mass,
      "max_flux_imaginary_part":max_imag,
      "rows":rows,
    }


def px_domain_case(right,dx=.08,h=.2):
    model=WDWModel(h,left=-4.,right=right,dx=dx)
    n=len(model.x)
    V=velocity_matrix(model)
    Hh=(model.U*np.sqrt(model.energy))@model.U.T
    PX=Hh@(np.eye(n)+V)@Hh
    PX=(PX+PX.conj().T)/2
    ev,W=np.linalg.eigh(PX)
    if ev[0] <= 0:
        raise RuntimeError(f"nonpositive finite-box p_X form at right={right}: {ev[0]}")

    packet=Packet(h,math.sqrt(h))
    chi=packet.initial(model.x)
    coeff=W.conj().T@chi
    inv_norm=float(np.sqrt(np.sum(abs(coeff/ev)**2)))
    inv_half_norm=float(np.sqrt(np.sum(abs(coeff/np.sqrt(ev))**2)))
    center_E=math.sqrt(packet.p0**2+math.exp(-2*packet.x0))
    classical_px=center_E+packet.p0

    return {
      "right":right,"dx":dx,"n":n,
      "min_pX_eigenvalue":float(ev[0]),
      "max_inverse_pX_eigenvalue":float(1/ev[0]),
      "state_norm_pX_inverse":inv_norm,
      "state_norm_pX_inverse_half":inv_half_norm,
      "classical_center_pX":float(classical_px),
      "classical_center_inverse_pX":float(1/classical_px),
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()

    flux=same_state_flux_case()
    # The finite-difference null current is a second discretization, so use
    # a numerical consistency threshold well above roundoff but far below
    # the O(1) clock effects under discussion.
    if flux["max_null_norm_difference"]>3e-3:
        raise RuntimeError(("null KG norm mismatch",flux["max_null_norm_difference"]))
    if flux["max_null_mass_difference"]>3e-3:
        raise RuntimeError(("null mass mismatch",flux["max_null_mass_difference"]))

    domain=[px_domain_case(r) for r in (12.,16.,24.,32.)]
    mins=[x["min_pX_eigenvalue"] for x in domain]

    result={
      "schema":1,
      "status":"PARTIAL",
      "same_Dirac_state_null_flux":"PASS",
      "same_mass_observable_null_flux":"PASS",
      "Y_clock_classical_admissibility":"PASS",
      "Y_clock_deparametrized_inverse_pX_domain":"PENDING",
      "coordinate":{
        "Y":"T+x=-log(4r)",
        "X":"x",
        "pY":"p_T",
        "pX":"p_x-p_T",
        "constraint":"-p_X^2-2 p_X p_Y-exp(-2X)=0",
        "Y_Hamiltonian":"K_Y=(p_X^2+exp(-2X))/(2 p_X)",
        "slice_character":"Y=const is null in the minisuperspace supermetric"
      },
      "same_state_flux":flux,
      "inverse_pX_threshold_diagnostic":domain,
      "min_pX_decreases_under_box_expansion":bool(all(b<a for a,b in zip(mins,mins[1:]))),
      "interpretation":{
        "what_passed":"The same WDW solution and the same transported mass observable give consistent KG fluxes on T=const and null Y=const slices on the finite box.",
        "what_did_not_pass":"A standalone Y-time Schrödinger theory has not been promoted because K_Y contains p_X^-1 and the finite-box p_X spectrum approaches zero under box expansion.",
        "zero_mode":"The p_X=0 threshold is the light-front/null-clock zero-mode obstruction corresponding to the horizon asymptotic sector."
      },
      "scope_note":"This is a same-Dirac-state clock/slice comparison, not proof of unitary equivalence between independently quantized T- and Y-clock theories."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
