"""Affine-mass selection audit for the KS positive-frequency sector.

This audit separates a classical affine identity from a quantum selection-rule
obstruction. In the reduced classical KS branch,

    E=sqrt(p^2+exp(-2x)),   u=asinh(p exp(x)),

obey {u,E}=1 and the T=0 Schwarzschild mass is

    mu0 = (1/4) exp(x) E(E+p) = (E/4) exp(u),

so {E,mu0}=-mu0. It is tempting to demand exact quantum covariance

    exp(-iHt/h) M exp(+iHt/h) = exp(-t) M.

If H is self-adjoint and semibounded while M is a nonzero positive
self-adjoint operator, this exact covariance is incompatible with the selected
positive-frequency Hilbert space: on the support of M, Q=log M would satisfy a
Weyl translation relation with H, forcing spec(H) to be invariant under all
real translations. That contradicts spec(H) subset [0,infinity).

This is a no-go for using exact affine covariance as the ordering-selection
principle inside the positive-frequency sector. It is not a proof that the
repository's positive quadratic-form candidate is non-closable.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))
from wdw_model import WDWModel
from ks_mass_dirac_observable_audit import mass0_apply


def classical_affine_symbolic() -> dict:
    x,p=sp.symbols("x p", real=True)
    E=sp.sqrt(p*p+sp.exp(-2*x))
    u=sp.asinh(p*sp.exp(x))
    mu=sp.exp(x)*E*(E+p)/4
    pb=lambda f,g: sp.simplify(sp.diff(f,x)*sp.diff(g,p)-sp.diff(f,p)*sp.diff(g,x))
    return {
        "poisson_u_E": str(sp.simplify(pb(u,E))),
        "poisson_E_mu_plus_mu": str(sp.simplify(pb(E,mu)+mu)),
        "poisson_mu_E_minus_mu": str(sp.simplify(pb(mu,E)-mu)),
    }


def null_canonical_symbolic() -> dict:
    # U,W are the null minisuperspace coordinates. W is used here instead of V
    # to avoid collision with the quantum velocity V=(i/h)[H,X].
    U,W,pU,pW=sp.symbols("U W p_U p_W", nonzero=True, real=True)
    pT=-U*pU+W*pW
    px=-U*pU-W*pW
    ckin=sp.expand(pT*pT-px*px)
    # exp(-2x)=-16 U W, so C=-4 U W (pU pW-4).
    C=sp.factor(ckin+16*U*W)
    mass=sp.simplify(sp.Rational(1,8)*pT*pW)
    mass_shell=sp.simplify(mass.subs(pU,4/pW))
    return {
        "pT2_minus_px2": str(sp.factor(ckin)),
        "constraint_factor": str(C),
        "mass_pT_pW_over_8": str(mass),
        "mass_on_shell": str(mass_shell),
    }


def affine_weyl_no_go() -> dict:
    return {
        "status":"PASS",
        "assumptions":[
            "H is self-adjoint and semibounded (the selected positive-frequency sector has spectrum [0,infinity)).",
            "M is a nonzero positive self-adjoint operator; restrict to its support so log(M) is defined by spectral calculus.",
            "Exact covariance U(t) M U(t)^*=exp(-t) M holds for all real t, U(t)=exp(-iHt/h).",
        ],
        "steps":[
            "Q=log(M) gives U(t) Q U(t)^*=Q-t I.",
            "W(s)=exp(i s Q) therefore obeys U(t)W(s)=exp(-i s t)W(s)U(t).",
            "Hence W(s)^* H W(s)=H+h s I for every real s.",
            "Unitary conjugation preserves spectrum, so spec(H)=spec(H)+h s for every real s.",
            "A nonempty semibounded spectrum cannot be invariant under all real translations.",
        ],
        "conclusion":"No nonzero positive self-adjoint M on the selected semibounded positive-frequency Hilbert space can obey the full classical affine covariance exactly for all real t.",
        "scope":"This blocks exact affine covariance as an ordering-uniqueness principle. It does not decide closability of the declared q_M form, which intentionally uses relational transport rather than an exact affine commutator law.",
    }


def finite_box_affine_residual(*,h=.2,left=-4.,right=12.,dx=.08) -> dict:
    model=WDWModel(h,left=left,right=right,dx=dx)
    n=len(model.x)
    Hm=(model.U*model.energy)@model.U.T
    eye=np.eye(n,dtype=complex)
    M=np.column_stack([mass0_apply(model,eye[:,j]) for j in range(n)])
    comm=Hm@M-M@Hm
    target=-1j*h*M
    residual=float(np.linalg.norm(comm-target)/max(np.linalg.norm(target),1e-30))
    antiherm=float(np.linalg.norm(comm+comm.conj().T)/max(np.linalg.norm(comm),1e-30))
    return {
        "h":h,"left":left,"right":right,"dx":dx,"n":n,
        "relative_residual_for_[H,M0]=-ihM0":residual,
        "commutator_antihermiticity_residual":antiherm,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    cl=classical_affine_symbolic()
    null=null_canonical_symbolic()
    nogo=affine_weyl_no_go()
    if cl != {"poisson_u_E":"1","poisson_E_mu_plus_mu":"0","poisson_mu_E_minus_mu":"0"}:
        raise RuntimeError(cl)
    if null["constraint_factor"] not in ("-4*U*W*(p_U*p_W - 4)", "-4*U*W*(p_U*p_W-4)"):
        raise RuntimeError(null)
    rows=[finite_box_affine_residual(right=r) for r in (8.,12.,16.)]
    if not all(r["commutator_antihermiticity_residual"]<2e-12 for r in rows):
        raise RuntimeError(rows)
    result={
      "schema":1,
      "status":"PARTIAL",
      "model":"Vacuum Schwarzschild/Kantowski-Sachs positive-frequency WDW minisuperspace",
      "units":"dimensionless repository conventions",
      "retained_variables":["x","p_x","T positive-frequency branch"],
      "omitted_variables":["matter","inhomogeneous modes","environment"],
      "classical_affine_structure":{
        "status":"PASS","symbolic":cl,
        "canonical_variables":"E=sqrt(p_x^2+exp(-2x)), u=asinh(p_x exp(x)); {u,E}=1",
        "reference_mass":"mu0=(1/4)exp(x)E(E+p_x)=(E/4)exp(u)",
        "affine_bracket":"{E,mu0}=-mu0",
      },
      "null_coordinate_bridge":{
        "status":"PASS","symbolic":null,
        "coordinates":"U=-exp(-x-T)/4, W=exp(-x+T)/4",
        "momenta":"p_T=-U p_U+W p_W, p_x=-U p_U-W p_W",
        "constraint":"C=-4 U W (p_U p_W-4)",
        "mass":"mu=(1/8)p_T p_W; on shell mu=(W p_W^2-4U)/8",
      },
      "quantum_affine_covariance_no_go":nogo,
      "declared_mass_candidate_negative_control":{
        "status":"NOT_EXACT_AFFINE_COVARIANT_ON_FINITE_BOX",
        "relation_tested":"[H,M0]=-i h M0",
        "finite_box_controls":rows,
        "interpretation":"The O(1) residual is expected to remain nonzero because exact affine covariance is incompatible with the semibounded positive-frequency Hilbert-space assumptions used by the no-go. This is not a closability failure.",
      },
      "literature_comparison":{
        "Cavaglia_deAlfaro_Filippov_1995":"A different Schwarzschild minisuperspace/gauge-fixed quantization finds an affine invariant algebra and a Hermitian mass operator J that has no self-adjoint extension because its conjugate variable has half-line support, while J^2 admits self-adjoint realizations. This is a comparison/warning, not an identification with repository M0.",
      },
      "mass_continuum_form":{
        "status":"PENDING","unchanged":True,
        "reason":"The no-go excludes exact affine covariance as a selection principle but neither proves nor disproves closability of q_M=(1/4)<H exp(X/2)psi,(I+V)H exp(X/2)psi>.",
      },
      "ordering_selection":{
        "status":"FAIL_FOR_EXACT_AFFINE_COVARIANCE_CRITERION",
        "finding":"Exact quantum preservation of the classical affine scaling law cannot be imposed together with semibounded H and nonzero positive self-adjoint mass in the selected sector.",
      },
      "limitations":[
        "The Weyl no-go assumes exact covariance for all real clock shifts and self-adjointness/positivity of M; weaker semiclassical covariance is not excluded.",
        "The Cavaglia et al. operator uses a different canonical representation and gauge-fixed Hilbert space; its deficiency indices are not transferred to M0.",
        "Finite-box commutator residuals are negative controls only and are not continuum proofs.",
        "No curvature or singularity claim follows from this audit.",
      ],
      "reproducible_command":"python research/ks_affine_mass_obstruction_audit.py --out artifacts/ks-affine-mass/summary.json",
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
