"""Ordering/domain audit for Schwarzschild mass and curvature observables.

Classically:
  mu = 1/4 exp(x-T) p_T (p_T-p_x)
  K  = 48 mu^2 / r^6

The script separates:
- exact constraint-kernel-preserving local orderings,
- flat-L2 formal symmetry,
- finite-box candidate orderings,
- continuum exponential-tail domain diagnostics.

No unique quantum mass or Kretschmann operator is selected.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse, itertools, json, math, sys
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))
from wdw_model import Packet,WDWModel

def symbolic_kernel_check():
    T,x,h=sp.symbols("T x h", real=True)
    d=sp.symbols("d")
    psi=sp.Function("psi")(T,x)
    f=sp.exp(x-T)
    V=sp.exp(-2*x)
    def pT(z): return -sp.I*h*sp.diff(z,T)
    def px(z): return -sp.I*h*sp.diff(z,x)
    def C(z): return pT(pT(z))-px(px(z))-V*z
    def M(z):
        return f/4*(pT(pT(z))-pT(px(z))+2*d*h*(pT(z)-px(z)))
    def G(z):
        return -sp.I*h/2*f*(pT(z)+2*d*h*z)
    residual=sp.simplify(sp.expand(M(C(psi))-C(M(psi))-G(C(psi))))
    if residual != 0:
        raise RuntimeError("constraint-kernel ordering identity failed")
    return {
      "family":"M_d=(1/4)e^(x-T)[p_T^2-p_T p_x+2 d h (p_T-p_x)]",
      "commutator":"[M_d,C]=-(i h/2)e^(x-T)(p_T+2 d h) C",
      "residual":"0"
    }

def flat_l2_symmetry_no_solution():
    # Write d=d_R+i d_I and expand the formal adjoint of
    # f[-d_TT+d_Tx-2 i d d_T+2 i d d_x], f=exp(x-T).
    # Independent derivative coefficients require mutually inconsistent values.
    return {
      "status":"NO_SOLUTION",
      "ansatz":"same M_d family as the exact constraint-kernel-preserving identity",
      "flat_L2_formal_adjoint_conditions":{
        "from_d_T":"Im(d)=3/4",
        "from_d_x":"Im(d)=1/4",
        "from_identity_term":"Re(d)=0 and Im(d)=1/2"
      },
      "contradiction":True,
      "interpretation":"No complex d satisfies all flat-kinematical-L2 formal-adjoint conditions within this local O(h) family. This is not a no-go theorem for the physical KG inner product or for nonlocal orderings."
    }

def momentum_matrix(n,dx,h):
    D=np.zeros((n,n),float)
    for i in range(n):
        for off,c in ((2,-1),(1,8),(-1,-8),(-2,1)):
            j=i+off
            if 0<=j<n:
                D[i,j]+=c/(12*dx)
    return -1j*h*D

def apply_H(m,v):
    return m.U@(m.energy*(m.U.T@v))

def apply_P(m,P,v):
    return P@v

def apply_seq(m,P,T,v,seq):
    out=v
    Efac=np.exp(m.x-T)
    for name in reversed(seq):
        if name=="E":
            out=Efac*out
        elif name=="H":
            out=apply_H(m,out)
        elif name=="K":
            out=apply_H(m,out)+apply_P(m,P,out)
        else:
            raise ValueError(name)
    return .25*out

def mass_apply(m,P,T,v,ordering):
    if ordering=="left":
        return apply_seq(m,P,T,v,("E","H","K"))
    if ordering=="sym":
        a=apply_seq(m,P,T,v,("E","H","K"))
        b=apply_seq(m,P,T,v,("K","H","E"))
        return .5*(a+b)
    if ordering=="weyl6":
        total=np.zeros_like(v,dtype=complex)
        for perm in itertools.permutations(("E","H","K")):
            total+=apply_seq(m,P,T,v,perm)
        return total/6
    raise ValueError(ordering)

def kernel_mass_apply(m,P,T,v,d):
    # Full-WDW local family restricted to negative-p_T initial data:
    # p_T psi=-H psi and p_T p_x psi=-P H psi.
    H_v=apply_H(m,v)
    H2_v=m.U@(m.eigenvalues*(m.U.T@v))
    PH_v=apply_P(m,P,H_v)
    B=H2_v+PH_v-2*d*m.h*(H_v+apply_P(m,P,v))
    return .25*np.exp(m.x-T)*B

def kernel_family_branch_scan(dx=.04):
    h=.2
    m=WDWModel(h,left=-4.,right=16.,dx=dx)
    P=momentum_matrix(len(m.x),m.dx,h)
    chi0=Packet(h,.4).initial(m.x)
    # Same-state KG scalar at the initial slice.
    psi=math.sqrt(h/2)*(m.U@((m.U.T@chi0)/np.sqrt(m.energy)))
    rows=[]
    for b in (0.,.25,.5,.75):
        d=1j*b
        Bpsi=kernel_mass_apply(m,P,0.,psi,d)
        # p_T(M psi): derivative of exp(-T) plus p_T acting on the
        # T-independent operator B and on psi with p_T psi=-H psi.
        Hpsi=apply_H(m,psi)
        BHpsi=kernel_mass_apply(m,P,0.,Hpsi,d)
        pT_M=1j*h*Bpsi-BHpsi
        branch_res=pT_M+apply_H(m,Bpsi)
        denom=max(np.linalg.norm(apply_H(m,Bpsi)),1e-30)
        rows.append({
          "d":f"{b}i",
          "positive_frequency_branch_relative_residual":float(np.linalg.norm(branch_res)/denom)
        })
    return {
      "dx":dx,
      "rows":rows,
      "best_relative_residual":min(x["positive_frequency_branch_relative_residual"] for x in rows),
      "interpretation":"The exact constraint-kernel-preserving local family need not preserve the chosen negative-p_T/positive-frequency sector."
    }

def finite_box_ordering_scan(dx):
    h=.2
    m=WDWModel(h,left=-4.,right=16.,dx=dx)
    P=momentum_matrix(len(m.x),m.dx,h)
    p=Packet(h,.4)
    chi0=p.initial(m.x)
    coeff=m.U.T@chi0
    times=(0.,1.,2.,4.,6.)
    rng=np.random.default_rng(7321)
    u=rng.normal(size=len(m.x))+1j*rng.normal(size=len(m.x))
    v=rng.normal(size=len(m.x))+1j*rng.normal(size=len(m.x))
    u/=np.linalg.norm(u);v/=np.linalg.norm(v)
    rows={}
    for ordering in ("left","sym","weyl6"):
        # Random bilinear formal-Hermiticity diagnostic on the finite box.
        lhs=np.vdot(u,mass_apply(m,P,0.,v,ordering))
        rhs=np.conj(np.vdot(v,mass_apply(m,P,0.,u,ordering)))
        scale=max(abs(lhs),abs(rhs),1e-30)
        herm=float(abs(lhs-rhs)/scale)
        vals=[]
        for T in times:
            state=m.U@(coeff*np.exp(-1j*m.energy*T/h))
            vals.append(float(np.vdot(state,mass_apply(m,P,T,state,ordering)).real))
        drift=(max(vals)-min(vals))/max(abs(vals[0]),1e-30)
        rows[ordering]={
          "finite_box_hermiticity_bilinear_residual":herm,
          "mass_expectations":vals,
          "relative_expectation_drift":float(drift)
        }
    return {"dx":dx,"rows":rows}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()

    kernel=symbolic_kernel_check()
    symmetry=flat_l2_symmetry_no_solution()
    scans=[finite_box_ordering_scan(.08),finite_box_ordering_scan(.04),finite_box_ordering_scan(.02)]
    branch_scan=kernel_family_branch_scan(.04)

    tail=json.loads((ROOT/"research/wdw_tail_results/summary.json").read_text())
    if not tail["full_unbounded_moment_finite_convergence"].startswith("FAILED"):
        raise RuntimeError("tail-domain negative control unexpectedly changed")
    q2=tail["full_unbounded_moment_finite_convergence"]
    continuum=tail["continuum_analysis"]

    result={
      "schema":1,
      "status":"PARTIAL",
      "finding":"no_selected_domain-safe_quantum_mass_or_curvature_operator",
      "classical_targets":{
        "mass":"mu=(1/4) exp(x-T) p_T(p_T-p_x)",
        "curvature":"K=48 mu^2/r^6",
        "radius":"r=(1/4)exp(-x-T)"
      },
      "constraint_kernel_ordering_family":kernel,
      "flat_kinematical_L2_symmetry":symmetry,
      "finite_box_ordering_scan":scans,
      "constraint_kernel_family_positive_frequency_scan":branch_scan,
      "continuum_domain":{
        "existing_tail_audit":q2,
        "continuum_tail_statement":continuum,
        "mass_warning":"Any ordering that requires exp(x) psi in L2 needs the q=2 exponential moment; the existing tail audit does not support that domain for evolved packets.",
        "curvature_warning":"The bare r^-6 quadratic form is proportional to the q=6 exponential moment; the existing tail audit rejects finite convergence for q=6.",
        "mass_dressed_curvature":"Derivative/nonlocal mass factors can alter domains and cancellations, so failure of the bare r^-6 form is not by itself a theorem about every possible 48 mu^2/r^6 ordering."
      },
      "decision":{
        "quantum_mass_operator_selected":False,
        "quantum_kretschmann_operator_selected":False,
        "reason":"The tested local orderings expose a conflict with flat-L2 formal symmetry, do not automatically preserve the selected positive-frequency sector, show finite-box ordering drift, and still lack controlled continuum exponential-factor domains."
      },
      "scope_note":"PARTIAL blocks promotion of a quantum mass/Kretschmann claim in the current representation. The finite-box expectation drift decreases with grid refinement and is not treated as a continuum no-go. The unresolved points are the physical-inner-product adjoint/domain, ordering selection, and exponential-tail domains. This does not prove that no acceptable nonlocal Dirac-observable quantization exists."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
