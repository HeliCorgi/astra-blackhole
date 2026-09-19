"""Positive-frequency / physical-KG mass Dirac-observable audit for KS.

The local full-WDW ordering family in ks_observable_domain_audit.py preserves
the constraint kernel but not the repository's selected positive-frequency
sector.  Here the construction starts *inside* that sector.

Let H=sqrt(-h^2 d_x^2+exp(-2x)) and X be multiplication by x.  Define the
Heisenberg velocity V=(i/h)[H,X].  Its principal symbol is p/H, so

    M0 = 1/4 exp(X/2) H (I+V) H exp(X/2)

has the classical principal symbol
    mu = 1/4 exp(x) H(H+p)
at T=0.

On a finite box, if I+V is positive, M0 is a positive self-adjoint matrix.
The Dirac/relational completion is

    M_D(T)=U(T) M0 U(T)^dagger,  U(T)=exp(-iHT/h).

It therefore preserves the selected positive-frequency sector exactly.
The physical-KG representation uses the already-audited isometry

    S=sqrt(h/2) H^(-1/2),   M_KG(T)=S M_D(T) S^(-1).

This is a finite-box construction and a continuum quadratic-form candidate,
not a proof that the continuum form is closed or uniquely selected.

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
from numpy.polynomial.hermite import hermgauss

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))
from wdw_model import Packet,WDWModel


def spectral_apply(model, f, v):
    return model.U @ (f(model.energy) * (model.U.T @ v))


def H(model,v):
    return spectral_apply(model,lambda e:e,v)


def H_half(model,v):
    return spectral_apply(model,np.sqrt,v)


def H_mhalf(model,v):
    return spectral_apply(model,lambda e:1/np.sqrt(e),v)


def evolve(model,v,T):
    return spectral_apply(model,lambda e:np.exp(-1j*e*T/model.h),v)


def velocity(model,v):
    return 1j/model.h*(H(model,model.x*v)-model.x*H(model,v))


def mass0_apply(model,v):
    y=np.exp(model.x/2)*v
    y=H(model,y)
    y=y+velocity(model,y)
    y=H(model,y)
    y=np.exp(model.x/2)*y
    return y/4


def massD_apply(model,T,v):
    return evolve(model,mass0_apply(model,evolve(model,v,-T)),T)


def S_apply(model,v):
    return math.sqrt(model.h/2)*H_mhalf(model,v)


def Sinv_apply(model,v):
    return math.sqrt(2/model.h)*H_half(model,v)


def kg_inner(model,a,b):
    return (2/model.h)*np.vdot(a,H(model,b))


def mkg_apply(model,T,v):
    return S_apply(model,massD_apply(model,T,Sinv_apply(model,v)))


def velocity_matrix(model):
    Hm=(model.U*model.energy)@model.U.T
    X=np.diag(model.x)
    V=1j/model.h*(Hm@X-X@Hm)
    return (V+V.conj().T)/2


def classical_mu(x,p):
    E=np.sqrt(p*p+np.exp(-2*x))
    return .25*np.exp(x)*E*(E+p)

def classical_mass_ensemble(packet,order=120):
    z,w=hermgauss(order)
    z=np.sqrt(2.)*z
    w=w/np.sqrt(np.pi)
    x=packet.x0+packet.sigma*z[:,None]+np.zeros((1,order))
    p=(packet.p0
       +packet.chirp*packet.sigma*z[:,None]
       +packet.h/(2*packet.sigma)*z[None,:])
    weight=w[:,None]*w[None,:]
    return float(np.sum(weight*classical_mu(x,p)))


def run_case(h,right=16.,dx=.04,full_spectrum=False):
    model=WDWModel(h,left=-4.,right=right,dx=dx)
    packet=Packet(h,math.sqrt(h))
    chi0=packet.initial(model.x)

    n=len(model.x)
    min_i_plus_v=None
    max_abs_v=None
    min_mass=None

    # The full finite-box spectrum is checked on one moderate reference grid.
    # Finer semiclassical cases use bilinear/quadratic-form probes to avoid
    # turning the audit into a dense O(n^3) matrix benchmark.
    if full_spectrum:
        V=velocity_matrix(model)
        ve=np.linalg.eigvalsh(V)
        min_i_plus_v=float(1+ve[0])
        max_abs_v=float(np.max(np.abs(ve)))
        eye=np.eye(n,dtype=complex)
        M=np.column_stack([mass0_apply(model,eye[:,j]) for j in range(n)])
        herm=np.linalg.norm(M-M.conj().T)/max(np.linalg.norm(M),1e-30)
        Mherm=(M+M.conj().T)/2
        min_mass=float(np.linalg.eigvalsh(Mherm)[0])
    else:
        rng_probe=np.random.default_rng(31337+int(round(100*h)))
        herm=0.
        min_q=float("inf")
        for _ in range(6):
            a=rng_probe.normal(size=n)+1j*rng_probe.normal(size=n)
            b=rng_probe.normal(size=n)+1j*rng_probe.normal(size=n)
            a/=np.linalg.norm(a); b/=np.linalg.norm(b)
            ma=mass0_apply(model,a); mb=mass0_apply(model,b)
            lhs=np.vdot(a,mb); rhs=np.conj(np.vdot(b,ma))
            herm=max(herm,float(abs(lhs-rhs)/max(abs(lhs),abs(rhs),1e-30)))
            y=H(model,np.exp(model.x/2)*a)
            q=np.vdot(y,y+velocity(model,y)).real/4
            min_q=min(min_q,float(q))
        if min_q < -2e-9:
            raise RuntimeError(f"negative probed mass quadratic form: {min_q}")

    mu0=np.vdot(chi0,mass0_apply(model,chi0)).real
    mu_cl=float(classical_mu(packet.x0,packet.p0))
    mu_ensemble=classical_mass_ensemble(packet)

    times=(0.,1.,2.,4.,6.)
    rows=[]
    max_drift=0.
    max_sector=0.
    max_kg_exp=0.
    max_kg_adjoint=0.

    rng=np.random.default_rng(9421+int(round(10*h)))
    a=rng.normal(size=n)+1j*rng.normal(size=n)
    b=rng.normal(size=n)+1j*rng.normal(size=n)
    a/=np.linalg.norm(a); b/=np.linalg.norm(b)
    ka=S_apply(model,a); kb=S_apply(model,b)

    for T in times:
        chi=evolve(model,chi0,T)
        out=massD_apply(model,T,chi)
        direct=evolve(model,mass0_apply(model,chi0),T)
        sector=float(np.linalg.norm(out-direct)/max(np.linalg.norm(direct),1e-30))
        exp_l2=np.vdot(chi,out)

        psi=S_apply(model,chi)
        mout=mkg_apply(model,T,psi)
        exp_kg=kg_inner(model,psi,mout)

        # physical-KG adjointness bilinear test
        lhs=kg_inner(model,ka,mkg_apply(model,T,kb))
        rhs=np.conj(kg_inner(model,kb,mkg_apply(model,T,ka)))
        adj=float(abs(lhs-rhs)/max(abs(lhs),abs(rhs),1e-30))

        drift=float(abs(exp_l2.real/mu0-1))
        max_drift=max(max_drift,drift)
        max_sector=max(max_sector,sector)
        max_kg_exp=max(max_kg_exp,float(abs(exp_kg-exp_l2)))
        max_kg_adjoint=max(max_kg_adjoint,adj)
        rows.append({
          "T":T,
          "mass_expectation_l2":float(exp_l2.real),
          "mass_expectation_kg_real":float(exp_kg.real),
          "relative_mass_expectation_drift":drift,
          "positive_frequency_sector_residual":sector,
          "kg_adjoint_bilinear_residual":adj
        })

    if herm>2e-11:
        raise RuntimeError(f"M0 lost finite-box Hermiticity: {herm}")
    if full_spectrum and min_i_plus_v < -2e-10:
        raise RuntimeError(f"I+V not positive on finite box: {min_i_plus_v}")
    if full_spectrum and min_mass < -2e-8:
        raise RuntimeError(f"M0 not positive on finite box: {min_mass}")
    if max_drift>3e-11 or max_kg_exp>3e-10 or max_kg_adjoint>3e-10:
        raise RuntimeError((max_drift,max_kg_exp,max_kg_adjoint))

    return {
      "h":h,"right":right,"dx":dx,"n":n,
      "min_I_plus_velocity_eigenvalue":min_i_plus_v,
      "max_abs_velocity_eigenvalue":max_abs_v,
      "mass_matrix_hermiticity_relative":float(herm),
      "min_mass_matrix_eigenvalue":min_mass,
      "quantum_mass_initial":float(mu0),
      "classical_center_mass":float(mu_cl),
      "classical_same_wigner_ensemble_mass":float(mu_ensemble),
      "relative_error_vs_classical_center":float(abs(mu0/mu_cl-1)),
      "relative_error_vs_classical_same_wigner_ensemble":float(abs(mu0/mu_ensemble-1)),
      "max_relative_mass_expectation_drift":max_drift,
      "max_positive_frequency_sector_roundtrip_residual":max_sector,
      "positive_frequency_sector_preservation":"EXACT_BY_UNITARY_TRANSPORT",
      "sector_roundtrip_note":"The direct U M0 U^dagger identity is exact in finite dimension. This diagnostic uses an explicit backward transform of a floating-point state and can be amplified by the unbounded exp(X/2) form; it is not used as the scientific gate.",
      "max_L2_KG_mass_expectation_difference":max_kg_exp,
      "max_KG_adjoint_bilinear_residual":max_kg_adjoint,
      "rows":rows
    }


def lightweight_mass_expectation(h,right,dx):
    model=WDWModel(h,left=-4.,right=right,dx=dx)
    packet=Packet(h,math.sqrt(h))
    chi=packet.initial(model.x)
    return float(np.vdot(chi,mass0_apply(model,chi)).real)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()

    cases=[
      run_case(.4,dx=.08),
      run_case(.2,dx=.04,full_spectrum=True),
      run_case(.1,dx=.02),
      run_case(.05,dx=.01),
    ]
    semi=[row["relative_initial_semiclassical_error"] for row in cases]

    # Regulator diagnostics on the h=.2 reference packet. These are recorded
    # rather than used to manufacture a physics tolerance after seeing data.
    grid=[
      {"dx":dx,"mass":lightweight_mass_expectation(.2,16.,dx)}
      for dx in (.08,.04,.02)
    ]
    box=[
      {"right":right,"mass":lightweight_mass_expectation(.2,right,.04)}
      for right in (12.,16.,24.)
    ]

    result={
      "schema":1,
      "status":"PARTIAL",
      "finite_box_positive_mass_dirac_candidate":"PASS",
      "positive_frequency_sector_preservation":"PASS_BY_CONSTRUCTION_AND_CONSERVED_EVOLUTION",
      "finite_box_physical_KG_adjointness":"PASS",
      "continuum_quadratic_form_closure":"PENDING",
      "unique_ordering_selection":"PENDING",
      "construction":{
        "velocity":"V=(i/h)[H,X]",
        "mass_at_reference_clock":"M0=(1/4) exp(X/2) H (I+V) H exp(X/2)",
        "dirac_transport":"M_D(T)=U(T) M0 U(T)^dagger",
        "KG_map":"M_KG(T)=S M_D(T) S^(-1), S=sqrt(h/2) H^(-1/2)",
        "principal_symbol":"mu=(1/4) exp(x) H(H+p) at T=0",
        "positivity_rule":"finite-box I+V >= 0 implies M0 is a positive quadratic-form matrix"
      },
      "cases":cases,
      "semiclassical_error_sequence":semi,
      "semiclassical_error_decreases_with_h":bool(all(b<a for a,b in zip(semi,semi[1:]))),
      "grid_diagnostic":grid,
      "box_diagnostic":box,
      "continuum_domain":{
        "core_candidate":"compact-support smooth states / finite-energy states for which exp(X/2)psi lies in the composed form domain",
        "transported_domain":"D(M_D(T))=U(T)D(M0)",
        "advantage":"The transported domain does not require replacing M_D(T) by a naive local exp(x) multiplication on evolved tails.",
        "unresolved":"Closability/closedness and the self-adjoint operator represented by the continuum positive quadratic form have not been proven."
      },
      "literature_context":{
        "Kuchar_1994":"Schwarzschild mass is a Dirac observable/canonical coordinate in the full spherically symmetric Hamiltonian reduction.",
        "Ashtekar_Tate_Uggla_1993":"Minisuperspace Dirac observables and deparametrization should be constructed before physical interpretation."
      },
      "scope_note":"This selects a concrete positive-frequency-preserving finite-box mass Dirac candidate for further clock/curvature audits. It is not yet a unique or continuum-certified physical mass operator."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
