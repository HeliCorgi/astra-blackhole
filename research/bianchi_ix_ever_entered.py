"""Finite-clock "ever entered B" class-operator audit for reduced quantum Bianchi IX.

Starting from the hard A-sector conditioned state at s1=4.6, define the
B-region as the union of the asymptotic B+ and B- wall sectors.  The no-entry
branch is approximated with a weak complex absorbing potential,

    G_eff(s) = G(s) - i V0 F_B,

where G=-sqrt(A) is the existing reduced generator and F_B is a smooth
approximation to the characteristic function of the B-region.

For the chosen finite clock interval, the complementary branch is defined by

    |psi_enter> = U |psi_A> - |psi_noB>.

This is a finite-clock complex-potential analogue of an "ever entered region"
class operator.  It is NOT the timeless Hamiltonian-constraint-commuting
S-matrix class operator of Halliwell.

The unrestricted and absorbed states share one block-Krylov approximation to
the unitary Hamiltonian substep at each time step.  Cubic grid enlargements are
linear and unnormalized for both columns.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from quantum_bianchi_ix_model import PacketSpec, initial_wigner_ensemble, _rhs_vec
from bianchi_ix_histories import (
    HistoryTimes, standard_grids, project_sector, propagate_linear,
    block_krylov_step, linear_regrid_matrix,
)


SQ3=np.sqrt(3.0)


@dataclass(frozen=True)
class CAPSpec:
    v0: float = 0.05
    width: float = 0.30
    maxdim: int = 66


def signed_B_distance(bp,bm):
    """Signed Euclidean distance from the A/B asymptotic-sector boundary.

    B={B+,B-} dominates A iff 3 beta_+ + sqrt(3)|beta_-| > 0.
    Positive distance points into B.
    """
    return (3.0*np.asarray(bp)+SQ3*np.abs(np.asarray(bm)))/np.sqrt(12.0)


def absorber_profile(model,width):
    """Compact C1 transition from 0 in A to 1 in B across width."""
    width=float(width)
    d=signed_B_distance(model.BP,model.BM)
    if width<0:
        raise ValueError("width must be nonnegative")
    if width==0:
        return (d>0).astype(float)
    x=np.clip(0.5+d/width,0.0,1.0)
    return x*x*(3.0-2.0*x)


def binary_sector(model):
    return signed_B_distance(model.BP,model.BM)>0


def initial_A_conditioned(grids=None,times=HistoryTimes()):
    if grids is None: grids=standard_grids()
    g1=grids[0]
    packet=PacketSpec(hbar=.2,sigma=.3)
    psi=propagate_linear(g1,g1.initial(packet),times.s0,times.s1,.05,48)
    a=project_sector(g1,psi,times.s1,0)
    w=float(np.vdot(a,a).real)
    if w<=0: raise RuntimeError("zero A conditioning weight")
    return a/np.sqrt(w),w,psi


def pair_step(model,full,noB,s_mid,ds,spec:CAPSpec):
    """One Strang CAP step with a common unitary block-Krylov substep."""
    prof=absorber_profile(model,spec.width).ravel()
    damp=np.exp(-0.5*spec.v0*ds*prof/model.hbar)
    pre=np.column_stack([np.asarray(full,complex),np.asarray(noB,complex)*damp])
    out,st=block_krylov_step(model,pre,float(s_mid),float(ds),int(spec.maxdim))
    return out[:,0],out[:,1]*damp,st


def propagate_pair_segment(model,full,noB,s0,s1,ds,spec:CAPSpec,save_times=()):
    if s1<s0: raise ValueError("increasing s required")
    save=sorted(float(x) for x in save_times if s0-1e-12<=x<=s1+1e-12)
    rows={};idx=0
    if save and abs(save[0]-s0)<1e-10:
        rows[save[0]]={"full_norm2":float(np.vdot(full,full).real),"noB_norm2":float(np.vdot(noB,noB).real)}
        idx=1
    n=max(1,int(round((s1-s0)/ds)));h=(s1-s0)/n;s=float(s0)
    max_gram=0.0;maxdim=0;min_ritz=np.inf
    for _ in range(n):
        full,noB,st=pair_step(model,full,noB,s+h/2,h,spec)
        max_gram=max(max_gram,float(st["gram_step_drift"]));maxdim=max(maxdim,int(st["dimension"]))
        min_ritz=min(min_ritz,float(st["min_eigenvalue"]));s+=h
        while idx<len(save) and s>=save[idx]-h/2:
            t=save[idx]
            rows[t]={"full_norm2":float(np.vdot(full,full).real),"noB_norm2":float(np.vdot(noB,noB).real)}
            idx+=1
    return full,noB,{"actual_ds":float(h),"max_unitary_gram_drift":max_gram,"max_dimension":maxdim,"min_ritz":float(min_ritz),"rows":rows}


def _regrid_pair(V,old,new,s):
    before=[float(np.vdot(V[:,i],V[:,i]).real) for i in range(2)]
    W=linear_regrid_matrix(V,old,new,3)
    after=[float(np.vdot(W[:,i],W[:,i]).real) for i in range(2)]
    return W,{"s":float(s),"full_norm_ratio":after[0]/before[0] if before[0] else None,
              "noB_norm_ratio":after[1]/before[1] if before[1] else None}


def run_quantum(spec=CAPSpec(),times=HistoryTimes()):
    if spec.v0<0: raise ValueError("v0 must be nonnegative")
    grids=standard_grids();g1,g2,g3,g4=grids
    init,a_weight,_=initial_A_conditioned(grids,times)
    full=init.copy();noB=init.copy();segments=[];regrids=[];curve={times.s1:{"full_norm2":1.0,"noB_norm2":1.0}}
    schedule=[
        (g1,times.s1,6.0,.05,[6.0]),
        (g2,6.0,12.0,.1,[10.385,12.0]),
        (g3,12.0,20.0,.15,[14.0,20.0]),
        (g4,20.0,times.s3,.2,[24.0,30.0,times.s3]),
    ]
    models=[g1,g2,g3,g4]
    for q,(model,s0,s1,ds,saves) in enumerate(schedule):
        full,noB,st=propagate_pair_segment(model,full,noB,s0,s1,ds,spec,saves)
        segments.append({"interval":[s0,s1],**{k:v for k,v in st.items() if k!="rows"}})
        curve.update(st["rows"])
        if q<3:
            V=np.column_stack([full,noB]);V,rg=_regrid_pair(V,models[q],models[q+1],s1)
            full,noB=V[:,0],V[:,1];regrids.append(rg)
            curve[s1]={"full_norm2":float(np.vdot(full,full).real),"noB_norm2":float(np.vdot(noB,noB).real)}

    entered=full-noB
    states=[noB,entered]
    norm=float(np.vdot(full,full).real)
    D=np.array([[np.vdot(states[j],states[i]) for j in range(2)] for i in range(2)],complex)/norm
    diag=D.diagonal().real
    off=D[0,1]
    den=np.sqrt(max(diag[0],0)*max(diag[1],0))
    out={
      "spec":{"v0":spec.v0,"width":spec.width,"maxdim":spec.maxdim},
      "conditioning":{"s":times.s1,"A_weight":a_weight},
      "interval":[times.s1,times.s3],
      "D":[[{"re":float(z.real),"im":float(z.imag)} for z in row] for row in D],
      "noB_diagonal_weight":float(diag[0]),
      "everB_diagonal_weight":float(diag[1]),
      "offdiagonal_abs":float(abs(off)),
      "offdiagonal_real":float(off.real),
      "normalized_coherence":float(abs(off)/den) if den>0 else 0.0,
      "normalized_real_coherence":float(abs(off.real)/den) if den>0 else 0.0,
      "diagonal_sum":float(diag.sum()),
      "additivity_defect":float(1.0-diag.sum()),
      "sum_all_D":{"re":float(D.sum().real),"im":float(D.sum().imag)},
      "closure_relative_norm":float(np.linalg.norm(noB+entered-full)/np.linalg.norm(full)),
      "unrestricted_final_norm2":norm,
      "segments":segments,"regrids":regrids,
      "survival_curve":[{"s":float(s),**curve[s]} for s in sorted(curve)],
      "note":"Diagonal weights become ordinary probabilities only if the two-history family decoheres. CAP parameters are regulator choices and must show a stable window."
    }
    return out


def classical_ever_B(n=4096,seed=20260917,ds=.01,times=HistoryTimes()):
    packet=PacketSpec(hbar=.2,sigma=.3);Y=initial_wigner_ensemble(packet,n,seed);s=2.0
    def advance(target):
        nonlocal Y,s
        while s<target-1e-12:
            h=min(float(ds),float(target-s))
            k1=_rhs_vec(s,Y);k2=_rhs_vec(s+h/2,Y+h*k1/2);k3=_rhs_vec(s+h/2,Y+h*k2/2);k4=_rhs_vec(s+h,Y+h*k3)
            Y += h*(k1+2*k2+2*k3+k4)/6;s+=h
    advance(times.s1)
    A=signed_B_distance(Y[0],Y[1])<=0
    count=int(A.sum())
    ever=np.zeros(n,bool);first=np.full(n,np.nan)
    while s<times.s3-1e-12:
        h=min(float(ds),float(times.s3-s))
        k1=_rhs_vec(s,Y);k2=_rhs_vec(s+h/2,Y+h*k1/2);k3=_rhs_vec(s+h/2,Y+h*k2/2);k4=_rhs_vec(s+h,Y+h*k3)
        Y += h*(k1+2*k2+2*k3+k4)/6;s+=h
        hit=A & (~ever) & (signed_B_distance(Y[0],Y[1])>0)
        first[hit]=s;ever |= hit
    frac=float(np.mean(ever[A])) if count else float("nan")
    q=np.nanquantile(first[A],[.1,.5,.9]).tolist() if np.any(np.isfinite(first[A])) else [None,None,None]
    return {"n":int(n),"seed":int(seed),"ds":float(ds),"first_A_count":count,"first_A_frequency":count/n,
            "everB_given_A":frac,"first_entry_s_quantiles_10_50_90":q}


def run_classical_control():
    a=classical_ever_B(ds=.01);b=classical_ever_B(ds=.005)
    return {"main":a,"refined":b,"everB_abs_difference":abs(a["everB_given_A"]-b["everB_given_A"])}
