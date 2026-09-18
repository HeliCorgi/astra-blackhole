"""Consistent-histories audit helpers for the reduced quantum Bianchi IX model.

For the A-conditioned family the branch state is

    |psi_{A,j,k}> = P_k(s3) U(s3,s2) P_j(s2) U(s2,s1)
                    P_A(s1) U(s1,s0)|psi0>,

and D(alpha,beta)=<psi_beta|psi_alpha>.

A critical numerical requirement is that all second-wall branches see the SAME
approximate evolution operator.  Independent single-vector Lanczos runs are not
strictly linear because their Krylov spaces depend on the input vector.  The
history audit therefore propagates the three second-wall branches together in
one common block-Krylov subspace at every step.

Staged cubic regrids are also kept linear and UNNORMALIZED.  Renormalizing
branches separately would destroy their relative amplitudes.
"""
from __future__ import annotations

from dataclasses import dataclass
import itertools
import numpy as np
from scipy.interpolate import RectBivariateSpline

from quantum_bianchi_ix_model import QuantumBianchiIX, PacketSpec
from quantum_bianchi_ix_long_model import wall_masks, WALL_NAMES


@dataclass(frozen=True)
class HistoryTimes:
    s0: float = 2.0
    s1: float = 4.6
    s2: float = 10.385
    s3: float = 32.46


def standard_grids(cap=30.0):
    return (
        QuantumBianchiIX(hbar=.2,bp_bounds=(-5,5),bm_bounds=(-8,5),nx=100,ny=130,potential_cap=cap),
        QuantumBianchiIX(hbar=.2,bp_bounds=(-12,12),bm_bounds=(-18,10),nx=144,ny=168,potential_cap=cap),
        QuantumBianchiIX(hbar=.2,bp_bounds=(-18,18),bm_bounds=(-24,24),nx=144,ny=192,potential_cap=cap),
        QuantumBianchiIX(hbar=.2,bp_bounds=(-32,32),bm_bounds=(-36,36),nx=214,ny=240,potential_cap=cap),
    )


def project_sector(model, v, s, sector):
    sector=int(sector)
    if sector not in (0,1,2):
        raise ValueError("sector must be 0, 1, or 2")
    mask=wall_masks(model,float(s))[sector].ravel()
    return np.where(mask,np.asarray(v,complex),0.0)


def linear_regrid(v, old, new, order=3):
    """Long-audit cubic regrid without branch-dependent normalization."""
    psi=np.asarray(v,complex).reshape(old.nx,old.ny)/np.sqrt(old.dx*old.dy)
    sr=RectBivariateSpline(old.bp,old.bm,psi.real,kx=order,ky=order,s=0)
    si=RectBivariateSpline(old.bp,old.bm,psi.imag,kx=order,ky=order,s=0)
    val=(sr(new.bp,new.bm)+1j*si(new.bp,new.bm))*np.sqrt(new.dx*new.dy)
    return val.ravel()


def linear_regrid_matrix(V, old, new, order=3):
    V=np.asarray(V,complex)
    if V.ndim==1:
        return linear_regrid(V,old,new,order)[:,None]
    return np.column_stack([linear_regrid(V[:,i],old,new,order) for i in range(V.shape[1])])


def propagate_linear(model, v, s0, s1, ds, krylov_dim):
    """Single-vector propagation used only before the history split."""
    v=np.asarray(v,complex).copy()
    if s1 < s0:
        raise ValueError("increasing s required")
    if s1 == s0 or np.linalg.norm(v)==0:
        return v
    n=max(1,int(round((s1-s0)/ds)))
    h=(s1-s0)/n
    s=float(s0)
    for _ in range(n):
        v,_=model.step(v,s+h/2,h,krylov_dim)
        s+=h
    return v


def _initial_span(V, tol=1e-13):
    U,sv,_=np.linalg.svd(np.asarray(V,complex),full_matrices=False)
    if len(sv)==0 or sv[0]==0:
        return np.empty((V.shape[0],0),complex)
    rank=int(np.sum(sv > tol*sv[0]))
    return U[:,:rank]


def block_krylov_step(model, V, s_mid, ds, maxdim=54, tol=1e-12):
    """Apply one common projected exp(+i ds sqrt(A)/hbar) to all columns.

    The block-Krylov subspace is generated from the span of every branch column.
    Once the projected Hermitian matrix H=Q^*AQ is built, the SAME matrix
    function Q exp(+i ds sqrt(H)/hbar) Q^* acts on all columns.  This restores
    linear recombination inside the represented branch span, unlike independent
    per-branch Lanczos approximations.
    """
    V=np.asarray(V,complex)
    if V.ndim!=2:
        raise ValueError("V must be N x nbranch")
    if maxdim < V.shape[1]:
        raise ValueError("maxdim must cover starting branch span")
    if np.linalg.norm(V)==0:
        return V.copy(),{"dimension":0,"min_eigenvalue":0.0,"gram_step_drift":0.0}

    Q0=_initial_span(V)
    basis=[Q0[:,i].copy() for i in range(Q0.shape[1])]
    frontier=list(basis)
    while frontier and len(basis)<maxdim:
        new=[]
        for q in frontier:
            w=model.apply_A(q,s_mid)
            if basis:
                Q=np.column_stack(basis)
                for _ in range(2):
                    w-=Q@(Q.conj().T@w)
            n=float(np.linalg.norm(w))
            if n>tol:
                w/=n
                basis.append(w);new.append(w)
            if len(basis)>=maxdim:
                break
        frontier=new

    Q=np.column_stack(basis)
    AQ=np.column_stack([model.apply_A(Q[:,i],s_mid) for i in range(Q.shape[1])])
    H=Q.conj().T@AQ
    H=(H+H.conj().T)/2
    ev,U=np.linalg.eigh(H)
    if ev.min() < -1e-8:
        raise RuntimeError(f"negative block-Krylov Ritz value {ev.min()}")
    coords=Q.conj().T@V
    phase=np.exp(1j*ds*np.sqrt(np.maximum(ev,0))/model.hbar)
    newcoords=U@(phase[:,None]*(U.conj().T@coords))
    out=Q@newcoords
    g0=V.conj().T@V;g1=out.conj().T@out
    den=max(float(np.linalg.norm(g0)),1e-30)
    return out,{
        "dimension":int(Q.shape[1]),
        "min_eigenvalue":float(ev.min()),
        "gram_step_drift":float(np.linalg.norm(g1-g0)/den),
    }


def propagate_block(model, V, s0, s1, ds, maxdim=54):
    V=np.asarray(V,complex).copy()
    if s1 < s0:
        raise ValueError("increasing s required")
    if s1==s0:
        return V,{"max_dimension":0,"max_gram_step_drift":0.0,"min_ritz":0.0}
    n=max(1,int(round((s1-s0)/ds)))
    h=(s1-s0)/n;s=float(s0)
    stats=[];min_ritz=np.inf
    for _ in range(n):
        V,st=block_krylov_step(model,V,s+h/2,h,maxdim)
        stats.append(st);min_ritz=min(min_ritz,st["min_eigenvalue"]);s+=h
    return V,{
        "max_dimension":max(x["dimension"] for x in stats),
        "max_gram_step_drift":max(x["gram_step_drift"] for x in stats),
        "min_ritz":float(min_ritz),
        "actual_ds":float(h),
    }


def to_second_time(v_at_s1, grids, times=HistoryTimes()):
    g1,g2,_,_=grids
    v=propagate_linear(g1,v_at_s1,times.s1,6.0,.05,48)
    v=linear_regrid(v,g1,g2,3)
    v=propagate_linear(g2,v,6.0,times.s2,.1,52)
    return v


def propagate_second_branches(V2, grids, times=HistoryTimes(), maxdim=54):
    _,g2,g3,g4=grids
    stats=[]
    V,st=propagate_block(g2,V2,times.s2,12.0,.1,maxdim);stats.append(("g2",st))
    V=linear_regrid_matrix(V,g2,g3,3)
    V,st=propagate_block(g3,V,12.0,20.0,.15,maxdim);stats.append(("g3",st))
    V=linear_regrid_matrix(V,g3,g4,3)
    V,st=propagate_block(g4,V,20.0,times.s3,.2,maxdim);stats.append(("g4",st))
    return V,stats


def decoherence_matrix(branches):
    keys=list(branches)
    states=[np.asarray(branches[k],complex) for k in keys]
    n=len(keys);D=np.empty((n,n),complex)
    for a in range(n):
        for b in range(n):
            D[a,b]=np.vdot(states[b],states[a])
    return keys,D


def decoherence_measures(D):
    D=np.asarray(D,complex)
    diag=np.maximum(D.diagonal().real,0.0)
    off=D.copy();np.fill_diagonal(off,0)
    ratios=[];real_ratios=[]
    for i,j in itertools.permutations(range(len(diag)),2):
        den=np.sqrt(diag[i]*diag[j])
        if den>1e-15:
            ratios.append(abs(D[i,j])/den)
            real_ratios.append(abs(D[i,j].real)/den)
    return {
        "trace_diagonal":float(diag.sum()),
        "max_offdiagonal_abs":float(np.max(np.abs(off))),
        "max_offdiagonal_real_abs":float(np.max(np.abs(off.real))),
        "offdiagonal_l1":float(np.sum(np.abs(off))),
        "max_normalized_pair_coherence":float(max(ratios,default=0.0)),
        "max_normalized_real_pair":float(max(real_ratios,default=0.0)),
        "hermiticity_max_abs":float(np.max(np.abs(D-D.conj().T))),
        "minimum_eigenvalue_hermitian_part":float(np.linalg.eigvalsh((D+D.conj().T)/2).min()),
    }


def _build_branches(V3, g4, times):
    branches={}
    for j in range(3):
        for k in range(3):
            branches[(j,k)]=project_sector(g4,V3[:,j],times.s3,k)
    return branches


def _normalized_history_data(branches, first_a_weight):
    keys,D=decoherence_matrix(branches)
    coherent=sum(branches.values())
    norm2=float(np.vdot(coherent,coherent).real)
    Dn=D/norm2
    rows=[];diag=Dn.diagonal().real
    for i,(j,k) in enumerate(keys):
        rows.append({
            "history":f"A->{WALL_NAMES[j]}->{WALL_NAMES[k]}",
            "second_sector":WALL_NAMES[j],"third_sector":WALL_NAMES[k],
            "conditional_diagonal_weight":float(diag[i]),
            "joint_A_diagonal_weight":float(first_a_weight*diag[i]),
        })
    coarse=[]
    for k in range(3):
        ids=[i for i,(j,kk) in enumerate(keys) if kk==k]
        v=sum(branches[keys[i]] for i in ids)
        pcoarse=float(np.vdot(v,v).real/norm2);pfine=float(np.sum(diag[ids]))
        coarse.append({"third_sector":WALL_NAMES[k],"coarse_probability":pcoarse,
                       "sum_fine_diagonal_weights":pfine,"additivity_defect":pcoarse-pfine})
    return keys,Dn,rows,coarse,coherent,norm2


def run_a_conditioned_histories(main_maxdim=54,control_maxdim=66):
    times=HistoryTimes();grids=standard_grids();g1,g2,_,g4=grids
    packet=PacketSpec(hbar=.2,sigma=.3)
    psi_s1=propagate_linear(g1,g1.initial(packet),times.s0,times.s1,.05,48)
    a_branch=project_sector(g1,psi_s1,times.s1,0)
    first_a_weight=float(np.vdot(a_branch,a_branch).real)
    first_non_a_weight=float(np.vdot(psi_s1-a_branch,psi_s1-a_branch).real)

    a_s2=to_second_time(a_branch,grids,times)
    V2=np.column_stack([project_sector(g2,a_s2,times.s2,j) for j in range(3)])
    second_weights=np.real(np.diag(V2.conj().T@V2)).tolist()
    gram2=V2.conj().T@V2

    V3,stats=propagate_second_branches(V2,grids,times,main_maxdim)
    branches=_build_branches(V3,g4,times)
    keys,Dn,rows,coarse,coherent,norm2=_normalized_history_data(branches,first_a_weight)
    direct=V3.sum(axis=1)
    gram3=V3.conj().T@V3
    projection_closure=float(np.linalg.norm(coherent-direct)/max(np.linalg.norm(direct),1e-30))

    cV3,cstats=propagate_second_branches(V2,grids,times,control_maxdim)
    cbranches=_build_branches(cV3,g4,times)
    _,cDn,_,_,_,_=_normalized_history_data(cbranches,first_a_weight)

    off=lambda M:M-np.diag(np.diag(M))
    control={
        "main_maxdim":int(main_maxdim),"control_maxdim":int(control_maxdim),
        "max_matrix_abs_difference":float(np.max(np.abs(Dn-cDn))),
        "max_diagonal_abs_difference":float(np.max(np.abs(Dn.diagonal().real-cDn.diagonal().real))),
        "max_offdiagonal_abs_difference":float(np.max(np.abs(off(Dn)-off(cDn)))),
    }
    return {
        "times":times,"first_a_weight":first_a_weight,"first_non_a_weight":first_non_a_weight,
        "second_projected_raw_weights":second_weights,"keys":keys,"D_conditional":Dn,
        "rows":rows,"coarse":coarse,
        "closure":{
            "projection_recombination_relative_norm":projection_closure,
            "coherent_final_norm2":norm2,
            "sum_all_D_normalized":complex(np.sum(Dn)),
            "second_branch_gram_relative_change":float(np.linalg.norm(gram3-gram2)/max(np.linalg.norm(gram2),1e-30)),
            "main_block_stats":stats,
            "control_block_stats":cstats,
        },
        "measures":decoherence_measures(Dn),"control":control,
    }
