"""Consistent-histories audit helpers for the reduced quantum Bianchi IX model.

The finite-time class branch state for an A-conditioned three-time history is

    |psi_{A,j,k}> = P_k(s3) U(s3,s2) P_j(s2) U(s2,s1) P_A(s1) U(s1,s0)|psi0>.

The decoherence functional is D(alpha,beta)=<psi_beta|psi_alpha>.
Projectors are the same A/B+/B- configuration-space wall sectors used by the
long-horizon audit.  Staged regridding is kept linear and UNNORMALIZED here:
renormalizing each branch would change history weights and is therefore
forbidden in this module.

This is a finite-time history diagnostic for one chosen reduced quantization,
not a unique Wheeler-DeWitt history construction or an observation.
"""
from __future__ import annotations

from dataclasses import dataclass
import itertools
import numpy as np
from scipy.interpolate import RectBivariateSpline

from quantum_bianchi_ix_model import QuantumBianchiIX, PacketSpec
from quantum_bianchi_ix_long_model import wall_masks, WALL_NAMES, fidelity


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
    """The long-audit cubic regrid without branch-dependent normalization."""
    psi=np.asarray(v,complex).reshape(old.nx,old.ny)/np.sqrt(old.dx*old.dy)
    sr=RectBivariateSpline(old.bp,old.bm,psi.real,kx=order,ky=order,s=0)
    si=RectBivariateSpline(old.bp,old.bm,psi.imag,kx=order,ky=order,s=0)
    val=(sr(new.bp,new.bm)+1j*si(new.bp,new.bm))*np.sqrt(new.dx*new.dy)
    return val.ravel()


def propagate_linear(model, v, s0, s1, ds, krylov_dim):
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


def to_second_time(v_at_s1, grids, times=HistoryTimes()):
    g1,g2,_,_=grids
    v=propagate_linear(g1,v_at_s1,times.s1,6.0,.05,48)
    v=linear_regrid(v,g1,g2,3)
    v=propagate_linear(g2,v,6.0,times.s2,.1,52)
    return v


def second_to_twenty(v_at_s2, grids, times=HistoryTimes()):
    _,g2,g3,g4=grids
    v=propagate_linear(g2,v_at_s2,times.s2,12.0,.1,52)
    v=linear_regrid(v,g2,g3,3)
    v=propagate_linear(g3,v,12.0,20.0,.15,52)
    v=linear_regrid(v,g3,g4,3)
    return v


def twenty_to_third(v_at_20, g4, times=HistoryTimes(), krylov_dim=52, ds=.2):
    return propagate_linear(g4,v_at_20,20.0,times.s3,ds,krylov_dim)


def decoherence_matrix(branches):
    """Rows/columns follow the input key order; D[a,b]=<psi_b|psi_a>."""
    keys=list(branches)
    states=[np.asarray(branches[k],complex) for k in keys]
    n=len(keys)
    D=np.empty((n,n),complex)
    for a in range(n):
        for b in range(n):
            D[a,b]=np.vdot(states[b],states[a])
    return keys,D


def decoherence_measures(D):
    D=np.asarray(D,complex)
    if D.ndim!=2 or D.shape[0]!=D.shape[1]:
        raise ValueError("square matrix required")
    diag=np.maximum(D.diagonal().real,0.0)
    off=D.copy()
    np.fill_diagonal(off,0)
    max_abs=float(np.max(np.abs(off))) if D.size else 0.0
    ratios=[]
    for i,j in itertools.permutations(range(len(diag)),2):
        den=np.sqrt(diag[i]*diag[j])
        if den>1e-15:
            ratios.append(abs(D[i,j])/den)
    return {
        "trace_diagonal":float(diag.sum()),
        "max_offdiagonal_abs":max_abs,
        "offdiagonal_l1":float(np.sum(np.abs(off))),
        "max_normalized_pair_coherence":float(max(ratios,default=0.0)),
        "hermiticity_max_abs":float(np.max(np.abs(D-D.conj().T))),
        "minimum_eigenvalue_hermitian_part":float(np.linalg.eigvalsh((D+D.conj().T)/2).min()),
    }


def run_a_conditioned_histories(control_krylov=None):
    """Compute the nine A -> j -> k histories through the third wall.

    Returns raw branch states plus closure/numerical metadata.  The first A
    projector is NOT renormalized.  A normalized conditional decoherence
    functional is formed only after all branch amplitudes have been computed.
    """
    times=HistoryTimes()
    grids=standard_grids()
    g1,g2,_,g4=grids
    packet=PacketSpec(hbar=.2,sigma=.3)

    psi=g1.initial(packet)
    psi_s1=propagate_linear(g1,psi,times.s0,times.s1,.05,48)
    a_branch=project_sector(g1,psi_s1,times.s1,0)
    first_a_weight=float(np.vdot(a_branch,a_branch).real)
    first_non_a_weight=float(np.vdot(psi_s1-a_branch,psi_s1-a_branch).real)

    a_s2=to_second_time(a_branch,grids,times)
    second_weights=[]
    at20={}
    for j in range(3):
        vj=project_sector(g2,a_s2,times.s2,j)
        second_weights.append(float(np.vdot(vj,vj).real))
        at20[j]=second_to_twenty(vj,grids,times)

    branches={}
    prefinal={}
    for j in range(3):
        v=twenty_to_third(at20[j],g4,times,52,.2)
        prefinal[j]=v
        for k in range(3):
            branches[(j,k)]=project_sector(g4,v,times.s3,k)

    keys,D=decoherence_matrix(branches)
    coherent=sum(branches.values())
    direct20=sum(at20.values())
    direct=twenty_to_third(direct20,g4,times,52,.2)
    norm2=float(np.vdot(coherent,coherent).real)
    Dn=D/norm2
    closure={
        "branch_sum_vs_direct_fidelity":fidelity(coherent,direct),
        "branch_sum_minus_direct_relative_norm":float(np.linalg.norm(coherent-direct)/np.linalg.norm(direct)),
        "coherent_final_norm2":norm2,
        "sum_all_D_normalized":complex(np.sum(Dn)),
    }

    rows=[]
    diag=Dn.diagonal().real
    for i,(j,k) in enumerate(keys):
        rows.append({
            "history":f"A->{WALL_NAMES[j]}->{WALL_NAMES[k]}",
            "second_sector":WALL_NAMES[j],
            "third_sector":WALL_NAMES[k],
            "conditional_diagonal_weight":float(diag[i]),
            "joint_A_diagonal_weight":float(first_a_weight*diag[i]),
        })

    coarse=[]
    for k in range(3):
        ids=[i for i,(j,kk) in enumerate(keys) if kk==k]
        v=sum(branches[keys[i]] for i in ids)
        pcoarse=float(np.vdot(v,v).real/norm2)
        pfine=float(np.sum(diag[ids]))
        coarse.append({
            "third_sector":WALL_NAMES[k],
            "coarse_probability":pcoarse,
            "sum_fine_diagonal_weights":pfine,
            "additivity_defect":pcoarse-pfine,
        })

    control=None
    if control_krylov is not None:
        cbranches={}
        for j in range(3):
            v=twenty_to_third(at20[j],g4,times,int(control_krylov),.2)
            for k in range(3):
                cbranches[(j,k)]=project_sector(g4,v,times.s3,k)
        ckeys,cD=decoherence_matrix(cbranches)
        cnorm=float(np.vdot(sum(cbranches.values()),sum(cbranches.values())).real)
        cDn=cD/cnorm
        control={
            "krylov_dim":int(control_krylov),
            "max_matrix_abs_difference":float(np.max(np.abs(Dn-cDn))),
            "max_diagonal_abs_difference":float(np.max(np.abs(Dn.diagonal().real-cDn.diagonal().real))),
            "max_offdiagonal_abs_difference":float(np.max(np.abs((Dn-np.diag(np.diag(Dn)))-(cDn-np.diag(np.diag(cDn)))))),
        }

    return {
        "times":times,
        "first_a_weight":first_a_weight,
        "first_non_a_weight":first_non_a_weight,
        "second_projected_raw_weights":second_weights,
        "keys":keys,
        "D_raw":D,
        "D_conditional":Dn,
        "rows":rows,
        "coarse":coarse,
        "closure":closure,
        "measures":decoherence_measures(Dn),
        "control":control,
    }
