"""Bridge audit between the existing square-root Bianchi IX branch and a naive WDW constraint.

Existing quantum Bianchi IX evolution uses

    i hbar d_s psi = G(s) psi,   G(s) = -sqrt(A(s)),

with A(s)=-hbar^2 Delta_beta + W(s,beta).

If one naively writes the second-order Wheeler-DeWitt-like constraint

    C_WDW = (i hbar d_s)^2 - A(s),

then a solution of the first-order branch obeys

    C_WDW psi = i hbar (d_s G) psi,

because G depends explicitly on s.  Therefore the current square-root branch is
not merely a coordinate rewrite of the naive second-order constraint.

This module quantifies that mismatch on the existing packet without claiming
that either quantization is uniquely physical.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh_tridiagonal

from quantum_bianchi_ix_model import QuantumBianchiIX, PacketSpec
from bianchi_ix_histories import propagate_linear


def apply_sqrt_A(model: QuantumBianchiIX, v, s: float, krylov_dim: int = 72):
    """Lanczos approximation to sqrt(A(s)) v with no state renormalization."""
    v=np.asarray(v,complex)
    beta=float(np.linalg.norm(v))
    if beta==0:
        return v.copy(), {"dimension":0,"min_ritz":0.0}
    q=v/beta
    qprev=np.zeros_like(q)
    bprev=0.0
    Q=[]; aa=[]; bb=[]
    for j in range(int(krylov_dim)):
        Q.append(q.copy())
        w=model.apply_A(q,float(s))
        a=float(np.vdot(q,w).real)
        w-=a*q+bprev*qprev
        aa.append(a)
        if j==int(krylov_dim)-1:
            break
        b=float(np.linalg.norm(w))
        bb.append(b)
        if b<1e-13:
            break
        qprev,q=q,w/b
        bprev=b
    Q=np.asarray(Q)
    aa=np.asarray(aa)
    bb=np.asarray(bb[:len(aa)-1])
    ev,U=eigh_tridiagonal(aa,bb)
    if ev.min() < -1e-8:
        raise RuntimeError(f"negative Ritz value {ev.min()}")
    coeff=U@(np.sqrt(np.maximum(ev,0.0))*U[0])
    return beta*(coeff@Q), {"dimension":int(len(aa)),"min_ritz":float(ev.min())}


def apply_G(model,v,s,krylov_dim=72):
    x,meta=apply_sqrt_A(model,v,s,krylov_dim)
    return -x,meta


def dotG_action(model,v,s,delta=0.01,krylov_dim=72):
    if delta<=0:
        raise ValueError("positive delta required")
    gp,_=apply_G(model,v,s+delta,krylov_dim)
    gm,_=apply_G(model,v,s-delta,krylov_dim)
    return (gp-gm)/(2*delta)


def state_at(model,s,ds=0.025):
    packet=PacketSpec(hbar=model.hbar,sigma=.3)
    return propagate_linear(model,model.initial(packet),2.0,float(s),float(ds),60)


def mismatch_row(model,s,delta=.01,krylov_dim=72,state=None):
    if state is None:
        state=state_at(model,s)
    dG=dotG_action(model,state,s,delta,krylov_dim)
    Av=model.apply_A(state,s)
    residual=1j*model.hbar*dG
    denom=float(np.linalg.norm(Av))
    return {
        "s":float(s),
        "delta":float(delta),
        "krylov_dim":int(krylov_dim),
        "state_norm":float(np.linalg.norm(state)),
        "Apsi_norm":denom,
        "ihbar_dotG_norm":float(np.linalg.norm(residual)),
        "relative_naive_constraint_mismatch":float(np.linalg.norm(residual)/denom) if denom else None,
        "overlap_phase_with_Apsi":{
            "re":float(np.vdot(Av,residual).real/(denom*np.linalg.norm(residual))) if denom and np.linalg.norm(residual) else 0.0,
            "im":float(np.vdot(Av,residual).imag/(denom*np.linalg.norm(residual))) if denom and np.linalg.norm(residual) else 0.0,
        },
    }


def run_audit(times=(3.0,4.6,5.5)):
    model=QuantumBianchiIX(hbar=.2,bp_bounds=(-5,5),bm_bounds=(-8,5),nx=100,ny=130,potential_cap=30.)
    rows=[];controls=[]
    for s in times:
        v=state_at(model,s)
        main=mismatch_row(model,s,.01,72,v)
        rows.append(main)
        delta_rows=[mismatch_row(model,s,d,72,v) for d in (.02,.01,.005)]
        k_rows=[mismatch_row(model,s,.01,k,v) for k in (60,72,84)]
        controls.append({"s":float(s),"delta_scan":delta_rows,"krylov_scan":k_rows})
    return {
        "schema":1,
        "model":"Existing reduced vacuum Bianchi IX square-root branch on the first finite beta box.",
        "identity":"For i*hbar*d_s psi=G(s)psi and G^2=A, [(i*hbar*d_s)^2-A]psi=i*hbar*(d_s G)psi.",
        "rows":rows,
        "controls":controls,
        "interpretation":[
            "A nonzero mismatch is structural: the time-dependent square-root branch is not exactly a solution of the naive second-order WDW constraint.",
            "This does not show that the square-root branch is wrong; it shows that switching to a timeless second-order constraint changes the quantization unless an extra term/order prescription is supplied.",
            "A Halliwell-style timeless class operator should therefore be implemented as a new explicitly declared constraint model, not advertised as a representation change of the existing finite-clock branch."
        ],
        "limitations":[
            "Only the first beta box and three packet times are sampled.",
            "d_s G is evaluated by finite differences and sqrt(A) by Lanczos; delta and Krylov scans are recorded.",
            "No timeless class operator or physical induced inner product is constructed in this audit.",
            "No black-hole observation, singularity-resolution result, or unique quantum-gravity claim."
        ]
    }
