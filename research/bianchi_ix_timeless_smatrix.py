"""Small-grid timeless Bianchi IX complex-potential S-matrix pilot.

This is a NEW second-order constraint quantization motivated by Halliwell's
complex-potential class operators.  It is deliberately separated from the
existing internal-clock square-root branch, whose non-equivalence to the naive
second-order constraint is recorded in bianchi_ix_timeless_bridge.py.

Constraint on a finite Dirichlet (s,beta+,beta-) box:

    C = P_s^2 - A(s)
      = hbar^2 L_s - hbar^2(L_+ + L_-) - W(s,beta),

where L=-d^2 is the positive 1D Dirichlet Laplacian.  An absolute potential cap
is used only to keep this first scattering pilot numerically finite.

For B={B+,B-}, define V=V0 F_B >=0.  The finite-window interaction-picture
approximation to the no-entry scattering class operator is

    S_T = exp(+i C T/hbar)
          exp[-i (C - i V) (2T)/hbar]
          exp(+i C T/hbar).

The exact infinite-window S-matrix is expected to commute with C.  This module
only tests whether finite-window commutator leakage and window dependence show
a convergent regime on a small regulated grid.  It does NOT assign physical
probabilities or implement the induced inner product.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, expm_multiply

from quantum_bianchi_ix_model import wall_grid
from bianchi_ix_ever_entered import signed_B_distance


@dataclass(frozen=True)
class WDWGridSpec:
    hbar: float = 0.2
    s_bounds: tuple = (2.0, 6.0)
    bp_bounds: tuple = (-3.5, 3.5)
    bm_bounds: tuple = (-4.0, 4.0)
    ns: int = 12
    nx: int = 14
    ny: int = 14
    potential_abs_cap: float = 8.0
    absorber_width: float = 0.30


def _interior(bounds,n):
    x=np.linspace(float(bounds[0]),float(bounds[1]),int(n)+2)[1:-1]
    return x,float(x[1]-x[0])


def positive_laplacian_1d(n,dx):
    main=np.full(int(n),2.0/dx**2)
    off=np.full(int(n)-1,-1.0/dx**2)
    return sparse.diags((off,main,off),(-1,0,1),format="csr")


def build_constraint(spec=WDWGridSpec()):
    s,ds=_interior(spec.s_bounds,spec.ns)
    bp,dx=_interior(spec.bp_bounds,spec.nx)
    bm,dy=_interior(spec.bm_bounds,spec.ny)
    Ls=positive_laplacian_1d(spec.ns,ds)
    Lp=positive_laplacian_1d(spec.nx,dx)
    Lm=positive_laplacian_1d(spec.ny,dy)
    Is=sparse.eye(spec.ns,format="csr")
    Ip=sparse.eye(spec.nx,format="csr")
    Im=sparse.eye(spec.ny,format="csr")
    Ipm=sparse.kron(Ip,Im,format="csr")
    Lbeta=sparse.kron(Lp,Im,format="csr")+sparse.kron(Ip,Lm,format="csr")
    kinetic=(spec.hbar**2)*(sparse.kron(Ls,Ipm,format="csr")-sparse.kron(Is,Lbeta,format="csr"))
    S,BP,BM=np.meshgrid(s,bp,bm,indexing="ij")
    W=wall_grid(S,BP,BM)
    Wreg=np.clip(W,-float(spec.potential_abs_cap),float(spec.potential_abs_cap))
    C=kinetic-sparse.diags(Wreg.ravel(),0,format="csr")
    d=signed_B_distance(BP,BM)
    width=float(spec.absorber_width)
    if width<0: raise ValueError("absorber width must be nonnegative")
    if width==0:
        F=(d>0).astype(float)
    else:
        x=np.clip(.5+d/width,0.,1.)
        F=x*x*(3.-2.*x)
    return {
        "spec":spec,"s":s,"bp":bp,"bm":bm,"S":S,"BP":BP,"BM":BM,
        "W":W,"Wreg":Wreg,"F_B":F,"C":C,
        "shape":(spec.ns,spec.nx,spec.ny),
        "spacing":(ds,dx,dy),
    }


def near_zero_modes(C,k=4):
    n=C.shape[0]
    k=min(int(k),n-2)
    vals,vecs=eigsh(C,k=k,sigma=0.0,which="LM")
    order=np.argsort(np.abs(vals))
    vals=vals[order];vecs=vecs[:,order]
    for j in range(vecs.shape[1]):
        vecs[:,j]/=np.linalg.norm(vecs[:,j])
    return vals,vecs


def smatrix_action(C,F_B,v,T,v0,hbar):
    """Finite-window complex-potential S_T action."""
    v=np.asarray(v,complex)
    T=float(T);v0=float(v0)
    if T<0 or v0<0: raise ValueError("nonnegative T and v0 required")
    if T==0: return v.copy()
    V=sparse.diags(v0*np.asarray(F_B,float).ravel(),0,format="csr")
    x=expm_multiply((1j*T/hbar)*C,v)
    x=expm_multiply((-1j*(2*T)/hbar)*(C-1j*V),x)
    x=expm_multiply((1j*T/hbar)*C,x)
    return x


def commutator_relative(C,F_B,v,T,v0,hbar):
    Sv=smatrix_action(C,F_B,v,T,v0,hbar)
    SCv=smatrix_action(C,F_B,C@v,T,v0,hbar)
    CSv=C@Sv
    num=float(np.linalg.norm(CSv-SCv))
    den=float(np.linalg.norm(CSv)+np.linalg.norm(SCv))
    return num/den if den else 0.0


def energy_leakage(C,psi,lam):
    psi=np.asarray(psi,complex)
    return float(np.linalg.norm(C@psi-lam*psi)/max(np.linalg.norm(psi),1e-30))


def edge_mass(grid,psi,width_cells=1):
    p=np.abs(np.asarray(psi).reshape(grid["shape"]))**2
    m=np.zeros(grid["shape"],bool)
    w=int(width_cells)
    m[:w,:,:]=True;m[-w:,:,:]=True;m[:,:w,:]=True;m[:,-w:,:]=True;m[:,:,:w]=True;m[:,:,-w:]=True
    return float(p[m].sum()/p.sum())


def capped_mass(grid,psi):
    p=np.abs(np.asarray(psi).reshape(grid["shape"]))**2
    m=np.abs(grid["W"])>=.999*grid["spec"].potential_abs_cap
    return float(p[m].sum()/p.sum())


def run_pilot(spec=WDWGridSpec()):
    g=build_constraint(spec);C=g["C"];F=g["F_B"]
    herm=float(sparse.linalg.norm(C-C.getH()))
    vals,U=near_zero_modes(C,4)
    lam=float(vals[0]);psi=U[:,0]
    combo=U[:,0]+(.45j)*U[:,1]-.25*U[:,2]
    combo/=np.linalg.norm(combo)
    Ts=(.25,.5,1.0)
    v0s=(0.0,.025,.05,.10)
    rows=[]
    states={}
    for v0 in v0s:
        for T in Ts:
            y=smatrix_action(C,F,psi,T,v0,spec.hbar)
            states[(v0,T)]=y
            rows.append({
                "v0":v0,"T":T,
                "norm2":float(np.vdot(y,y).real),
                "constraint_energy_leakage":energy_leakage(C,y,lam),
                "low_energy_commutator_relative":commutator_relative(C,F,combo,T,v0,spec.hbar),
            })
    window=[]
    for v0 in (.025,.05,.10):
        prev=None
        for T in Ts:
            y=states[(v0,T)]
            if prev is not None:
                window.append({"v0":v0,"from_T":prev[0],"to_T":T,
                    "relative_state_change":float(np.linalg.norm(y-prev[1])/max(np.linalg.norm(y),1e-30))})
            prev=(T,y)
    zero=max(r["low_energy_commutator_relative"] for r in rows if r["v0"]==0.0)
    zero_id=max(float(np.linalg.norm(states[(0.0,T)]-psi)) for T in Ts)
    return {
      "schema":1,
      "constraint":"C=hbar^2 L_s-hbar^2(L_++L_-)-W_reg on a finite Dirichlet 3D minisuperspace box.",
      "grid":{"s_bounds":spec.s_bounds,"bp_bounds":spec.bp_bounds,"bm_bounds":spec.bm_bounds,
              "ns":spec.ns,"nx":spec.nx,"ny":spec.ny,"hbar":spec.hbar,
              "potential_abs_cap":spec.potential_abs_cap,"absorber_width":spec.absorber_width},
      "matrix":{"dimension":int(C.shape[0]),"hermiticity_residual":herm,
                "near_zero_eigenvalues":[float(x) for x in vals]},
      "reference_mode":{"lambda":lam,"constraint_residual":energy_leakage(C,psi,lam),
                        "edge_mass_one_cell":edge_mass(g,psi,1),"capped_region_mass":capped_mass(g,psi)},
      "rows":rows,"window_increment":window,
      "negative_control":{"v0_zero_max_commutator_relative":zero,"v0_zero_max_identity_error":zero_id},
      "interpretation_rule":"Do not call C_no or 1-C_no probabilities. First require a finite-window regime where commutator leakage and S_T changes decrease and numerical box/cap dependence is controlled.",
      "limitations":[
        "This is a new naive second-order WDW constraint quantization, not the existing finite-clock square-root model.",
        "Ordinary Euclidean grid norms are numerical diagnostics only; the induced physical inner product is not implemented.",
        "The scattering window is finite and the minisuperspace box has Dirichlet boundaries and an absolute potential cap.",
        "No claim of a converged infinite-window S-matrix, physical history probability, black-hole observation, or singularity resolution."
      ]
    }
