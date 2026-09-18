"""Stationary resolvent/T-operator pilot for the timeless Bianchi IX constraint.

This module follows the negative finite-window Dirichlet S-matrix pilot.  It
keeps the same second-order constraint discretization but switches the
numerical question from a long-time box evolution to a stationary resolvent.

For the physical-region absorber V >= 0,

    H_eff = C - i V,
    G0(z) = (z - C)^(-1),
    GV(z) = (z - C + i V)^(-1),   z = E + i eta, eta > 0.

The finite-dimensional identities tested here are

    GV = G0 + G0 U GV,            U = -i V
    T  = U + U GV U,

with no claim that a finite-box eta-regularized T is already Halliwell's
infinite-volume on-shell S-matrix.  The pilot asks only whether the T-action
has a numerically stable trend as eta is reduced while boundary support stays
controlled.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from bianchi_ix_timeless_smatrix import (
    WDWGridSpec,
    build_constraint,
    edge_mass,
    low_energy_edge_min_packet,
    near_zero_modes,
)


@dataclass(frozen=True)
class ResolventSpec:
    energy: float = 0.0
    etas: tuple = (0.10, 0.05, 0.025)
    v0s: tuple = (0.025, 0.05)
    source_modes: int = 16
    shell_modes: int = 12


def _factor(C, energy, eta, v_diag=None):
    eta=float(eta)
    if eta <= 0:
        raise ValueError("eta must be positive")
    n=C.shape[0]
    A=(complex(energy,eta)*sparse.eye(n,format="csc")-C.astype(complex).tocsc())
    if v_diag is not None:
        vd=np.asarray(v_diag,float).ravel()
        if vd.size != n:
            raise ValueError("v_diag shape mismatch")
        A=A+1j*sparse.diags(vd,0,format="csc")
    return splu(A)


def resolvent_action(C,b,energy,eta,v_diag=None):
    b=np.asarray(b,complex)
    return _factor(C,energy,eta,v_diag).solve(b)


def dyson_relative(C,b,energy,eta,v_diag):
    """Check GV b = G0 b + G0 U GV b for U=-iV."""
    vd=np.asarray(v_diag,float).ravel()
    free=_factor(C,energy,eta,None)
    full=_factor(C,energy,eta,vd)
    b=np.asarray(b,complex)
    x0=free.solve(b)
    x=full.solve(b)
    rhs=x0+free.solve((-1j*vd)*x)
    return float(np.linalg.norm(x-rhs)/max(np.linalg.norm(x),1e-30))


def t_action(C,phi,energy,eta,v_diag):
    """Apply T(z)=U+U GV(z) U, U=-iV, to one vector."""
    vd=np.asarray(v_diag,float).ravel()
    phi=np.asarray(phi,complex)
    uphi=(-1j*vd)*phi
    if not np.any(vd):
        return np.zeros_like(phi),np.zeros_like(phi)
    x=_factor(C,energy,eta,vd).solve(uphi)
    return uphi+(-1j*vd)*x,x


def shell_fraction(shell_vectors,v):
    v=np.asarray(v,complex)
    den=float(np.vdot(v,v).real)
    if den <= 0:
        return 0.0
    a=shell_vectors.conj().T@v
    return float(np.vdot(a,a).real/den)


def run_pilot(grid_spec=WDWGridSpec(), spec=ResolventSpec()):
    g=build_constraint(grid_spec)
    C=g["C"]
    shell_vals,shell=near_zero_modes(C,spec.shell_modes)
    _,psi,packet=low_energy_edge_min_packet(g,spec.source_modes)
    source_residual=float(np.linalg.norm(C@psi)/max(np.linalg.norm(psi),1e-30))

    rows=[]
    tstates={}
    for v0 in spec.v0s:
        vd=float(v0)*np.asarray(g["F_B"],float).ravel()
        born=(-1j*vd)*psi
        born_norm=float(np.linalg.norm(born))
        for eta in spec.etas:
            t,response=t_action(C,psi,spec.energy,eta,vd)
            tstates[(float(v0),float(eta))]=t
            rows.append({
                "v0":float(v0),
                "eta":float(eta),
                "dyson_relative":dyson_relative(C,psi,spec.energy,eta,vd),
                "free_resolvent_norm":float(np.linalg.norm(resolvent_action(C,psi,spec.energy,eta))),
                "t_norm":float(np.linalg.norm(t)),
                "t_to_born_norm_ratio":float(np.linalg.norm(t)/max(born_norm,1e-30)),
                "t_shell_fraction":shell_fraction(shell,t),
                "response_norm":float(np.linalg.norm(response)),
                "response_edge_mass_one_cell":edge_mass(g,response,1),
            })

    eta_change=[]
    etas=tuple(float(x) for x in spec.etas)
    for v0 in spec.v0s:
        for hi,lo in zip(etas[:-1],etas[1:]):
            a=tstates[(float(v0),hi)]
            b=tstates[(float(v0),lo)]
            eta_change.append({
                "v0":float(v0),
                "from_eta":hi,
                "to_eta":lo,
                "relative_t_change":float(np.linalg.norm(b-a)/max(np.linalg.norm(b),1e-30)),
            })

    zero=np.zeros(C.shape[0],float)
    t0,_=t_action(C,psi,spec.energy,etas[0],zero)
    return {
        "schema":1,
        "method":"Stationary eta-regularized resolvent/T-operator pilot; not yet an on-shell S-matrix.",
        "grid":{
            "s_bounds":grid_spec.s_bounds,"bp_bounds":grid_spec.bp_bounds,"bm_bounds":grid_spec.bm_bounds,
            "ns":grid_spec.ns,"nx":grid_spec.nx,"ny":grid_spec.ny,"hbar":grid_spec.hbar,
            "potential_abs_cap":grid_spec.potential_abs_cap,"absorber_width":grid_spec.absorber_width,
            "dimension":int(C.shape[0]),
        },
        "resolvent":{
            "energy":float(spec.energy),"etas":[float(x) for x in spec.etas],
            "v0s":[float(x) for x in spec.v0s],
            "source_modes":int(spec.source_modes),"shell_modes":int(spec.shell_modes),
        },
        "reference_packet":{
            **packet,
            "constraint_zero_residual":source_residual,
            "shell_eigenvalues":[float(x) for x in shell_vals[:8]],
        },
        "negative_control":{
            "v0_zero_t_norm":float(np.linalg.norm(t0)),
        },
        "rows":rows,
        "eta_change":eta_change,
        "decision_rule":"Do not construct physical history probabilities. A next-stage scattering calculation requires an eta trend with decreasing T-action changes, small algebraic Dyson residuals, controlled response edge mass, and later independent box/cap scans.",
        "limitations":[
            "Same new second-order WDW constraint quantization as the preceding timeless pilot; not the finite-clock square-root model.",
            "The finite Dirichlet box and absolute potential cap remain regulators.",
            "eta is a stationary outgoing-resolvent regulator, not a physical parameter.",
            "Euclidean grid norms and shell projections are numerical diagnostics, not the induced physical inner product.",
            "No converged infinite-volume S-matrix, decoherence functional, history probability, black-hole observation, or singularity-resolution claim is made.",
        ],
    }
