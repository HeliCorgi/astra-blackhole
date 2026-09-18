"""Outer-boundary absorbing-layer pilot for the timeless Bianchi IX resolvent.

The physical history absorber V_B is kept distinct from a numerical outer
absorber W_out. The latter is introduced only to test whether the unstable
eta->0 trend of the finite-Dirichlet stationary resolvent is dominated by the
reflecting outer box.

Use the common background

    H_bg = C - i W_out

for both background and B-absorbing resolvents, and define

    G_bg(z) = (z - C + i W_out)^(-1)
    G_full(z) = (z - C + i W_out + i V_B)^(-1).

With U_B=-i V_B,

    G_full = G_bg + G_bg U_B G_full
    T_B = U_B + U_B G_full U_B.

W_out is never interpreted as a physical potential.
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
class OuterCAPSpec:
    energy: float = 0.0
    etas: tuple = (0.10, 0.05, 0.025)
    v0s: tuple = (0.025, 0.05)
    outer_strengths: tuple = (0.0, 0.05, 0.10, 0.20)
    layer_cells: int = 3
    source_modes: int = 20
    shell_modes: int = 12


def _smoothstep(x):
    x=np.clip(np.asarray(x,float),0.0,1.0)
    return x*x*(3.0-2.0*x)


def outer_cap_profile(shape,layer_cells=3):
    """Dimensionless 0..1 absorber profile supported near every outer face."""
    shape=tuple(int(n) for n in shape)
    w=int(layer_cells)
    if w < 1:
        raise ValueError("layer_cells must be >= 1")
    axes=[]
    for n in shape:
        if n < 2:
            raise ValueError("each dimension must contain at least two cells")
        i=np.arange(n)
        d=np.minimum(i,n-1-i)
        r=np.clip((w-d)/float(w),0.0,1.0)
        axes.append(_smoothstep(r))
    S,P,M=np.meshgrid(axes[0],axes[1],axes[2],indexing="ij")
    return np.maximum(np.maximum(S,P),M)


def _factor(C,energy,eta,outer_diag=None,v_diag=None):
    eta=float(eta)
    if eta <= 0:
        raise ValueError("eta must be positive")
    n=C.shape[0]
    A=complex(energy,eta)*sparse.eye(n,format="csc")-C.astype(complex).tocsc()
    if outer_diag is not None:
        od=np.asarray(outer_diag,float).ravel()
        if od.size != n or np.min(od) < 0:
            raise ValueError("invalid outer absorber")
        A=A+1j*sparse.diags(od,0,format="csc")
    if v_diag is not None:
        vd=np.asarray(v_diag,float).ravel()
        if vd.size != n or np.min(vd) < 0:
            raise ValueError("invalid B absorber")
        A=A+1j*sparse.diags(vd,0,format="csc")
    return splu(A)


def dyson_relative(C,b,energy,eta,outer_diag,v_diag):
    b=np.asarray(b,complex)
    vd=np.asarray(v_diag,float).ravel()
    bg=_factor(C,energy,eta,outer_diag,None)
    full=_factor(C,energy,eta,outer_diag,vd)
    x0=bg.solve(b)
    x=full.solve(b)
    rhs=x0+bg.solve((-1j*vd)*x)
    return float(np.linalg.norm(x-rhs)/max(np.linalg.norm(x),1e-30))


def t_action(C,phi,energy,eta,outer_diag,v_diag):
    phi=np.asarray(phi,complex)
    vd=np.asarray(v_diag,float).ravel()
    uphi=(-1j*vd)*phi
    if not np.any(vd):
        return np.zeros_like(phi),np.zeros_like(phi)
    x=_factor(C,energy,eta,outer_diag,vd).solve(uphi)
    return uphi+(-1j*vd)*x,x


def background_resolvent_action(C,phi,energy,eta,outer_diag):
    return _factor(C,energy,eta,outer_diag,None).solve(np.asarray(phi,complex))


def weighted_mass(weight,psi):
    p=np.abs(np.asarray(psi,complex).ravel())**2
    den=float(p.sum())
    if den <= 0:
        return 0.0
    w=np.asarray(weight,float).ravel()
    return float(np.dot(w,p)/den)


def shell_fraction(shell_vectors,v):
    v=np.asarray(v,complex)
    den=float(np.vdot(v,v).real)
    if den <= 0:
        return 0.0
    a=shell_vectors.conj().T@v
    return float(np.vdot(a,a).real/den)


def run_pilot(grid_spec=WDWGridSpec(),spec=OuterCAPSpec()):
    g=build_constraint(grid_spec)
    C=g["C"]
    profile=outer_cap_profile(g["shape"],spec.layer_cells).ravel()
    shell_vals,shell=near_zero_modes(C,spec.shell_modes)
    _,psi,packet=low_energy_edge_min_packet(g,spec.source_modes)

    rows=[]
    tstates={}
    for gamma in spec.outer_strengths:
        outer=float(gamma)*profile
        for v0 in spec.v0s:
            vd=float(v0)*np.asarray(g["F_B"],float).ravel()
            born=(-1j*vd)*psi
            born_norm=float(np.linalg.norm(born))
            for eta in spec.etas:
                bg=background_resolvent_action(C,psi,spec.energy,eta,outer)
                t,response=t_action(C,psi,spec.energy,eta,outer,vd)
                key=(float(gamma),float(v0),float(eta))
                tstates[key]=t
                rows.append({
                    "outer_strength":float(gamma),
                    "v0":float(v0),
                    "eta":float(eta),
                    "dyson_relative":dyson_relative(C,psi,spec.energy,eta,outer,vd),
                    "background_resolvent_norm":float(np.linalg.norm(bg)),
                    "background_edge_mass_one_cell":edge_mass(g,bg,1),
                    "background_outer_layer_weight":weighted_mass(profile,bg),
                    "t_norm":float(np.linalg.norm(t)),
                    "t_to_born_norm_ratio":float(np.linalg.norm(t)/max(born_norm,1e-30)),
                    "t_shell_fraction":shell_fraction(shell,t),
                    "response_norm":float(np.linalg.norm(response)),
                    "response_edge_mass_one_cell":edge_mass(g,response,1),
                    "response_outer_layer_weight":weighted_mass(profile,response),
                })

    eta_change=[]
    for gamma in spec.outer_strengths:
        for v0 in spec.v0s:
            es=tuple(float(x) for x in spec.etas)
            for hi,lo in zip(es[:-1],es[1:]):
                a=tstates[(float(gamma),float(v0),hi)]
                b=tstates[(float(gamma),float(v0),lo)]
                eta_change.append({
                    "outer_strength":float(gamma),
                    "v0":float(v0),
                    "from_eta":hi,
                    "to_eta":lo,
                    "relative_t_change":float(np.linalg.norm(b-a)/max(np.linalg.norm(b),1e-30)),
                })

    outer_change=[]
    gs=tuple(float(x) for x in spec.outer_strengths)
    for v0 in spec.v0s:
        for eta in spec.etas:
            for g0,g1 in zip(gs[:-1],gs[1:]):
                a=tstates[(g0,float(v0),float(eta))]
                b=tstates[(g1,float(v0),float(eta))]
                outer_change.append({
                    "v0":float(v0),
                    "eta":float(eta),
                    "from_outer_strength":g0,
                    "to_outer_strength":g1,
                    "relative_t_change":float(np.linalg.norm(b-a)/max(np.linalg.norm(b),1e-30)),
                })

    zero=np.zeros(C.shape[0],float)
    t0,_=t_action(
        C,psi,spec.energy,float(spec.etas[0]),
        float(spec.outer_strengths[-1])*profile,zero)
    return {
        "schema":1,
        "method":"Stationary resolvent with a separate numerical outer absorbing layer; not a PML and not an on-shell S-matrix.",
        "grid":{
            "s_bounds":grid_spec.s_bounds,
            "bp_bounds":grid_spec.bp_bounds,
            "bm_bounds":grid_spec.bm_bounds,
            "ns":grid_spec.ns,"nx":grid_spec.nx,"ny":grid_spec.ny,
            "hbar":grid_spec.hbar,
            "potential_abs_cap":grid_spec.potential_abs_cap,
            "absorber_width":grid_spec.absorber_width,
            "dimension":int(C.shape[0]),
        },
        "scan":{
            "energy":float(spec.energy),
            "etas":[float(x) for x in spec.etas],
            "v0s":[float(x) for x in spec.v0s],
            "outer_strengths":[float(x) for x in spec.outer_strengths],
            "layer_cells":int(spec.layer_cells),
            "source_modes":int(spec.source_modes),
            "shell_modes":int(spec.shell_modes),
        },
        "reference_packet":{
            **packet,
            "constraint_zero_residual":float(np.linalg.norm(C@psi)/max(np.linalg.norm(psi),1e-30)),
            "shell_eigenvalues":[float(x) for x in shell_vals[:8]],
            "outer_profile_weight":weighted_mass(profile,psi),
        },
        "negative_control":{"v0_zero_t_norm":float(np.linalg.norm(t0))},
        "rows":rows,
        "eta_change":eta_change,
        "outer_change":outer_change,
        "decision_rule":"The outer absorber is useful only if eta-successive T-action changes decrease over a non-tuned range of outer strengths and the result is not strongly dependent on outer strength/layer settings. Otherwise keep it as a negative numerical boundary control.",
        "limitations":[
            "The outer absorber is a numerical boundary CAP, distinct from the B-region complex potential and not a physical interaction.",
            "This is not a perfectly matched layer and no reflection coefficient has yet been measured independently.",
            "The finite Dirichlet box remains behind the absorbing layer.",
            "Euclidean grid norms and shell projections are numerical diagnostics, not the induced physical inner product.",
            "No converged infinite-volume S-matrix, decoherence functional, history probability, black-hole observation, or singularity-resolution claim is made.",
        ],
    }
