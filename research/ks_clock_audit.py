"""Classical clock audit for the Schwarzschild/Kantowski-Sachs WDW constraint.

This audits admissibility of intrinsic clocks on the classical constraint surface.
It does not claim quantum equivalence between different clock quantizations.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()

    T,x,pT,px=sp.symbols("T x p_T p_x", real=True)
    C=pT**2-px**2-sp.exp(-2*x)

    def pb(f,g):
        return sp.simplify(
            sp.diff(f,T)*sp.diff(g,pT)-sp.diff(f,pT)*sp.diff(g,T)
            +sp.diff(f,x)*sp.diff(g,px)-sp.diff(f,px)*sp.diff(g,x)
        )

    Y=x+T
    r=sp.exp(-Y)/4
    pb_T=sp.simplify(pb(T,C))
    pb_x=sp.simplify(pb(x,C))
    pb_Y=sp.simplify(pb(Y,C))
    pb_r=sp.simplify(pb(r,C))

    exact={
      "{T,C}":str(pb_T),
      "{x,C}":str(pb_x),
      "{Y,C}":str(pb_Y),
      "{r,C}":str(pb_r),
    }
    if pb_T != 2*pT or pb_x != -2*px or pb_Y != 2*(pT-px):
        raise RuntimeError(exact)

    # Negative-frequency branch used by ih d_T psi = +H psi has p_T=-H.
    # H>|p_x| for every finite x because exp(-2x)>0.
    probe_x=np.array([-4.,0.,2.,8.,16.])
    probe_p=np.array([-5.,-1.,0.,1.,5.])
    min_minus_T=float("inf")
    min_minus_Y=float("inf")
    for xv in probe_x:
        for pv in probe_p:
            H=math.sqrt(pv*pv+math.exp(-2*xv))
            # Magnitudes of negative branch brackets.
            min_minus_T=min(min_minus_T,2*H)
            min_minus_Y=min(min_minus_Y,2*(H+pv))
            if not H>abs(pv):
                raise RuntimeError("finite-x positivity failure")
    if min_minus_T<=0 or min_minus_Y<=0:
        raise RuntimeError("clock bracket vanished on finite-x probe")

    # Existing exact classical paths: x has a turning point, while
    # T and Y=x+T remain monotone and r remains monotone decreasing.
    sys_path=ROOT/"research"
    import sys
    sys.path.insert(0,str(sys_path))
    from wdw_model import classical_path

    times=np.linspace(0.,10.,401)
    path_checks=[]
    for x0,p0 in [(2.,-1.),(1.2,-.6),(2.5,-1.4)]:
        xx,pp=classical_path(x0,p0,times)
        yy=xx+times
        rr=np.exp(-yy)/4
        energy=math.sqrt(p0*p0+math.exp(-2*x0))
        # Branch convention pT=-energy.
        t_bracket=-2*energy*np.ones_like(times)
        y_bracket=2*(-energy-pp)
        x_bracket=-2*pp
        if not np.all(t_bracket<0):
            raise RuntimeError("T clock lost monotonicity")
        if not np.all(y_bracket<0):
            raise RuntimeError("areal-radius clock lost monotonicity")
        if not np.all(np.diff(yy)>0) or not np.all(np.diff(rr)<0):
            raise RuntimeError("Y/r monotonicity failure")
        if np.min(np.abs(x_bracket))>5e-2:
            # Samples need only bracket the x turning point closely.
            raise RuntimeError("x-clock turning point not resolved")
        path_checks.append({
          "x0":x0,"p0":p0,
          "min_abs_T_bracket":float(np.min(np.abs(t_bracket))),
          "min_abs_Y_bracket":float(np.min(np.abs(y_bracket))),
          "min_abs_x_bracket":float(np.min(np.abs(x_bracket))),
          "Y_monotone_increasing":bool(np.all(np.diff(yy)>0)),
          "r_monotone_decreasing":bool(np.all(np.diff(rr)<0))
        })

    result={
      "schema":1,
      "status":"PARTIAL",
      "classical_T_clock":"PASS",
      "classical_areal_clock_Y":"PASS",
      "x_as_global_clock":"FAIL",
      "quantum_multi_clock_equivalence":"PENDING",
      "branch":{
        "repo_positive_frequency_equation":"i h d_T psi = H psi",
        "canonical_momentum_branch":"p_T=-H",
        "H":"sqrt(p_x^2+exp(-2x))",
        "finite_x_branch_gap":"H>|p_x|, so p_T never crosses zero"
      },
      "exact_poisson_brackets":exact,
      "alternative_clock":{
        "Y":"x+T=-log(4r)",
        "negative_branch_bracket":"{Y,C}=-2(H+p_x)<0 at finite x",
        "interpretation":"Y is the logarithmic inverse areal radius and is classically monotone through the Schwarzschild interior."
      },
      "negative_control":"x is not a global clock because {x,C}=-2 p_x vanishes at the classical reflection point p_x=0.",
      "path_checks":path_checks,
      "scope_note":"Classical clock admissibility is closed for T and Y. The audit remains PARTIAL because the same quantum state/observable has not yet been quantized and compared using Y as the clock."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
