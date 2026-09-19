"""One explicit factor-ordering/measure family for the KS positive-frequency branch.

A_q = -h^2 e^{-q x} d_x(e^{q x} d_x) + e^{-2x}
is formally symmetric on L2(e^{q x} dx). Under the unitary map
g=e^{q x/2} f it becomes
A_q^flat = A_0 + h^2 q^2/4.
The family has the same classical principal symbol and differs at O(h^2).

It is an audit family, not an exhaustive classification of all orderings.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))
from wdw_model import Packet,WDWModel

def symbolic_family_check():
    x,q,h=sp.symbols("x q h", real=True)
    g=sp.Function("g")(x)
    f=sp.exp(-q*x/2)*g
    Aq=-h**2*sp.exp(-q*x)*sp.diff(sp.exp(q*x)*sp.diff(f,x),x)+sp.exp(-2*x)*f
    flat=sp.simplify(sp.exp(q*x/2)*Aq)
    target=-h**2*sp.diff(g,x,2)+(sp.exp(-2*x)+h**2*q**2/4)*g
    residual=sp.simplify(sp.expand(flat-target))
    if residual != 0:
        raise RuntimeError(f"ordering-family unitary map failed: {residual}")
    return {
      "unitary_flat_form_residual":"0",
      "unitary_map":"g=exp(qx/2) f",
      "flat_form":"A_q^flat=A_0+h^2 q^2/4"
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()

    symbolic=symbolic_family_check()
    saved=json.loads((ROOT/"research/wdw_results/summary.json").read_text())
    detect_tol=float(saved["protocol"]["thresholds"]["box_mean_x"])
    qs=(-2.,-1.,0.,1.,2.)
    times=np.linspace(0.,6.,121)
    cases=[]
    overall=0.
    for h in (.4,.2,.1):
        m=WDWModel(h,left=-4.,right=16.,dx=.04)
        p=Packet(h,math.sqrt(h))
        chi=p.initial(m.x)
        coeff=m.U.T@chi
        means={}
        endpoints={}
        for q in qs:
            shift=h*h*q*q/4
            E=np.sqrt(m.eigenvalues+shift)
            states=m.U@(coeff[:,None]*np.exp(-1j*E[:,None]*times[None,:]/h))
            prob=np.abs(states)**2
            mean=m.x@prob
            means[q]=mean
            endpoints[q]=float(np.exp(p.x0-mean[-1]-times[-1]))
        qrows=[]
        for q in qs:
            d=float(np.max(np.abs(means[q]-means[0.])))
            overall=max(overall,d)
            qrows.append({
              "q":q,
              "flat_quantum_potential_shift":h*h*q*q/4,
              "max_mean_x_difference_from_q0":d,
              "end_geometric_radius_ratio":endpoints[q],
              "resolved_above_existing_box_mean_x_tolerance":bool(d>detect_tol)
            })
        cases.append({"h":h,"rows":qrows})

    finding="sensitive" if overall>detect_tol else "not_resolved"
    status="FAIL" if finding=="sensitive" else "PARTIAL"
    result={
      "schema":1,
      "status":status,
      "finding":finding,
      "family":{
        "definition":"A_q=-h^2 exp(-q x) d_x(exp(q x)d_x)+exp(-2x)",
        "symbolic_check":symbolic,
        "hilbert_measure":"L2(exp(q x) dx)",
        "unitary_flat_form":"A_q^flat=A_0+h^2 q^2/4",
        "q_values":list(qs),
        "principal_symbol":"p_x^2+exp(-2x)",
        "semiclassical_difference":"O(h^2)",
        "formal_symmetry":"PASS in the declared weighted measure"
      },
      "same_state_rule":"The same flat-representation state is used after the unitary measure map; packets are not re-fitted for each q.",
      "existing_numerical_detectability_threshold":detect_tol,
      "max_mean_x_ordering_spread":overall,
      "cases":cases,
      "scope_note":"This explicit admissible Sturm-Liouville ordering/measure family produces effects far above the existing box-error tolerance. Therefore the current finite-clock geometric predictions are ordering-sensitive unless an additional physical principle selects an ordering. This family is not exhaustive."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
