"""Parse the independent Maxima Bianchi IX GR reduction and compare repo code.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

REQUIRED_ZERO_MARKERS=[
    "ASTRA_MAXIMA_R3_CLOSED_DELTA",
    "ASTRA_MAXIMA_R3_MISNER_DELTA",
    "ASTRA_MAXIMA_KINETIC_DELTA",
    "ASTRA_MAXIMA_PA_DELTA",
    "ASTRA_MAXIMA_PP_DELTA",
    "ASTRA_MAXIMA_PM_DELTA",
    "ASTRA_MAXIMA_CONSTRAINT_DELTA",
]

def load_module(path: Path, name: str):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def parse_markers(text: str) -> dict[str,str]:
    out={}
    for line in text.splitlines():
        line=line.strip()
        if line.startswith("ASTRA_MAXIMA_") and "=" in line:
            key,value=line.split("=",1)
            out[key]=value.strip()
    if "ASTRA_MAXIMA_BIANCHI_IX_OK" not in text:
        raise ValueError("Maxima success marker missing")
    for key in REQUIRED_ZERO_MARKERS:
        actual=out.get(key)
        if actual != "0":
            raise ValueError(f"Maxima symbolic check failed {key}={actual!r}")
    return out

def compare_repository() -> dict:
    classical=load_module(ROOT/"research/bianchi_ix_model.py","astra_bianchi_ix_classical_maxima")
    quantum=load_module(ROOT/"research/quantum_bianchi_ix_model.py","astra_bianchi_ix_quantum_maxima")
    points=[
        (2.0,0.0,0.0),
        (2.3,0.17,-0.08),
        (3.1,-0.22,0.11),
        (4.6,0.35,0.19),
        (5.5,-0.41,-0.23),
    ]
    max_wall=0.0
    max_grad=0.0
    for s,bp,bm in points:
        V=classical.potential(bp,bm)
        expected_w=2*np.exp(-4*s)*V
        w,gp,gm=classical.wall_term_and_gradient(s,bp,bm)
        q_w,q_gp,q_gm=quantum._wall_vec(s,bp,bm)
        max_wall=max(max_wall,abs(w-expected_w),abs(q_w-expected_w))
        max_grad=max(max_grad,abs(w-q_w),abs(gp-q_gp),abs(gm-q_gm))
    if max_wall > 2e-13 or max_grad > 2e-13:
        raise ValueError(f"repository wall implementation mismatch wall={max_wall} grad={max_grad}")

    max_rhs=0.0
    for s,y in [
        (2.2,np.array([0.08,-0.03,1.1,0.27])),
        (4.0,np.array([-0.16,0.12,0.86,-0.41])),
    ]:
        bp,bm,pp,pm=y
        W,gp,gm=classical.wall_term_and_gradient(s,bp,bm)
        H=np.sqrt(pp*pp+pm*pm+W)
        expected=np.array([-pp/H,-pm/H,0.5*gp/H,0.5*gm/H])
        got=classical.rhs(s,y)
        max_rhs=max(max_rhs,float(np.max(np.abs(got-expected))))
    if max_rhs > 2e-13:
        raise ValueError(f"classical reduced Hamilton equations mismatch {max_rhs}")

    return {
        "max_wall_relation_abs":max_wall,
        "max_classical_quantum_wall_gradient_abs":max_grad,
        "max_reduced_hamilton_rhs_abs":max_rhs,
        "source_sha256":{
            "research/bianchi_ix_model.py":sha256(ROOT/"research/bianchi_ix_model.py"),
            "research/quantum_bianchi_ix_model.py":sha256(ROOT/"research/quantum_bianchi_ix_model.py"),
        },
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--maxima-log",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()

    markers=parse_markers(args.maxima_log.read_text(errors="replace"))
    repo=compare_repository()
    result={
        "schema":1,
        "status":"PASS",
        "backend":"Maxima/ctensor",
        "version":"5.46.0",
        "comparison":"PASS",
        "symbolic_zero_markers":{k:markers[k] for k in REQUIRED_ZERO_MARKERS},
        **repo,
        "verified_relations":[
            "R3=-12 exp(-2 alpha) V",
            "KijKij-K^2=6/N^2(-alpha_dot^2+beta_plus_dot^2+beta_minus_dot^2)",
            "C=1/2 exp(-3 alpha)(-p_alpha^2+p_plus^2+p_minus^2)+exp(alpha)V",
            "W=2 exp(-4s)V",
            "G_s=-sqrt(p_plus^2+p_minus^2+W)",
        ],
        "scope_note":"Agreement provides a second executed CAS backend independent of the Cadabra/SymPy path for the stated classical reduction convention. It does not establish unique quantization or black-hole singularity resolution."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
