"""Check the independent Bianchi IX GR reduction against repository code.

The executed Cadabra result must independently recover the spatial-curvature,
ADM kinetic and Hamiltonian-constraint coefficients.  This script then checks
those coefficients against the current classical/quantum implementation.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]


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


def check_result(reference: dict, result: dict) -> dict:
    if result.get("schema") != 1:
        raise ValueError("unsupported Cadabra result schema")
    for key,value in result["checks"].items():
        if value != "0":
            raise ValueError(f"Cadabra symbolic check failed {key}={value}")
    if result["derived"]["spatial_curvature_relation"] != "R3=-12*exp(-2*alpha)*V":
        raise ValueError("unexpected spatial-curvature relation")
    if result["derived"]["adm_kinetic_invariant"] != "6/N**2*(-alpha_dot**2+beta_plus_dot**2+beta_minus_dot**2)":
        raise ValueError("unexpected ADM kinetic invariant")
    if result["derived"]["reduced_s_generator"] != "-sqrt(p_plus**2+p_minus**2+W)":
        raise ValueError("unexpected reduced generator sign")
    if result["derived"]["wall_terms"] != reference["expected"]["wall_exponential_terms"]:
        raise ValueError("Cadabra wall coefficient table does not match the declared repository convention")

    classical=load_module(ROOT/"research/bianchi_ix_model.py","astra_bianchi_ix_classical")
    quantum=load_module(ROOT/"research/quantum_bianchi_ix_model.py","astra_bianchi_ix_quantum")

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

    # Check reduced Hamilton equations including the sign from s=-alpha.
    states=[
        (2.2,np.array([0.08,-0.03,1.1,0.27])),
        (4.0,np.array([-0.16,0.12,0.86,-0.41])),
    ]
    max_rhs=0.0
    for s,y in states:
        bp,bm,pp,pm=y
        W,gp,gm=classical.wall_term_and_gradient(s,bp,bm)
        H=np.sqrt(pp*pp+pm*pm+W)
        expected=np.array([-pp/H,-pm/H,0.5*gp/H,0.5*gm/H])
        got=classical.rhs(s,y)
        max_rhs=max(max_rhs,float(np.max(np.abs(got-expected))))
    if max_rhs > 2e-13:
        raise ValueError(f"classical reduced Hamilton equations mismatch {max_rhs}")

    files={p:sha256(ROOT/p) for p in reference["source_files"]}
    return {
        "schema":1,
        "status":"PASS",
        "backend":"Cadabra",
        "comparison":"PASS",
        "source_sha256":files,
        "max_wall_relation_abs":max_wall,
        "max_classical_quantum_wall_gradient_abs":max_grad,
        "max_reduced_hamilton_rhs_abs":max_rhs,
        "verified_relations":[
            "R3=-12 exp(-2 alpha) V",
            "KijKij-K^2=6/N^2(-alpha_dot^2+beta_plus_dot^2+beta_minus_dot^2)",
            "C=1/2 exp(-3 alpha)(-p_alpha^2+p_plus^2+p_minus^2)+exp(alpha)V",
            "W=2 exp(-4s)V",
            "G_s=-sqrt(p_plus^2+p_minus^2+W)",
        ],
        "scope_note":"Agreement verifies the repository's stated classical reduction convention against one executed independent CAS backend. It does not establish uniqueness of the convention or validity of the subsequent quantization."
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--reference",type=Path,required=True)
    ap.add_argument("--cadabra-result",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    reference=json.loads(args.reference.read_text())
    result=json.loads(args.cadabra_result.read_text())
    summary=check_result(reference,result)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
