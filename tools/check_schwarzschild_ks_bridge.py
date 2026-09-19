"""Check the Schwarzschild-interior -> KS -> repo WDW bridge.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def load_module(path: Path, name: str):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    mod=importlib.util.module_from_spec(spec)
    sys.modules[name]=mod
    spec.loader.exec_module(mod)
    return mod

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def parse_maxima(text: str):
    required=[
      "ASTRA_KS_MAXIMA_R3_DELTA",
      "ASTRA_KS_MAXIMA_KINETIC_DELTA",
      "ASTRA_KS_MAXIMA_PO_DELTA",
      "ASTRA_KS_MAXIMA_PB_DELTA",
      "ASTRA_KS_MAXIMA_HAMILTONIAN_DELTA",
      "ASTRA_KS_MAXIMA_REPO_CONSTRAINT_DELTA",
      "ASTRA_KS_MAXIMA_PB_X_DELTA",
      "ASTRA_KS_MAXIMA_PB_T_DELTA",
      "ASTRA_KS_MAXIMA_SCHWARZSCHILD_CONSTRAINT_DELTA",
      "ASTRA_KS_MAXIMA_RADIUS_DELTA",
      "ASTRA_KS_MAXIMA_MASS_DELTA",
      "ASTRA_KS_MAXIMA_MASS_PB_DELTA",
    ]
    vals={}
    for line in text.splitlines():
        line=line.strip()
        if line.startswith("ASTRA_KS_MAXIMA_") and "=" in line:
            k,v=line.split("=",1); vals[k]=v.strip()
    if "ASTRA_KS_MAXIMA_BRIDGE_OK" not in text:
        raise ValueError("Maxima bridge success marker missing")
    for k in required:
        if vals.get(k)!="0":
            raise ValueError(f"Maxima check failed {k}={vals.get(k)!r}")
    return vals

def check_cadabra(result: dict):
    if result.get("schema")!=1 or result.get("backend",{}).get("name")!="Cadabra":
        raise ValueError("unexpected Cadabra bridge result")
    for k,v in result.get("checks",{}).items():
        if v!="0":
            raise ValueError(f"Cadabra check failed {k}={v}")

def repo_checks():
    wdw=load_module(ROOT/"research/wdw_model.py","astra_wdw_bridge")
    cases=[
      (2.0,-1.0),
      (1.2,-0.6),
      (2.5,-1.4),
    ]
    max_constraint=0.0
    max_radius_relation=0.0
    max_mass_drift=0.0
    max_schwarzschild_relation=0.0
    times=np.linspace(0,6,61)
    for x0,p0 in cases:
        x,p=wdw.classical_path(x0,p0,times)
        E=math.sqrt(p0*p0+math.exp(-2*x0))
        u0=math.asinh(p0*math.exp(x0))
        mu=E*math.exp(u0)/4
        r=np.exp(-x-times)/4
        constraint=p*p+np.exp(-2*x)-E*E
        max_constraint=max(max_constraint,float(np.max(np.abs(constraint))))
        schwarz=1/(1+np.exp(2*(times+u0)))
        max_schwarzschild_relation=max(
            max_schwarzschild_relation,
            float(np.max(np.abs(r/(2*mu)-schwarz)))
        )
        pT=-E
        mass=np.exp(x-times)*pT*(pT-p)/4
        max_mass_drift=max(
            max_mass_drift,
            float(np.max(np.abs(mass-mu)))
        )
        max_radius_relation=max(
            max_radius_relation,
            float(np.max(np.abs(r-0.25*np.exp(-x-times))))
        )
    tol=5e-13
    if max(max_constraint,max_radius_relation,max_mass_drift,max_schwarzschild_relation)>tol:
        raise ValueError(
          f"repo bridge mismatch C={max_constraint} r={max_radius_relation} "
          f"M={max_mass_drift} Schw={max_schwarzschild_relation}"
        )
    return {
      "max_classical_constraint_abs":max_constraint,
      "max_radius_identity_abs":max_radius_relation,
      "max_mass_observable_drift_abs":max_mass_drift,
      "max_schwarzschild_radius_relation_abs":max_schwarzschild_relation,
      "wdw_model_sha256":sha256(ROOT/"research/wdw_model.py")
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--reference",type=Path,required=True)
    ap.add_argument("--maxima-log",type=Path,required=True)
    ap.add_argument("--cadabra-result",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    reference=json.loads(args.reference.read_text())
    maxima=parse_maxima(args.maxima_log.read_text(errors="replace"))
    cadabra=json.loads(args.cadabra_result.read_text())
    check_cadabra(cadabra)
    repo=repo_checks()
    out={
      "schema":1,
      "status":"PASS",
      "target":reference["target"],
      "maxima":"PASS",
      "cadabra":"PASS",
      "repository_comparison":"PASS",
      **repo,
      "verified_relations":[
        "Schwarzschild interior is KS with N=F^-1/2, A=lambda*F^1/2, R=rho",
        "R3=2/R^2",
        "P_Omega^2-P_beta^2+48 exp(-2 sqrt(3) Omega)=0",
        "C_repo=-(1/3) C_Misner",
        "r=(1/4) exp(-x-T)",
        "mu_D=(1/4) exp(x-T) p_T (p_T-p_x)",
        "{mu_D,C_repo}=-(p_T/2) exp(x-T) C_repo"
      ],
      "scope_note":"This closes the classical Schwarzschild/KS bridge only. Quantum clock, inner product, operator domain, ordering and curvature-operator obligations remain separate."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
