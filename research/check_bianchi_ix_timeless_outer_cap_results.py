from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np


def _key(gamma, v0):
    return f"g={gamma:.3f}|v={v0:.3f}"


def _vkey(v0):
    return f"v={v0:.3f}"


def _lkey(layer, gamma, v0):
    return f"L={int(layer)}|g={gamma:.3f}|v={v0:.3f}"


def compact_from_artifact(artifact_dir: Path):
    primary=json.loads((artifact_dir/"summary.json").read_text())
    layers=json.loads((artifact_dir/"layer_scan.json").read_text())

    eta={}
    for row in primary["eta_change"]:
        k=_key(float(row["outer_strength"]),float(row["v0"]))
        eta.setdefault(k,[]).append((float(row["from_eta"]),float(row["relative_t_change"])))
    eta={k:[v for _,v in sorted(vals,reverse=True)] for k,vals in eta.items()}

    strength={}
    for row in primary["outer_change"]:
        if abs(float(row["eta"])-0.025)>1e-12:
            continue
        k=_vkey(float(row["v0"]))
        strength.setdefault(k,[]).append((float(row["from_outer_strength"]),float(row["relative_t_change"])))
    strength={k:[v for _,v in sorted(vals)] for k,vals in strength.items()}

    selected={}
    for row in primary["rows"]:
        eta0=float(row["eta"])
        gamma=float(row["outer_strength"])
        if abs(eta0-0.025)>1e-12 or gamma not in (0.0,0.2):
            continue
        k=_key(gamma,float(row["v0"]))
        selected[k]={
            "background_resolvent_norm":float(row["background_resolvent_norm"]),
            "background_edge_mass_one_cell":float(row["background_edge_mass_one_cell"]),
            "response_edge_mass_one_cell":float(row["response_edge_mass_one_cell"]),
            "t_norm":float(row["t_norm"]),
        }

    layer_eta={}
    for row in layers["eta_change"]:
        k=_lkey(int(row["layer_cells"]),float(row["outer_strength"]),float(row["v0"]))
        layer_eta.setdefault(k,[]).append((float(row["from_eta"]),float(row["relative_t_change"])))
    layer_eta={k:[v for _,v in sorted(vals,reverse=True)] for k,vals in layer_eta.items()}

    layer_change={}
    for row in layers["layer_change"]:
        if abs(float(row["eta"])-0.025)>1e-12:
            continue
        k=_key(float(row["outer_strength"]),float(row["v0"]))
        layer_change.setdefault(k,[]).append((int(row["from_layer_cells"]),float(row["relative_t_change"])))
    layer_change={k:[v for _,v in sorted(vals)] for k,vals in layer_change.items()}

    return {
        "reference":{
            "edge_mass_one_cell":float(primary["reference_packet"]["edge_mass_one_cell"]),
            "constraint_zero_residual":float(primary["reference_packet"]["constraint_zero_residual"]),
            "layer_profile_weights":{str(k):float(v) for k,v in layers["reference_packet"]["profile_weights"].items()},
        },
        "primary_eta_change":eta,
        "primary_strength_change_eta_0p025":strength,
        "selected_response_edge_eta_0p025":selected,
        "layer_eta_change":layer_eta,
        "layer_change_eta_0p025":layer_change,
    }


def _flatten(obj, prefix=""):
    out=[]
    if isinstance(obj,dict):
        for key in sorted(obj):
            out.extend(_flatten(obj[key],f"{prefix}/{key}"))
    elif isinstance(obj,list):
        for i,value in enumerate(obj):
            out.extend(_flatten(value,f"{prefix}/{i}"))
    elif isinstance(obj,(int,float)):
        out.append((prefix,float(obj)))
    else:
        raise TypeError(f"unsupported value at {prefix}: {type(obj)!r}")
    return out


def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument("--artifact-dir",type=Path,required=True)
    ap.add_argument("--baseline",type=Path,required=True)
    args=ap.parse_args(argv)

    baseline=json.loads(args.baseline.read_text())
    expected={
        k:baseline[k]
        for k in (
            "reference",
            "primary_eta_change",
            "primary_strength_change_eta_0p025",
            "selected_response_edge_eta_0p025",
            "layer_eta_change",
            "layer_change_eta_0p025",
        )
    }
    actual=compact_from_artifact(args.artifact_dir)
    exp=_flatten(expected)
    got=_flatten(actual)
    exp_paths=[p for p,_ in exp]
    got_paths=[p for p,_ in got]
    if exp_paths != got_paths:
        missing=sorted(set(exp_paths)-set(got_paths))
        extra=sorted(set(got_paths)-set(exp_paths))
        raise RuntimeError(f"outer CAP result structure mismatch missing={missing} extra={extra}")
    x=np.asarray([v for _,v in exp],float)
    y=np.asarray([v for _,v in got],float)
    if not np.allclose(x,y,rtol=3e-5,atol=3e-8):
        i=int(np.argmax(np.abs(x-y)))
        raise RuntimeError(
            f"outer CAP regression mismatch at {exp_paths[i]}: expected={x[i]} got={y[i]}"
        )
    print(json.dumps({
        "checked_numeric_values":int(x.size),
        "max_abs_difference":float(np.max(np.abs(x-y))),
        "status":"ok",
    },indent=2))


if __name__=="__main__":
    main()
