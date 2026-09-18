from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from bianchi_ix_timeless_resolvent import run_pilot

def _signature(d):
    r=d["reference_packet"]
    v=[
        r["edge_mass_one_cell"],r["energy_mean"],r["energy_spread"],
        r["capped_region_mass"],r["constraint_zero_residual"],
        d["negative_control"]["v0_zero_t_norm"],
    ]
    for x in d["rows"]:
        v += [
            x["dyson_relative"],x["free_resolvent_norm"],x["t_norm"],
            x["t_to_born_norm_ratio"],x["t_shell_fraction"],
            x["response_norm"],x["response_edge_mass_one_cell"],
        ]
    for x in d["eta_change"]:
        v.append(x["relative_t_change"])
    return np.asarray(v,float)

def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--check-against",type=Path)
    a=ap.parse_args(argv)
    r=run_pilot()
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"summary.json").write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n")
    if a.check_against:
        old=json.loads(a.check_against.read_text())
        x=_signature(r);y=_signature(old)
        if x.shape!=y.shape or not np.allclose(x,y,rtol=3e-5,atol=3e-8):
            raise RuntimeError("timeless resolvent regression mismatch")
    print(json.dumps({
        "grid":r["grid"],
        "reference_packet":r["reference_packet"],
        "negative_control":r["negative_control"],
        "rows":r["rows"],
        "eta_change":r["eta_change"],
    },indent=2))
if __name__=="__main__":
    main()
