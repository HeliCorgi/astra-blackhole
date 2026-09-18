from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from bianchi_ix_timeless_smatrix import run_pilot
def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True);ap.add_argument("--check-against",type=Path);a=ap.parse_args(argv)
    r=run_pilot();a.out.mkdir(parents=True,exist_ok=True);(a.out/"summary.json").write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n")
    if a.check_against:
        old=json.loads(a.check_against.read_text())
        import numpy as np
        def sig(d):
            v=[d["reference_packet"]["edge_mass_one_cell"],d["reference_packet"]["energy_spread"],d["negative_control"]["v0_zero_max_commutator_relative"]]
            for x in d["rows"]:
                v += [x["norm2"],x["constraint_energy_spread"],x["commutator_absolute"],x["commutator_relative"]]
            for x in d["window_increment"]:
                v.append(x["relative_state_change"])
            return np.asarray(v,float)
        x=sig(r);y=sig(old)
        if x.shape!=y.shape or not np.allclose(x,y,rtol=3e-5,atol=3e-8):
            raise RuntimeError("timeless S-matrix regression mismatch")
    print(json.dumps({"matrix":r["matrix"],"reference_packet":r["reference_packet"],
      "negative_control":r["negative_control"],"rows":r["rows"],"window_increment":r["window_increment"]},indent=2))
if __name__=="__main__":main()
