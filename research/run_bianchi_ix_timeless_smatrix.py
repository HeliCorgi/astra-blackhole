from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from bianchi_ix_timeless_smatrix import run_pilot

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True);a=ap.parse_args(argv)
    r=run_pilot();a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"summary.json").write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"matrix":r["matrix"],"reference_mode":r["reference_mode"],
      "negative_control":r["negative_control"],"rows":r["rows"],"window_increment":r["window_increment"]},indent=2))
if __name__=="__main__":main()
