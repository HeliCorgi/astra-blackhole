from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from bianchi_ix_classical_histories import run_control

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True);a=ap.parse_args(argv)
    r=run_control();a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"summary.json").write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"first_A_frequency":r["main"]["first_A_frequency"],
      "A_conditioned_9":r["main"]["A_conditioned_9"],
      "refined_max_abs":r["refined_max_A_conditioned_abs_difference"],
      "label_change_fraction":r["label_change_fraction"]},ensure_ascii=False,indent=2))
if __name__=="__main__":main()
