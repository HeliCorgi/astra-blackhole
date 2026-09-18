from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from bianchi_ix_timeless_bridge import run_audit

def signature(d):
    out=[]
    for r in d["rows"]:
        out += [r["relative_naive_constraint_mismatch"],r["ihbar_dotG_norm"]]
    return out

def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--check-against",type=Path)
    a=ap.parse_args(argv)
    r=run_audit()
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"summary.json").write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n")
    if a.check_against:
        old=json.loads(a.check_against.read_text())
        import numpy as np
        x=np.asarray(signature(r));y=np.asarray(signature(old))
        if x.shape!=y.shape or not np.allclose(x,y,rtol=5e-3,atol=5e-5):
            raise RuntimeError("timeless bridge regression mismatch")
    print(json.dumps({"rows":r["rows"]},ensure_ascii=False,indent=2))
if __name__=="__main__":main()
