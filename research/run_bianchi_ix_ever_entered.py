from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from bianchi_ix_ever_entered import CAPSpec,run_quantum,run_classical_control

def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--v0",type=float,default=.05)
    ap.add_argument("--width",type=float,default=.30)
    ap.add_argument("--maxdim",type=int,default=66)
    ap.add_argument("--skip-classical",action="store_true")
    a=ap.parse_args(argv);t=time.time()
    q=run_quantum(CAPSpec(a.v0,a.width,a.maxdim))
    out={"schema":1,
      "scope":"Finite-clock CAP ever-entered-B audit, conditioned on A at s=4.6, in the existing reduced quantum Bianchi IX square-root model.",
      "class_operators":{"noB":"time-ordered evolution with G_eff=G-i V0 F_B","everB":"U-noB"},
      "region":"B is the union of B+ and B- asymptotic sectors; A is its complement. F_B uses a compact smooth transition around the fixed beta-plane boundary.",
      "quantum":q,
      "classical_control":None if a.skip_classical else run_classical_control(),
      "limitations":["Finite internal-clock construction; not the timeless constraint-commuting Halliwell S-matrix class operator.",
        "Complex-potential strength and smoothing width are regulators; a plateau is required before physical interpretation.",
        "Finite boxes, potential cap, cubic regrids and one square-root ordering are inherited from the parent Bianchi IX audit.",
        "Classical control is a deterministic Wigner ensemble, not an observation or hidden-variable interpretation.",
        "No claim of black-hole singularity resolution or verified quantum gravity."],
      "elapsed_s":time.time()-t}
    a.out.mkdir(parents=True,exist_ok=True);(a.out/"summary.json").write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"spec":q["spec"],"A_weight":q["conditioning"]["A_weight"],"noB":q["noB_diagonal_weight"],
      "everB":q["everB_diagonal_weight"],"offdiag":q["offdiagonal_abs"],"real_offdiag":q["offdiagonal_real"],
      "coherence":q["normalized_coherence"],"additivity_defect":q["additivity_defect"],
      "classical":out["classical_control"],"elapsed_s":out["elapsed_s"]},ensure_ascii=False,indent=2))
if __name__=="__main__":main()
