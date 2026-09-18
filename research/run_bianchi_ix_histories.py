from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from bianchi_ix_histories import run_a_conditioned_histories

def cpair(z):
    z=complex(z);return {"re":float(z.real),"im":float(z.imag)}
def matrix_json(M): return [[cpair(x) for x in row] for row in np.asarray(M)]

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--main-maxdim",type=int,default=54);ap.add_argument("--control-maxdim",type=int,default=66)
    a=ap.parse_args(argv);t0=time.time();r=run_a_conditioned_histories(a.main_maxdim,a.control_maxdim)
    ranked=sorted(r["rows"],key=lambda x:x["conditional_diagonal_weight"],reverse=True)
    out={
      "schema":2,
      "scope":"Nine A-conditioned three-time histories in the same reduced quantum Bianchi IX square-root quantization. Second-wall branches are propagated with one common block-Krylov projected operator at every step.",
      "class_operator":"C_{A,j,k}=P_k(s3) U(s3,s2) P_j(s2) U(s2,s1) P_A(s1) U(s1,s0)",
      "decoherence_functional":"D(alpha,beta)=<psi_beta|psi_alpha>; conditional D is normalized by the coherent final A-branch norm after the same linear staged regrids.",
      "history_times":{"s0":r["times"].s0,"s1":r["times"].s1,"s2":r["times"].s2,"s3":r["times"].s3},
      "first_wall":{"A_weight":r["first_a_weight"],"non_A_weight":r["first_non_a_weight"],
        "conditioning_note":"Primary family conditions on A. The remaining first-wall non-A mass is reported but not expanded into all 27 histories."},
      "matrix_order":[f"A->{j}->{k}" for j,k in r["keys"]],
      "decoherence_matrix_conditional":matrix_json(r["D_conditional"]),
      "decoherence_measures":r["measures"],"history_diagonal_weights":r["rows"],
      "ranked_history_diagonal_weights":ranked,"third_wall_coarse_graining":r["coarse"],
      "closure":{**r["closure"],"sum_all_D_normalized":cpair(r["closure"]["sum_all_D_normalized"])},
      "block_krylov_control":r["control"],"elapsed_s":time.time()-t0,
      "rejected_prototype":{"reason":"Independent per-branch Lanczos approximations did not define one common linear numerical propagator.","branch_sum_minus_direct_relative_norm":0.023884798136052206,"status":"rejected; not used for scientific interpretation."},
      "interpretation_rules":[
        "Diagonal entries are history weights, not ordinary probabilities unless the chosen family decoheres to the stated tolerance.",
        "Coarse/fine additivity defects and off-diagonal D entries measure interference between alternative second-wall histories ending in the same third-wall sector.",
        "Different final sectors are orthogonal by construction; nontrivial consistency concerns earlier alternatives with a common final sector.",
        "The first A result is conditioned because its one-time weight is about 98%; omitted first-wall non-A outcomes are reported separately.",
        "Cubic regrids are numerical and intentionally unnormalized branch maps. They are not physical measurements.",
        "This finite-clock construction is not the invariant timeless class-operator construction for a Wheeler-DeWitt constraint."
      ]}
    a.out.mkdir(parents=True,exist_ok=True);(a.out/"summary.json").write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"first_A_weight":out["first_wall"]["A_weight"],"measures":out["decoherence_measures"],
      "coarse":out["third_wall_coarse_graining"],"top_histories":ranked[:4],"closure":out["closure"],
      "control":out["block_krylov_control"],"elapsed_s":out["elapsed_s"]},ensure_ascii=False,indent=2))
if __name__=="__main__":main()
