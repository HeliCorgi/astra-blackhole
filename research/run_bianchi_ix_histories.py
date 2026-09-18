from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from bianchi_ix_histories import run_a_conditioned_histories


def cpair(z):
    z=complex(z)
    return {"re":float(z.real),"im":float(z.imag)}


def matrix_json(M):
    return [[cpair(x) for x in row] for row in np.asarray(M)]


def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--control-krylov",type=int,default=60)
    a=ap.parse_args(argv)
    t0=time.time()
    r=run_a_conditioned_histories(a.control_krylov)

    # Sort only for presentation; matrix order remains explicitly recorded.
    ranked=sorted(r["rows"],key=lambda x:x["conditional_diagonal_weight"],reverse=True)
    out={
        "schema":1,
        "scope":"Nine A-conditioned three-time histories in the same reduced quantum Bianchi IX square-root quantization. Projectors are A/B+/B- beta-space wall sectors at s=4.6, 10.385, 32.46.",
        "class_operator":"C_{A,j,k}=P_k(s3) U(s3,s2) P_j(s2) U(s2,s1) P_A(s1) U(s1,s0)",
        "decoherence_functional":"D(alpha,beta)=<psi_beta|psi_alpha>; reported conditional matrix is normalized by the coherent final A-branch norm after the same linear staged regrids.",
        "history_times":{"s0":r["times"].s0,"s1":r["times"].s1,"s2":r["times"].s2,"s3":r["times"].s3},
        "first_wall":{"A_weight":r["first_a_weight"],"non_A_weight":r["first_non_a_weight"],"conditioning_note":"The nine-history primary family conditions on A at the first wall. Non-A first-wall histories are not expanded into the full 27-history family in this audit."},
        "matrix_order":[f"A->{j}->{k}" for j,k in r["keys"]],
        "decoherence_matrix_conditional":matrix_json(r["D_conditional"]),
        "decoherence_measures":r["measures"],
        "history_diagonal_weights":r["rows"],
        "ranked_history_diagonal_weights":ranked,
        "third_wall_coarse_graining":r["coarse"],
        "closure":{**r["closure"],"sum_all_D_normalized":cpair(r["closure"]["sum_all_D_normalized"])},
        "krylov_control":r["control"],
        "elapsed_s":time.time()-t0,
        "interpretation_rules":[
            "Diagonal entries are history weights. They are probabilities only when the chosen history family decoheres to the required accuracy.",
            "Nonzero additivity defects between a coarse final event and the sum of fine-history diagonals are direct interference diagnostics.",
            "Different final-sector histories are orthogonal by construction; the nontrivial consistency test is between different earlier histories ending in the same final sector.",
            "The first A outcome is conditioned because the published first-wall run assigns about 98% one-time mass to A; the omitted non-A first outcomes are reported separately.",
            "Staged cubic regrids are numerical, linear and intentionally unnormalized in the history calculation. They are not physical measurements.",
            "This is not a unique timeless Wheeler-DeWitt class operator, a consistent-histories proof for quantum gravity, or a black-hole observation."
        ]
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"summary.json").write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({
        "first_A_weight":out["first_wall"]["A_weight"],
        "measures":out["decoherence_measures"],
        "coarse":out["third_wall_coarse_graining"],
        "top_histories":out["ranked_history_diagonal_weights"][:4],
        "elapsed_s":out["elapsed_s"]
    },ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
