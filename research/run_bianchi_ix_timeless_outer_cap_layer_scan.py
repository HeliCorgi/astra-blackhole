from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))

from bianchi_ix_timeless_smatrix import WDWGridSpec,build_constraint,edge_mass,low_energy_edge_min_packet
from bianchi_ix_timeless_outer_cap import outer_cap_profile,t_action,weighted_mass

def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args(argv)

    grid_spec=WDWGridSpec()
    g=build_constraint(grid_spec)
    C=g["C"]
    _,psi,packet=low_energy_edge_min_packet(g,20)
    etas=(0.10,0.05,0.025)
    v0s=(0.025,0.05)
    gammas=(0.10,0.20)
    layers=(1,2,3)

    states={}
    rows=[]
    profile_weights={}
    for layer in layers:
        profile=outer_cap_profile(g["shape"],layer).ravel()
        profile_weights[str(layer)]=weighted_mass(profile,psi)
        for gamma in gammas:
            outer=gamma*profile
            for v0 in v0s:
                vd=v0*np.asarray(g["F_B"],float).ravel()
                for eta in etas:
                    t,response=t_action(C,psi,0.0,eta,outer,vd)
                    states[(layer,gamma,v0,eta)]=t
                    rows.append({
                        "layer_cells":layer,
                        "outer_strength":gamma,
                        "v0":v0,
                        "eta":eta,
                        "t_norm":float(np.linalg.norm(t)),
                        "response_edge_mass_one_cell":edge_mass(g,response,1),
                        "response_outer_layer_weight":weighted_mass(profile,response),
                    })

    eta_change=[]
    for layer in layers:
        for gamma in gammas:
            for v0 in v0s:
                for hi,lo in zip(etas[:-1],etas[1:]):
                    x=states[(layer,gamma,v0,hi)]
                    y=states[(layer,gamma,v0,lo)]
                    eta_change.append({
                        "layer_cells":layer,"outer_strength":gamma,"v0":v0,
                        "from_eta":hi,"to_eta":lo,
                        "relative_t_change":float(np.linalg.norm(y-x)/max(np.linalg.norm(y),1e-30)),
                    })

    layer_change=[]
    for gamma in gammas:
        for v0 in v0s:
            for eta in etas:
                for l0,l1 in zip(layers[:-1],layers[1:]):
                    x=states[(l0,gamma,v0,eta)]
                    y=states[(l1,gamma,v0,eta)]
                    layer_change.append({
                        "outer_strength":gamma,"v0":v0,"eta":eta,
                        "from_layer_cells":l0,"to_layer_cells":l1,
                        "relative_t_change":float(np.linalg.norm(y-x)/max(np.linalg.norm(y),1e-30)),
                    })

    result={
        "schema":1,
        "method":"Outer-boundary CAP layer-width control; compares 1/2/3-cell numerical layers at fixed strengths.",
        "grid":{"dimension":int(C.shape[0]),"shape":g["shape"]},
        "scan":{"etas":etas,"v0s":v0s,"outer_strengths":gammas,"layer_cells":layers},
        "reference_packet":{**packet,"profile_weights":profile_weights},
        "rows":rows,
        "eta_change":eta_change,
        "layer_change":layer_change,
        "decision_rule":"Improvement is not accepted if it requires a broad layer with large reference-state overlap or if T-action changes materially with layer width.",
        "limitations":[
            "Numerical outer CAP only; not PML or exterior complex scaling.",
            "No independent reflection coefficient is measured.",
            "No physical inner product, decoherence functional, or history probability is computed.",
        ],
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"layer_scan.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
