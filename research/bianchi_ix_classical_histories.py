"""Classical three-time wall-sector histories for the same Bianchi IX Wigner control.

Each sample follows one deterministic classical Hamilton trajectory, so joint
A/B+/B- labels at fixed times define ordinary empirical frequencies.  This is
kept separate from quantum branch diagonal weights, which require decoherence
before they can be interpreted as probabilities.
"""
from __future__ import annotations
import numpy as np
from quantum_bianchi_ix_model import PacketSpec, initial_wigner_ensemble, _rhs_vec
from quantum_bianchi_ix_long_model import positive_wall_terms, WALL_NAMES


def history_labels(packet=PacketSpec(hbar=.2,sigma=.3),times=(4.6,10.385,32.46),
                   n=4096,seed=20260917,ds=.01):
    times=np.asarray(times,float)
    if np.any(np.diff(times)<=0) or times[0]<=2:
        raise ValueError("ordered history times after s=2 required")
    Y=initial_wigner_ensemble(packet,n,seed);s=2.0;labels=[]
    for target in times:
        while s < target-1e-12:
            h=min(float(ds),float(target-s))
            k1=_rhs_vec(s,Y);k2=_rhs_vec(s+h/2,Y+h*k1/2)
            k3=_rhs_vec(s+h/2,Y+h*k2/2);k4=_rhs_vec(s+h,Y+h*k3)
            Y += h*(k1+2*k2+2*k3+k4)/6;s+=h
        labels.append(np.argmax(positive_wall_terms(target,Y[0],Y[1]),axis=0))
    return np.asarray(labels,int)


def summarize(labels):
    labels=np.asarray(labels,int)
    if labels.ndim!=2 or labels.shape[0]!=3:
        raise ValueError("expected 3 x N labels")
    n=labels.shape[1]
    full=[]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                count=int(np.sum((labels[0]==i)&(labels[1]==j)&(labels[2]==k)))
                full.append({"history":f"{WALL_NAMES[i]}->{WALL_NAMES[j]}->{WALL_NAMES[k]}",
                             "count":count,"frequency":count/n})
    A=(labels[0]==0);na=int(A.sum())
    cond=[]
    for j in range(3):
        for k in range(3):
            count=int(np.sum(A&(labels[1]==j)&(labels[2]==k)))
            cond.append({"history":f"A->{WALL_NAMES[j]}->{WALL_NAMES[k]}",
                         "count":count,
                         "conditional_frequency":count/na if na else float("nan"),
                         "unconditional_frequency":count/n})
    return {"sample_count":int(n),"first_A_count":na,"first_A_frequency":na/n,
            "full_27":full,"A_conditioned_9":cond}


def run_control(ds=.01,refined_ds=.005):
    a=history_labels(ds=ds);b=history_labels(ds=refined_ds)
    sa=summarize(a);sb=summarize(b)
    pa=np.array([x["conditional_frequency"] for x in sa["A_conditioned_9"]])
    pb=np.array([x["conditional_frequency"] for x in sb["A_conditioned_9"]])
    return {"schema":1,"times":[4.6,10.385,32.46],"main_ds":ds,"refined_ds":refined_ds,
            "main":sa,"refined_max_A_conditioned_abs_difference":float(np.max(abs(pa-pb))),
            "label_change_fraction":float(np.mean(np.any(a!=b,axis=0))),
            "scope":"Same positive Gaussian Wigner ensemble as the quantum Bianchi IX audit; deterministic classical Hamilton trajectories.",
            "limitations":["4096 fixed Sobol samples; empirical ensemble control, not an observation.",
                           "Sector labels are fixed-time beta-space labels, not continuous first-hit times.",
                           "Some larger Sobol ensembles can include samples that leave the chosen positive reduced branch; this audit intentionally preserves the existing 4096-sample control."]}


if __name__=="__main__":
    import json
    print(json.dumps(run_control(),indent=2))
