from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--out',type=Path);a=ap.parse_args(argv)
    d=json.loads((a.results/'summary.json').read_text());rows=d['selected_comparisons'];s=np.array([x['s'] for x in rows])
    fig,ax=plt.subplots(figsize=(8.8,5.2));ax.plot(s,[x['quantum_to_single_trajectory'] for x in rows],marker='o',label='quantum centroid vs single classical orbit');ax.plot(s,[x['quantum_to_classical_ensemble_mean'] for x in rows],marker='s',label='quantum centroid vs classical ensemble mean')
    for b in d['central_bounces']:ax.axvline(b['s'],linestyle='--',alpha=.55,label=f"classical wall {b['bounce_number']} ({b['central_wall']})")
    ax.set_xlabel('s=-alpha');ax.set_ylabel('distance in anisotropy plane');ax.set_title('Same Bianchi IX model through the third classical wall encounter\nCentroid departure is not the same thing as a quantum-only effect');ax.grid(True,alpha=.25);ax.legend(fontsize=8);fig.tight_layout()
    if a.out:a.out.mkdir(parents=True,exist_ok=True);fig.savefig(a.out/'quantum_bianchi_ix_long_distance.png',dpi=170)
    return fig
if __name__=='__main__':main()
