from __future__ import annotations
import argparse,json
from pathlib import Path
import matplotlib.pyplot as plt

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--out',type=Path);a=ap.parse_args(argv)
    d=json.loads((a.results/'curves.json').read_text());fig,ax=plt.subplots(figsize=(8.8,5.2))
    for row in d['ks']:
        label=f"KS perturb {row['perturbation']}"
        ax.plot(d['ks_times'],[abs(x) for x in row['lyapunov']],label=label)
    ax.axhline(2.373138220831251,linestyle='--',label='Gauss/BKL invariant-measure Lyapunov')
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Evolution time / map-control scale');ax.set_ylabel('|finite-time Lyapunov diagnostic|')
    ax.set_title('Current KS reduction is integrable; BKL control is chaotic\nDifferent systems, same diagnostic family; not a black-hole observation')
    ax.grid(True,which='both',alpha=.25);ax.legend(fontsize=8);fig.tight_layout()
    if a.out:
        a.out.mkdir(parents=True,exist_ok=True);fig.savefig(a.out/'trajectory_chaos.png',dpi=170)
    return fig
if __name__=='__main__':main()
