"""Plot measured cutoff dependence, not a finite continuum curvature claim.
Copyright 2026 HeliCorgi. Apache-2.0.
"""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--out',type=Path)
    a=ap.parse_args(argv);out=a.out or a.results;out.mkdir(parents=True,exist_ok=True)
    rows=json.loads((a.results/'cutoffs.json').read_text())
    rows=[r for r in rows if r['grid']=='R256' and r['packet']=='h0.2' and r['T']==2]
    fig,ax=plt.subplots(figsize=(9,5.3))
    for q in (2,4,6):ax.plot([r['cutoff'] for r in rows],[r[f'log10_partial_J{q}'] for r in rows],marker='o',label=f'q={q}')
    ax.set_xlabel('Observable cutoff L in x (solver boundary fixed at x=256)')
    ax.set_ylabel('log10 of partial expectation of (r_ref/r)^q')
    ax.set_title('Including more of the tail prevents a finite plateau\nSame reduced WDW state, h=0.2, T=2; not full 4D curvature')
    ax.grid(True,alpha=.25);ax.legend();fig.tight_layout();fig.savefig(out/'wdw_tail_moments.png',dpi=170)
    data=json.loads((a.results/'plot_data.json').read_text())['R256']['h0.2']
    fig2,ax2=plt.subplots(figsize=(9,5.3));ax2.loglog(data['x'],data['density'],label='Numerical wavefunction density')
    ax2.loglog(data['x'],data['asymptotic_density'],'--',label='Threshold prediction: C / (x + shift)^6')
    ax2.set_xlabel('x (larger x means smaller areal radius at fixed clock)');ax2.set_ylabel('Probability density in x')
    ax2.set_title('Non-Gaussian tail: coefficient computed from the initial state\nNo fit to the future tail; restricted WDW representation')
    ax2.grid(True,which='both',alpha=.25);ax2.legend();fig2.tight_layout();fig2.savefig(out/'wdw_tail_density.png',dpi=170)
    return fig,fig2
if __name__=='__main__':main()
