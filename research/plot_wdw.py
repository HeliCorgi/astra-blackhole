"""Quantum/classical radius comparison. Copyright 2026 HeliCorgi. Apache-2.0."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib.pyplot as plt

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);a=ap.parse_args(argv)
    d=json.loads((a.results/'curves.json').read_text());t=np.asarray(d['time'])
    fig,ax=plt.subplots(figsize=(9,5.2))
    ax.plot(t,np.exp(2-np.asarray(d['point_mean_x'])-t),'--',label='Classical point trajectory')
    for name in ['h0.4','h0.2','h0.1']:
        ax.plot(t,d['cases'][name]['geometric_mean_radius_ratio'],label='Quantum '+name)
    ax.set_yscale('log');ax.set_xlabel('Intrinsic clock T = log(a), NOT proper time')
    ax.set_ylabel('Geometric mean areal radius / reference initial radius')
    ax.set_title('Restricted WDW geometry: reflection in x is NOT a radius bounce\nVacuum minisuperspace, chosen clock and inner product; no real-world validation')
    ax.grid(True,which='both',alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(a.results/'wdw_radius.png',dpi=170)
    return fig
if __name__=='__main__':main();plt.close('all')
