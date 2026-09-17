from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--out',type=Path);a=ap.parse_args(argv)
    d=json.loads((a.results/'curves.json').read_text());s=np.array(d['s'])
    fig,ax=plt.subplots(figsize=(9,5.4))
    ax.plot(d['beta_plus'],d['beta_minus'],label='continuous Bianchi IX anisotropy path')
    ax.set_xlabel('beta+');ax.set_ylabel('beta-');ax.set_title('Bianchi IX: continuous wall bounces before the asymptotic BKL map\nClassical homogeneous model; not a black-hole observation')
    ax.grid(True,alpha=.25);ax.legend();fig.tight_layout()
    if a.out:a.out.mkdir(parents=True,exist_ok=True);fig.savefig(a.out/'bianchi_ix_path.png',dpi=170)
    return fig
if __name__=='__main__':main()
