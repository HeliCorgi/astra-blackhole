"""Plot the fixed held-out closure comparison; no observational data.
Copyright 2026 HeliCorgi. Apache-2.0.
"""
import argparse,json
from pathlib import Path
import matplotlib.pyplot as plt

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=Path('artifacts/closure'))
    args=parser.parse_args(argv)
    result=json.loads((args.results/'summary.json').read_text())
    fig,ax=plt.subplots(figsize=(8.5,5.2))
    for model,label in [('A','A: perfect fluid'),('B','B: density + pressure difference'),('C','C: also retain c4')]:
        rows=[r for r in result['aggregates'] if r['split']=='heldout' and r['model']==model]
        ax.plot([r['rate'] for r in rows],[100*r['worst_Pi_error'] for r in rows],marker='o',label=label)
    ax.axhline(1,linestyle=':',label='1% pressure-error criterion')
    ax.set_xscale('symlog',linthresh=1);ax.set_yscale('log')
    ax.set_xlabel('Relaxation rate (model units; 0 = collisionless)')
    ax.set_ylabel('Maximum |Pi_model - Pi_ref| / rho_ref (%)')
    ax.set_title('Same interval [0, 0.30]; 60 held-out classical cases\nAccuracy alone does not certify a positive reconstructed distribution')
    ax.grid(True,which='both',alpha=.25);ax.legend(fontsize=9)
    fig.tight_layout();fig.savefig(args.results/'closure_pressure_errors.png',dpi=170)
    return fig
if __name__=='__main__':main()
