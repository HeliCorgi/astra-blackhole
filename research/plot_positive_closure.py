"""Plot paired closure errors; no change in state count or time horizon.
Copyright 2026 HeliCorgi. Apache-2.0.
"""
from pathlib import Path
import argparse,json
import matplotlib.pyplot as plt

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=Path(__file__).resolve().parent/'positive_closure_results')
    p.add_argument('--out',type=Path)
    a=p.parse_args(argv);out=a.out or a.results;out.mkdir(parents=True,exist_ok=True)
    data=json.loads((a.results/'summary.json').read_text())
    fig,ax=plt.subplots(figsize=(10,5.9))
    for method,label,linestyle in [
        ('B_zero','Density + pressure difference: zero tail','--'),
        ('B_positive','Same state: positive angular closure','-'),
        ('C_zero','Also retain c4: zero tail',':'),
        ('C_positive','Same state + c4: positive angular closure','-.')]:
        sub=[x for x in data['aggregates'] if x['split']=='heldout' and x['method']==method]
        ax.plot([x['rate'] for x in sub],[100*x['worst_Pi_error'] for x in sub],
                linestyle,marker='o',label=label)
    ax.axhline(1.,linestyle=':',linewidth=1,label='Fixed 1% pressure criterion')
    ax.set_xscale('symlog',linthresh=1);ax.set_yscale('log')
    ax.set_xlim(left=-.3)
    ax.set_xlabel('Relaxation rate (model units; 0 = collisionless)')
    ax.set_ylabel('Maximum |Pi model - Pi reference| / rho reference (%)')
    ax.set_title('Same 60 previously evaluated geometry cases, interval [0, 0.30]\n'
                 'Classical benchmark: positivity is not a guarantee of greater accuracy')
    ax.grid(True,which='both',alpha=.25);ax.legend(fontsize=9)
    fig.tight_layout();fig.savefig(out/'positive_closure_comparison.png',dpi=170)
    return fig
if __name__=='__main__':main()
