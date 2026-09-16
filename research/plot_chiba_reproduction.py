"""Plot newly computed observables, not digitized author data."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def main(argv=None):
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--out',type=Path);a=p.parse_args(argv)
    curves=json.loads((a.results/'curves.json').read_text());out=a.out or a.results;out.mkdir(parents=True,exist_ok=True)
    fig,ax=plt.subplots(figsize=(8.7,5.5))
    for kappa,rows in curves.items():
        use=[r for r in rows if r['mean_U']<=.02]
        ax.plot([r['mean_V'] for r in use],[r['mean_U'] for r in use],label=f'Published prescription, kappa={kappa}')
    v=np.linspace(0,.32,120);ax.plot(v,v-.3,'--',label='Classical trajectory')
    ax.axhline(0,ls=':',lw=1)
    ax.set(xlabel='V = Cartesian clock + signed mean X',ylabel='U = Cartesian clock - signed mean X',
           title='Same published WDW initial data and Klein-Gordon prescription\nRecomputed equations, not a point-by-point match to author data',ylim=(-.33,.025))
    ax.grid(True,alpha=.25);ax.legend(fontsize=9);fig.tight_layout();fig.savefig(out/'chiba_expectation.png',dpi=170)
    return fig
if __name__=='__main__':main()
