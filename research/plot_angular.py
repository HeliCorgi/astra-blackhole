"""Plot an angular-state audit. This is NOT observational data."""
from pathlib import Path
import argparse,json
import matplotlib.pyplot as plt


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results',type=Path,default=Path(__file__).parent/'angular_results')
    args=ap.parse_args();data=json.loads((args.results/'curves.json').read_text())
    for mass in (1,0):
        fig,ax=plt.subplots(figsize=(8.5,5.0))
        for ell,style in ((4,'-'),(6,'--'),(8,'-.')):
            ax.plot(data['time'][1:],data[f'm{mass}_ell{ell}_separation_percent'][1:],style,label=f'Angular P{ell} pair')
        ax.set_yscale('log');ax.set_xlabel('Proper time (model units)')
        ax.set_ylabel('Symmetric relative separation of K (%)')
        ax.set_title(f'Omitted angular structure, particle mass = {mass}\nClassical homogeneous model; not forecast error or observations')
        ax.grid(True,which='both',alpha=.25);ax.legend();fig.tight_layout()
        fig.savefig(args.results/f'angular_mass{mass}.png',dpi=160)
    return args.results

if __name__=='__main__':main()
