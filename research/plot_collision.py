"""Plot recorded model comparisons; not observational data."""
from pathlib import Path
import argparse,json
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results',type=Path,default=Path(__file__).parent/'collision_results')
    args=ap.parse_args();r=json.loads((args.results/'results.json').read_text())
    rows=[x for x in r['cases'] if x['ell']==4 and x['rate']>0]
    fig,ax=plt.subplots(figsize=(8.5,5.2))
    ax.loglog([x['rate'] for x in rows],[x['pair_K_gap_percent'] for x in rows], 'o-',label='Difference between matched states (symmetric %)')
    ax.loglog([x['rate'] for x in rows],[x['max_fluid_K_error_percent'] for x in rows], 's--',label='Departure from perfect fluid (fluid denominator, %)')
    ax.set_xlabel('Constant relaxation rate / inverse model-time unit')
    ax.set_ylabel('Curvature discrepancy at t = 0.45 (%)')
    ax.set_title('Initial-state ambiguity and fluid-closure discrepancy are different\nMassless classical RTA benchmark; P4 matched-state pair')
    ax.grid(True,which='both',alpha=.25);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(args.results/'collision_comparison.png',dpi=180)
    return fig
if __name__=='__main__':main()
