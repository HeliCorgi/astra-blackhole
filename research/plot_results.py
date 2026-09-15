"""Plot saved results without prescribing colors or styles."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def main():
    folder=Path(__file__).resolve().parent/'results'
    d=np.load(folder/'trajectories.npz')
    fig,ax=plt.subplots(figsize=(9,5.4))
    labels={1:'Same density, pressure and initial geometry',
            2:'Also match the fourth-order momentum moment',
            3:'Also match the sixth-order momentum moment'}
    for level in (1,2,3):
        gap=d[f'level{level}_gap_percent'];good=(d['time']>=.015)&(gap>1e-10)
        ax.plot(d['time'][good],gap[good],label=labels[level],linewidth=2)
    ax.set_yscale('log')
    ax.set_xlabel('Elapsed proper time (model units)')
    ax.set_ylabel('Symmetric separation in curvature (%)')
    ax.set_title('What the initial averaged state fails to distinguish\n'
                 'Classical Einstein–Vlasov benchmark — not observational data')
    ax.grid(True,which='both',alpha=.25)
    ax.legend(fontsize=9,loc='lower right')
    fig.text(.5,.01,'Each curve compares a DIFFERENT matched pair; this is not a predictor-accuracy comparison.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.04,1,1))
    fig.savefig(folder/'omitted_moment_separation.png',dpi=180)
    return fig
if __name__=='__main__':
    main();plt.show()
