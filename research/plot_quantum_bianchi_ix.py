from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from bianchi_ix_model import evolve

def ellipse(ax,mean,cov,scale=1.,**kw):
    val,vec=np.linalg.eigh(np.asarray(cov));order=np.argsort(val)[::-1];val=val[order];vec=vec[:,order]
    angle=np.degrees(np.arctan2(vec[1,0],vec[0,0]));e=Ellipse(mean,2*scale*np.sqrt(val[0]),2*scale*np.sqrt(val[1]),angle=angle,fill=False,**kw);ax.add_patch(e)

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--out',type=Path);a=ap.parse_args(argv)
    details=a.results/'details.json';d=json.loads(details.read_text()) if details.exists() else json.loads((a.results/'summary.json').read_text())
    if 'main' not in d: raise ValueError('plot requires regenerated details.json artifact')
    q=d['main']['rows'];e=d['classical_ensemble']['rows']
    s=np.linspace(2,6,1001);cl=evolve(np.array([0.,0.,1.,.31]),2,6,t_eval=s,rtol=3e-10,atol=3e-12,max_step=.01)
    fig,ax=plt.subplots(figsize=(8.7,6.2));ax.plot(cl.y[0],cl.y[1],linestyle='--',label='single classical trajectory')
    qm=np.array([r['mean'] for r in q]);em=np.array([r['mean'] for r in e]);ax.plot(qm[:,0],qm[:,1],marker='o',label='quantum packet mean');ax.plot(em[:,0],em[:,1],marker='s',label='classical Wigner-ensemble mean')
    for target in (4.6,6.0):
        qr=min(q,key=lambda x:abs(x['s']-target));er=min(e,key=lambda x:abs(x['s']-target));ellipse(ax,qr['mean'],qr['covariance'],1.0,linestyle='-');ellipse(ax,er['mean'],er['covariance'],1.0,linestyle=':');ax.annotate(f"s={target:g}",qr['mean'])
    ax.set_xlabel(r'$\beta_+$');ax.set_ylabel(r'$\beta_-$');ax.set_title('Same Bianchi IX model: trajectory, classical ensemble, quantum packet\nFirst wall encounter only; ellipses are 1-sigma covariance diagnostics');ax.grid(True,alpha=.25);ax.legend();fig.tight_layout()
    if a.out:a.out.mkdir(parents=True,exist_ok=True);fig.savefig(a.out/'quantum_bianchi_ix_compare.png',dpi=170)
    return fig
if __name__=='__main__':main()
