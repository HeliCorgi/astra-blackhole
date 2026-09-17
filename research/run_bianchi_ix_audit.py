from __future__ import annotations
import argparse,json,hashlib,platform
from pathlib import Path
import numpy as np, scipy
from bianchi_ix_model import evolve,kasner_from_state,bkl_step,bkl_orbit,first_bkl_divergence,wall_term_and_gradient
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    initial=np.array([0.,0.,1.,.31]); epochs=np.array([2.,6.,15.,25.])
    sol=evolve(initial,2,25,t_eval=epochs)
    K=[kasner_from_state(s,y) for s,y in zip(epochs,sol.y.T)]
    u=[x['u'] for x in K]
    transitions=[]
    for i in range(2):
        pred=bkl_step(u[i]);err=abs(u[i+1]/pred-1)
        transitions.append({'from_s':float(epochs[i]),'to_s':float(epochs[i+1]),
                            'u_before':u[i],'u_after':u[i+1],'bkl_prediction':pred,
                            'relative_error':float(err)})
    s=np.linspace(2,25,4601); dense=evolve(initial,2,25,t_eval=s,rtol=4e-10,atol=4e-12,max_step=.015)
    y=dense.y; ang=np.unwrap(np.arctan2(y[3],y[2])); rate=np.gradient(ang,s)
    wall=[]
    for ss,yy in zip(s,y.T):
        W,_,_=wall_term_and_gradient(ss,yy[0],yy[1]);wall.append(abs(W)/(yy[2]**2+yy[3]**2))
    wall=np.array(wall)
    peaks=[]
    for i in range(1,len(s)-1):
        if abs(rate[i])>.5 and abs(rate[i])>=abs(rate[i-1]) and abs(rate[i])>abs(rate[i+1]):
            if not peaks or s[i]-peaks[-1]['s']>1:
                peaks.append({'s':float(s[i]),'abs_d_momentum_angle_ds':float(abs(rate[i])),
                              'wall_ratio':float(wall[i])})
    sens=first_bkl_divergence(u[0],1e-12,.1,200)
    pert=initial.copy();pert[1]+=1e-8
    alt=evolve(pert,2,25,t_eval=s,rtol=4e-10,atol=4e-12,max_step=.015)
    sep=np.linalg.norm(alt.y-y,axis=0)
    result={
      'scope':'Vacuum homogeneous Bianchi IX classical minisuperspace; s=-alpha. BKL map is an asymptotic control. No stochastic forcing and no quantum Bianchi IX solver.',
      'initial_state':initial.tolist(),'epoch_samples':epochs.tolist(),'kasner_epochs':[
          {'s':float(ss),**{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in kk.items()}}
          for ss,kk in zip(epochs,K)],
      'bkl_transition_checks':transitions,'detected_bounce_peaks':peaks,
      'continuous_nearby_orbit':{'initial_beta_minus_offset':1e-8,
          'final_phase_space_separation':float(sep[-1]),'max_separation':float(sep.max()),
          'growth_factor':float(sep[-1]/1e-8),
          'interpretation':'Only a finite two-bounce segment; lack of exponential separation here does not refute asymptotic Mixmaster chaos.'},
      'bkl_sensitivity':{'initial_u':u[0],'initial_delta':1e-12,'threshold':.1,
          'first_crossing_step':int(sens[-1][0]) if sens[-1][3]>=.1 else None,
          'gap_at_last_step':float(sens[-1][3]),'rows':[list(map(float,r)) for r in sens]},
      'controls':{'kasner_sum_max_abs':float(max(abs(x['sum']-1) for x in K[:3])),
          'kasner_square_sum_max_abs_in_free_epochs':float(max(abs(x['sum_squares']-1) for x in K[:3])),
          'max_wall_ratio_epoch_samples':float(max(x['wall_ratio'] for x in K[:3]))},
      'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
      'source_sha256':{n:sha(ROOT/n) for n in ['bianchi_ix_model.py','run_bianchi_ix_audit.py']},
      'quantum_status':'Not yet a full quantum Bianchi IX calculation. Existing repository WDW work shows how a quantum geometry state can lack a single classical trajectory, but it is a different reduced model. A same-Hamiltonian 2D Bianchi IX WDW audit remains next.'
    }
    (out/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    curves={'s':s.tolist(),'beta_plus':y[0].tolist(),'beta_minus':y[1].tolist(),
            'p_plus':y[2].tolist(),'p_minus':y[3].tolist(),'momentum_angle':ang.tolist(),
            'wall_ratio':wall.tolist(),'nearby_separation':sep.tolist(),
            'bkl_u':bkl_orbit(u[0],60).tolist()}
    (out/'curves.json').write_text(json.dumps(curves)+'\n')
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'bianchi_ix_results');a=ap.parse_args();r=run(a.out)
    print(json.dumps({'u_transitions':r['bkl_transition_checks'],'bounces':r['detected_bounce_peaks'],'bkl_divergence_step':r['bkl_sensitivity']['first_crossing_step']},indent=2))
if __name__=='__main__':main()
