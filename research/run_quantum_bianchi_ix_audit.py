from __future__ import annotations
import argparse,json,hashlib,platform,time
from pathlib import Path
import numpy as np, scipy
from bianchi_ix_model import evolve
from quantum_bianchi_ix_model import QuantumBianchiIX,PacketSpec,propagate_quantum,propagate_classical_ensemble
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def by_s(rows):return {round(float(x['s']),6):x for x in rows}
def mean_diff(a,b):return float(np.linalg.norm(np.asarray(a)-np.asarray(b)))

def qcase(name,**kw):
    grid=QuantumBianchiIX(hbar=.2,bp_bounds=(-5,5),bm_bounds=(-8,5),nx=kw.pop('nx',100),ny=kw.pop('ny',130),potential_cap=kw.pop('cap',30))
    packet=PacketSpec(hbar=.2,sigma=.3)
    s1=kw.pop('s1',6.); save=kw.pop('save_times',(2,2.5,3,3.5,4,4.3,4.6,5,5.5,6) if s1>=6 else (2,2.5,3,3.5,4,4.3,4.6))
    t0=time.time();run=propagate_quantum(grid,packet,2,s1,kw.pop('ds',.05),kw.pop('krylov',48),save,keep_density=kw.pop('density',False))
    if kw:raise ValueError(kw)
    return {'name':name,'grid':{'nx':grid.nx,'ny':grid.ny,'dx':grid.dx,'dy':grid.dy,'cap':grid.cap},'elapsed_s':time.time()-t0,**run}

def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);packet=PacketSpec(hbar=.2,sigma=.3);times=np.array([2.,2.5,3.,3.5,4.,4.3,4.6,5.,5.5,6.])
    main=qcase('main',density=True,krylov=48,ds=.05)
    checks=[qcase('krylov56',krylov=56,ds=.05),
            qcase('cap50_first_bounce',krylov=48,ds=.05,cap=50,s1=4.6),
            qcase('dt04_first_bounce',krylov=48,ds=.04,s1=4.6),
            qcase('coarse_grid_first_bounce',nx=80,ny=104,krylov=48,ds=.05,s1=4.6)]
    ens=propagate_classical_ensemble(packet,2,6,.004,times,n=4096,seed=20260917)
    central=evolve(np.array([0.,0.,1.,.31]),2,6,t_eval=times,rtol=2e-11,atol=2e-13,max_step=.005)
    central_rows=[{'s':float(s),'mean':central.y[:2,i].tolist()} for i,s in enumerate(times)]
    qm=by_s(main['rows']); em=by_s(ens['rows']); cm=by_s(central_rows)
    comparisons=[]
    for s in times:
        q=qm[round(float(s),6)];e=em[round(float(s),6)];c=cm[round(float(s),6)]
        comparisons.append({'s':float(s),'quantum_to_single_trajectory':mean_diff(q['mean'],c['mean']),
                            'quantum_to_classical_ensemble_mean':mean_diff(q['mean'],e['mean']),
                            'quantum_rms_width':q['rms_width'],'classical_ensemble_rms_width':e['rms_width'],
                            'width_ratio_quantum_over_ensemble':float(q['rms_width']/e['rms_width'])})
    ref=by_s(checks[0]['rows']);numerical=[]
    for case in [main,*checks[1:]]:
        rr=by_s(case['rows']);d=[];common=sorted(set(rr).intersection(ref))
        for ss in common:
            a=rr[ss];b=ref[ss]
            d.append({'s':float(ss),'mean_distance':mean_diff(a['mean'],b['mean']),
                      'rms_width_abs_difference':abs(a['rms_width']-b['rms_width']),
                      'edge_mass':a['edge_mass'],'capped_region_mass':a['capped_region_mass']})
        numerical.append({'case':case['name'],'rows':d,'max_mean_distance':max(x['mean_distance'] for x in d),
                          'max_rms_width_abs_difference':max(x['rms_width_abs_difference'] for x in d)})
    densities={k:v.tolist() for k,v in main['densities'].items()}
    result={'scope':'Vacuum homogeneous Bianchi IX, same reduced classical principal symbol and a chosen direct square-root quantization. First resolved wall encounter only; no long-time quantum-chaos claim.',
      'quantization':{'generator':'G(s)=-sqrt(-hbar^2 Delta_beta + W(s,beta))','hbar':.2,'sigma':.3,
        'initial_center':[0,0],'initial_momentum':[1,.31],'ordering':'flat beta-plane finite-difference Laplacian plus multiplicative W',
        'boundary':'finite Dirichlet box','potential_cap':'numerical regulator only; accepted runs require small capped-region mass'},
      'main':{k:v for k,v in main.items() if k not in ('final_state','densities')},'classical_ensemble':{k:v for k,v in ens.items() if k!='final'},
      'central_trajectory':central_rows,'comparisons':comparisons,'numerical_checks':numerical,
      'headline':{'at_s_4p6':next(x for x in comparisons if x['s']==4.6),'at_s_6':next(x for x in comparisons if x['s']==6.0),
        'main_max_edge_mass':max(x['edge_mass'] for x in main['rows']),'main_max_capped_region_mass':max(x['capped_region_mass'] for x in main['rows'])},
      'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
      'sources':{n:sha(ROOT/n) for n in ['bianchi_ix_model.py','quantum_bianchi_ix_model.py','run_quantum_bianchi_ix_audit.py']},
      'limitations':['This is a quantization choice, not a unique Wheeler-DeWitt physical Hilbert space.','Only s=2..6, covering the first wall reflection, is accepted.','The classical Wigner ensemble is a control for finite initial spread, not a hidden-variable interpretation.','Grid Shannon entropy is diagnostic only.','No claim that a wave-packet centroid is a physical trajectory.']}
    compact={'scope':result['scope'],'quantization':result['quantization'],'headline':result['headline'],
      'selected_comparisons':[x for x in comparisons if x['s'] in (2.0,4.6,6.0)],
      'numerical_check_maxima':[{'case':x['case'],'max_mean_distance':x['max_mean_distance'],'max_rms_width_abs_difference':x['max_rms_width_abs_difference']} for x in numerical],
      'environment':result['environment'],'sources':result['sources'],'limitations':result['limitations']}
    (out/'summary.json').write_text(json.dumps(compact,ensure_ascii=False,indent=2)+'\n')
    (out/'details.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    (out/'densities.json').write_text(json.dumps({'bp':QuantumBianchiIX().bp.tolist(),'bm':QuantumBianchiIX().bm.tolist(),'densities':densities})+'\n')
    return result

def regression_signature(d):
    if 'main' not in d:
        vals=[]
        for row in d['selected_comparisons']:
            vals += [row['quantum_to_single_trajectory'],row['quantum_to_classical_ensemble_mean'],row['quantum_rms_width'],row['classical_ensemble_rms_width']]
        for row in d['numerical_check_maxima']: vals += [row['max_mean_distance'],row['max_rms_width_abs_difference']]
        vals += [d['headline']['main_max_edge_mass'],d['headline']['main_max_capped_region_mass']]
        return np.asarray(vals,float)
    vals=[]
    for row in d['comparisons']:
        if row['s'] in (2.0,4.6,6.0): vals += [row['quantum_to_single_trajectory'],row['quantum_to_classical_ensemble_mean'],row['quantum_rms_width'],row['classical_ensemble_rms_width']]
    for row in d['numerical_checks']: vals += [row['max_mean_distance'],row['max_rms_width_abs_difference']]
    vals += [d['headline']['main_max_edge_mass'],d['headline']['main_max_capped_region_mass']]
    return np.asarray(vals,float)

def check_against(result,path):
    old=json.loads(Path(path).read_text());a=regression_signature(result);b=regression_signature(old)
    if a.shape!=b.shape or not np.allclose(a,b,rtol=8e-5,atol=2e-7):
        raise RuntimeError(f'regression mismatch: max_abs={np.max(abs(a-b)) if a.shape==b.shape else "shape"}')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'quantum_bianchi_ix_results');ap.add_argument('--check-against',type=Path);a=ap.parse_args();r=run(a.out)
    if a.check_against: check_against(r,a.check_against)
    print(json.dumps(r['headline'],indent=2))
if __name__=='__main__':main()
