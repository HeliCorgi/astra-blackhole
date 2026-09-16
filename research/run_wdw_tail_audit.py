"""Tail-sensitive observables: fixed model, regulator failures retained.
Copyright 2026 HeliCorgi. Apache-2.0.
"""
from __future__ import annotations
import argparse,hashlib,json,platform,time
from pathlib import Path
import numpy as np
import scipy
from wdw_model import Packet,WDWModel
from wdw_tail_model import states_at,log_moment,threshold_overlap,tail_density,classical_log_moments
ROOT=Path(__file__).resolve().parent
LN10=np.log(10.)


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def packets(protocol):
    return [(p['id'],Packet(**{k:v for k,v in p.items() if k!='id'})) for p in protocol['packets']]


def run_grid(protocol,grid,out,only_h02=False):
    started=time.perf_counter();rows=[];cuts=[];tails=[];errors=[];plots={}
    selected=packets(protocol)
    if only_h02:selected=[p for p in selected if p[0]=='h0.2']
    for h in sorted(set(p.h for _,p in selected)):
        model=WDWModel(h,grid['left'],grid['right'],grid['dx'])
        for pid,packet in selected:
            if packet.h!=h:continue
            states,checks=states_at(model,packet,protocol['times'])
            roundtrip=checks.pop('initial_roundtrip');errors.append({'packet':pid,**checks})
            initial=abs(packet.initial(model.x))**2
            rt=abs(roundtrip)**2
            for j,t in enumerate(protocol['times']):
                mass=abs(states[:,j])**2;x=model.x
                row={'grid':grid['id'],'packet':pid,'T':t,'right':grid['right'],'dx':grid['dx'],
                     'mean_x':float(x@mass),'variance_x':float((x*x)@mass-(x@mass)**2),
                     'rg_over_ref':float(np.exp(packet.x0-x@mass-t)),
                     'last_unit_probability':float(mass[x>grid['right']-1].sum())}
                for q in protocol['powers']:
                    lg=log_moment(x,mass,q,t,packet.x0)
                    row[f'log10_J{q}']=lg/LN10
                    edge=x>grid['right']-1
                    row[f'last_unit_J{q}_fraction']=float(np.exp(log_moment(x[edge],mass[edge],q,t,packet.x0)-lg))
                for threshold in protocol['bounded_tail_thresholds']:
                    row[f'P_x_above_{threshold}']=float(mass[x>=threshold].sum())
                if t==0:
                    for q in protocol['powers']:
                        row[f'exact_initial_log10_J{q}']=q*q*packet.sigma**2/2/LN10
                        row[f'roundtrip_log10_J{q}']=log_moment(x,rt,q,0,packet.x0)/LN10
                rows.append(row)
                if t in (0.5,2,6):
                    for L in protocol['cutoffs']:
                        if L>=grid['right']-1:continue
                        item={'grid':grid['id'],'packet':pid,'T':t,'cutoff':L}
                        for q in protocol['powers']:
                            item[f'log10_partial_J{q}']=log_moment(x,mass,q,t,packet.x0,L)/LN10
                            item[f'log10_capped_J{q}']=log_moment(x,mass,q,t,packet.x0,L,True)/LN10
                        cuts.append(item)
                if t==2:
                    B,_=threshold_overlap(packet)
                    for xx in protocol['tail_points']:
                        if xx>=grid['right']-1:continue
                        idx=np.argmin(abs(x-xx));actual=float(mass[idx]/grid['dx'])
                        asym=float(tail_density(xx,t,h,B))
                        tails.append({'grid':grid['id'],'packet':pid,'x':xx,'density':actual,'asymptotic_density':asym,'ratio':actual/asym})
                    if grid['id'] in ('R128','R256'):
                        keep=(x>=8)&(x<=50)
                        plots[pid]={'x':x[keep].tolist(),'density':(mass[keep]/grid['dx']).tolist(),
                                    'asymptotic_density':tail_density(x[keep],t,h,B).tolist()}
    data={'grid':grid,'rows':rows,'cutoffs':cuts,'tail_points':tails,'checks':errors,'plot_data':plots,
          'seconds':time.perf_counter()-started}
    save(out/'parts'/f"{grid['id']}.json",data)
    print(grid['id'],'states',len(errors),'seconds',round(data['seconds'],2),flush=True)
    return data


def assemble(protocol,out):
    grids=protocol['grids']+protocol['extra_grids_h02_only']
    data=[json.loads((out/'parts'/f"{g['id']}.json").read_text()) for g in grids]
    rows=sum((a['rows'] for a in data),[]);cuts=sum((a['cutoffs'] for a in data),[]);tails=sum((a['tail_points'] for a in data),[])
    checks=sum((a['checks'] for a in data),[])
    initial_error=max(abs(r[f'log10_J{q}']-r[f'exact_initial_log10_J{q}'])*LN10 for r in rows if r['T']==0 for q in protocol['powers'])
    overlaps=[];classical=[]
    for pid,p in packets(protocol):
        B,error=threshold_overlap(p);B2,error2=threshold_overlap(p,14.,1e-14)
        overlaps.append({'packet':pid,'real':B.real,'imag':B.imag,'abs':abs(B),'quad_reported_error':error,
                         'refinement_difference':abs(B-B2),'nonzero_numerically_resolved':bool(abs(B)>100*max(error,error2,abs(B-B2)))})
        for q in protocol['powers']:
            a=classical_log_moments(p,protocol['times'],q,96);b=classical_log_moments(p,protocol['times'],q,144)
            classical.append({'packet':pid,'q':q,'log10_J':(a/LN10).tolist(),'log_refinement_max_difference':float(np.max(abs(a-b))),
                 'quadrature_resolved_to_1e4_in_log':bool(np.max(abs(a-b))<1e-4),
                 'finite_expectation_basis':'Analytic speed bound, not quadrature convergence',
                 'upper_bound_log10':((q*q*p.sigma**2/2+2*q*np.array(protocol['times']))/LN10).tolist()})
    def get(g,p,t):return next(r for r in rows if r['grid']==g and r['packet']==p and r['T']==t)
    def density(g,p,x):return next(r['density'] for r in tails if r['grid']==g and r['packet']==p and r['x']==x)
    tail_grid=max(abs(density('R128','h0.2',x)/density('R128_fine','h0.2',x)-1) for x in protocol['tail_points'])
    slopes=[]
    for g in ('R128','R256'):
        for pid,p in packets(protocol):
            sub=[r for r in tails if r['grid']==g and r['packet']==pid and protocol['tail_fit_interval'][0]<=r['x']<=protocol['tail_fit_interval'][1]]
            if len(sub)<3:continue
            X=np.array([r['x'] for r in sub])+np.log(2*p.h)-np.euler_gamma
            slope=np.polyfit(np.log(X),np.log([r['density'] for r in sub]),1)[0]
            slopes.append({'grid':g,'packet':pid,'density_power_slope':float(slope),'ratios':[r['ratio'] for r in sub]})
    core=[get(g,'h0.2',2) for g in ('R12','R16','R24','R32')]
    fixed_cutoff_comparisons=[]
    for L in (8,12,16,24,32,40,48):
        a=next(r for r in cuts if r['grid']=='R128' and r['packet']=='h0.2' and r['T']==2 and r['cutoff']==L)
        b=next(r for r in cuts if r['grid']=='R256' and r['packet']=='h0.2' and r['T']==2 and r['cutoff']==L)
        fixed_cutoff_comparisons.append({'cutoff':L,'partial_log10_difference':a['log10_partial_J6']-b['log10_partial_J6'],
                                        'capped_log10_difference':a['log10_capped_J6']-b['log10_capped_J6']})
    quality={'max_norm_error':max(r['max_norm_error'] for r in checks),
             'max_direct_state_l2_difference':max(r['max_direct_state_l2_difference'] for r in checks),
             'initial_log_moment_max_error':initial_error,'tail_density_grid_relative_difference_h02':tail_grid,
             'threshold_overlap_max_refinement_difference':max(r['refinement_difference'] for r in overlaps)}
    limits=protocol['checks'];valid=quality['max_norm_error']<limits['norm'] and quality['max_direct_state_l2_difference']<limits['direct_state_l2'] and initial_error<limits['initial_log_moment'] and tail_grid<limits['tail_grid_relative']
    summary={'schema':1,'scope':protocol['scope'],'grid_state_count':len(checks),'time_rows':len(rows),'cutoff_rows':len(cuts),
             'numerical_implementation_checks_pass':bool(valid),'quality':quality,
             'full_unbounded_moment_finite_convergence':'FAILED for q=2,4,6 in this regulator study; no finite value reported',
             'four_dimensional_K_operator':'NOT constructed: quantum mass observable and ordering remain unspecified',
             'continuum_analysis':'Threshold-scattering asymptotic, under stated smooth spectral remainder assumptions: psi~4i B T/(pi X^3). B is numerically nonzero in all 5 inputs. This asymptotic implies infinite positive exponential moments at fixed T>0. Not a Lean theorem or a certified global remainder bound.',
             'core_h02_T2':core,'threshold_overlaps':overlaps,'tail_slopes_diagnostic_only':slopes,
             'fixed_cutoff_domain_comparison':fixed_cutoff_comparisons,
             'classical_ensemble':classical,
             'regulator_source_and_regression':'All detailed rows retained; compare flags/counts and toleranced moderate-range diagnostics, not platform-sensitive enormous noise-amplified roundtrip moments.',
             'source_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['wdw_model.py','wdw_tail_model.py','run_wdw_tail_audit.py','wdw_tail_protocol.json']},
             'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__}}
    save(out/'summary.json',summary);save(out/'observables.json',rows);save(out/'cutoffs.json',cuts);save(out/'tail_points.json',tails)
    save(out/'plot_data.json',{a['grid']['id']:a['plot_data'] for a in data if a['plot_data']})
    save(out/'protocol.json',protocol)
    if not valid:raise RuntimeError('Numerical implementation check failed; raw results retained.')
    return summary


def compare(actual,saved):
    for key in ['schema','grid_state_count','time_rows','cutoff_rows','numerical_implementation_checks_pass','full_unbounded_moment_finite_convergence','four_dimensional_K_operator','source_sha256']:
        if actual[key]!=saved[key]:raise AssertionError(key)
    for a,b in zip(actual['core_h02_T2'],saved['core_h02_T2'],strict=True):
        for key in ['mean_x','log10_J2','log10_J4','log10_J6','last_unit_probability']:
            np.testing.assert_allclose(a[key],b[key],rtol=5e-4,atol=1e-7,err_msg=key)
    for a,b in zip(actual['threshold_overlaps'],saved['threshold_overlaps'],strict=True):
        np.testing.assert_allclose([a['real'],a['imag']],[b['real'],b['imag']],rtol=1e-8,atol=1e-12)
    print('Saved selected-observable regression passed; not a proof of continuum convergence.')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,default=ROOT/'wdw_tail_results')
    ap.add_argument('--phase',choices=['all','core','wide','fine','large','assemble'],default='all');ap.add_argument('--check-against',type=Path)
    args=ap.parse_args();p=json.loads((ROOT/'wdw_tail_protocol.json').read_text());phase=args.phase
    if phase in ('all','core'):
        for g in p['grids'][:-1]:run_grid(p,g,args.out)
    if phase in ('all','wide'):run_grid(p,p['grids'][-1],args.out)
    if phase in ('all','fine'):run_grid(p,p['extra_grids_h02_only'][0],args.out,True)
    if phase in ('all','large'):run_grid(p,p['extra_grids_h02_only'][1],args.out,True)
    if phase in ('all','assemble'):
        result=assemble(p,args.out)
        if args.check_against:compare(result,json.loads(args.check_against.read_text()))
        print(json.dumps(result['quality'],indent=2))
if __name__=='__main__':main()
