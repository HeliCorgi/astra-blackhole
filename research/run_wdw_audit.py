"""Reproduce the finite-clock quantum geometry audit. Apache-2.0; HeliCorgi 2026."""
from __future__ import annotations
from pathlib import Path
import argparse,hashlib,json,platform,time
import numpy as np
import scipy
from wdw_model import Packet,WDWModel,classical_path,classical_ensemble
ROOT=Path(__file__).resolve().parent

def serial(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x

def execute(out):
    start=time.perf_counter();p=json.loads((ROOT/'wdw_protocol.json').read_text());out.mkdir(parents=True,exist_ok=True)
    ts=np.linspace(0,p['time_end'],p['time_samples']);rows=[];curves={};checks=[]
    packets=[(f'h{h}',Packet(h,np.sqrt(h),x0=p['x0'],p0=p['p0'])) for h in p['h_values']]
    pair=p['phase_pair']
    packets += [(f'phase{g:+}',Packet(pair['h'],pair['sigma'],g,p['x0'],p['p0'])) for g in pair['chirps']]
    px,pp=classical_path(p['x0'],p['p0'],ts)
    for name,packet in packets:
        m=WDWModel(packet.h,**p['grid']);o=m.evolve(packet,ts,True)
        f=WDWModel(packet.h,**{**p['grid'],'dx':p['refinement_dx']}).evolve(packet,ts)
        b=WDWModel(packet.h,**p['expanded_box']).evolve(packet,ts)
        e=classical_ensemble(packet,ts,p['ensemble_orders'][0]);ef=classical_ensemble(packet,ts,p['ensemble_orders'][1])
        ck=dict(case=name,norm=o['norm_error'],energy=o['energy_relative_drift'],discrete_constraint=o['discrete_constraint_relative_residual'],
                grid_mean_x=float(max(abs(o['mean_x']-f['mean_x']))),box_mean_x=float(max(abs(o['mean_x']-b['mean_x']))),
                box_variance=float(max(abs(o['var_x']-b['var_x']))),ensemble_mean_x=float(max(abs(e['mean_x']-ef['mean_x']))))
        ck['passed']=all(ck[k]<v for k,v in p['thresholds'].items());checks.append(ck)
        row=dict(case=name,h=packet.h,sigma=packet.sigma,chirp=packet.chirp,initial=o['initial'],analytic_initial=packet.continuum_moments(),
                 end_mean_x=float(o['mean_x'][-1]),end_var_x=float(o['var_x'][-1]),
                 end_radius_geometric_mean_ratio=float(o['geometric_mean_radius_ratio'][-1]),
                 max_log_radius_difference_from_point=float(max(abs(o['mean_x']-px))),
                 max_log_radius_difference_from_ensemble=float(max(abs(o['mean_x']-e['mean_x']))),
                 sampled_radius_geometric_mean_monotone_decreasing=bool(np.all(np.diff(o['mean_log_radius_ratio'])<0)),
                 x_turn_sample_time=float(ts[np.argmin(o['mean_x'])]),max_edge_mass=o['max_edge_mass'])
        rows.append(row)
        curves[name]={k:o[k] for k in ['mean_x','var_x','mean_log_radius_ratio','geometric_mean_radius_ratio']}
        curves[name]['ensemble_mean_x']=e['mean_x'];curves[name]['ensemble_var_x']=e['var_x']
        # Wave density only for one example, generated artifact rather than source control.
        if name=='h0.2':np.savez_compressed(out/'wave_example.npz',x=o['x'],time=ts,density=o['density'])
        print(name,json.dumps(row,allow_nan=False),json.dumps(ck),flush=True)
        del m,o,f,b
    # Deliberately too-short box: test that an apparently stable unitary solution can still be wrong.
    packet=packets[-1][1];bad=WDWModel(packet.h,**p['negative_control']).evolve(packet,ts)
    true=curves[packets[-1][0]]
    negative=dict(max_mean_x_disagreement=float(max(abs(bad['mean_x']-true['mean_x']))),norm_error=bad['norm_error'],max_edge_mass=bad['max_edge_mass'])
    negative['detected']=negative['max_mean_x_disagreement']>p['thresholds']['box_mean_x']
    a,b=rows[-2:];pair_result=dict(initial_matched=['x marginal density','continuum p marginal density','mean_x','mean_p','var_x','var_p'],
       not_matched=['x-p covariance','mean_H','all metric time derivatives'],
       end_log_radius_difference=abs(a['end_mean_x']-b['end_mean_x']),
       end_radius_geometric_mean_symmetric_gap=2*abs(a['end_radius_geometric_mean_ratio']-b['end_radius_geometric_mean_ratio'])/(a['end_radius_geometric_mean_ratio']+b['end_radius_geometric_mean_ratio']),
       end_variances=[a['end_var_x'],b['end_var_x']],
       classical_ensemble_end_mean_x_gap=abs(curves[packets[-2][0]]['ensemble_mean_x'][-1]-curves[packets[-1][0]]['ensemble_mean_x'][-1]))
    summary=dict(schema=1,protocol=p,cases=rows,comparisons=checks,phase_pair=pair_result,negative_control=negative,
      all_checks_passed=all(c['passed'] for c in checks) and negative['detected'],
      classical=dict(end_x=float(px[-1]),end_radius_ratio=float(np.exp(p['x0']-px[-1]-ts[-1])),
                     x_turn_time=float(-np.arcsinh(p['p0']*np.exp(p['x0']))),radius_log_derivative_at_x_turn=-1.),
      sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['wdw_model.py','run_wdw_audit.py','wdw_protocol.json']},
      environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),elapsed_seconds=time.perf_counter()-start,
      conclusion='No bounce of the sampled geometric mean areal radius. No singularity resolution or quantum curvature finiteness demonstrated. Reduced quantum dynamics, not a claim about the full theory.')
    (out/'summary.json').write_text(json.dumps(serial(summary),indent=2,allow_nan=False)+'\n')
    (out/'curves.json').write_text(json.dumps(serial(dict(time=ts,point_mean_x=px,cases=curves)),allow_nan=False)+'\n')
    if not summary['all_checks_passed']:raise RuntimeError('Audit failed; results preserved without relaxing acceptance criteria.')
    return summary

def compare(actual,baseline):
    b=json.loads(baseline.read_text())
    assert actual['protocol']==b['protocol'] and actual['sources']==b['sources']
    assert actual['all_checks_passed'] and len(actual['cases'])==len(b['cases'])
    for x,y in zip(actual['cases'],b['cases']):
        assert (x['case'],x['sampled_radius_geometric_mean_monotone_decreasing'])==(y['case'],y['sampled_radius_geometric_mean_monotone_decreasing'])
        for k in ['end_mean_x','end_var_x','end_radius_geometric_mean_ratio','max_log_radius_difference_from_point','max_log_radius_difference_from_ensemble']:
            np.testing.assert_allclose(x[k],y[k],rtol=2e-5,atol=2e-8)
    print('Saved physical diagnostics and protocol/source fingerprints agree; not a physical validity guarantee.')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,default=ROOT/'wdw_results');ap.add_argument('--check-against',type=Path)
    args=ap.parse_args();result=execute(args.out)
    if args.check_against:compare(result,args.check_against)
