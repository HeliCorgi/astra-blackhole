"""Compare three fixed closures without fitting or changing forecast horizons.
Copyright 2026 HeliCorgi. Apache-2.0.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,platform,time
from pathlib import Path
import numpy as np
import scipy
from numpy.polynomial import legendre as leg
from collision_model import RelaxationModel
from closure_model import ClosureModel,initial_state,moment_diagnostics,observe_series
ROOT=Path(__file__).resolve().parent


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def first_time(mask,times):
    ids=np.flatnonzero(mask)
    return None if not len(ids) else float(times[ids[0]])


def execute(protocol: dict,out: Path) -> dict:
    out.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter(); end=protocol['horizon']; th=protocol['thresholds']
    times=np.unique(np.r_[np.linspace(0,end,protocol['linear_time_samples']),
                         np.geomspace(1e-6,end,protocol['early_time_samples'])])
    mu=np.linspace(-1,1,protocol['reference_mu_samples'])
    rows=[];quality=[];examples={};reference_time=0.;reduced_time=0.
    for g in protocol['geometries']:
        for distribution in protocol['distributions']:
            for rate in protocol['rates']:
                cid=f"{g['id']}/{distribution['id']}/nu{rate}"
                config=(g['density'],g['Hb'],distribution['coefficients'])
                reference_obs=[];br=[]
                t0=time.perf_counter()
                for order,method,rtol,atol in [
                    (protocol['reference_lmax'],protocol['reference_method'],protocol['rtol'],protocol['atol']),
                    (protocol['verification_lmax'],protocol['verification_method'],protocol['verification_rtol'],protocol['verification_atol'])]:
                    m=RelaxationModel(order,rate); y0=initial_state(order,*config)
                    sol=m.evolve(end=end,y0=y0,method=method,rtol=rtol,atol=atol)
                    o=observe_series(m,sol,times); reference_obs.append(o)
                    br.append(leg.legvander(mu,order)[:,::2]@o['c'])
                reference_time+=time.perf_counter()-t0
                r,rh=reference_obs
                q=dict(case=cid,K_agreement=float(np.max(abs(r['K']/rh['K']-1))),
                    Pi_agreement=float(np.max(abs(r['Pi']-rh['Pi'])/rh['rho'])),
                    brightness_agreement=float(np.max(abs(br[0]-br[1])/rh['rho'])),
                    minimum_brightness=float(min(np.min(br[0]/r['rho']),np.min(br[1]/rh['rho']))),
                    constraint=float(max(np.max(abs(r['constraint'])),np.max(abs(rh['constraint'])))),
                    number_drift=float(max(np.max(abs(r['number_drift'])),np.max(abs(rh['number_drift'])))))
                q['valid']=(q['K_agreement']<=th['reference_relative_curvature_agreement']
                    and q['Pi_agreement']<=th['reference_pressure_agreement_over_density']
                    and q['brightness_agreement']<=th['reference_brightness_agreement_over_density']
                    and q['minimum_brightness']>=-th['realizability_roundoff_tolerance']
                    and q['constraint']<=th['scaled_constraint']
                    and q['number_drift']<=th['comoving_number_relative_drift'])
                quality.append(q)
                case_curves={'time':times.tolist(),'reference_K':r['K'].tolist(),
                             'reference_Pi_over_rho':(r['Pi']/r['rho']).tolist()}
                for label,order in protocol['closures'].items():
                    t0=time.perf_counter(); m=ClosureModel(order,rate)
                    sol=m.evolve(initial_state(order,*config),end,rtol=protocol['rtol'],atol=protocol['atol'])
                    o=observe_series(m,sol,times); reduced_time+=time.perf_counter()-t0
                    dg=[moment_diagnostics(c) for c in o['c'].T]
                    margins=np.array([d['moment_margin'] for d in dg]); pmin=np.array([d['polynomial_minimum'] for d in dg])
                    ke=abs(o['K']/r['K']-1);pe=abs(o['Pi']-r['Pi'])/r['rho']
                    bad=(ke>th['curvature_relative_error'])|(pe>th['pressure_gap_error_over_reference_density'])
                    aok=not np.any(bad); mok=np.min(margins)>=-th['realizability_roundoff_tolerance']
                    pok=np.min(pmin)>=-th['realizability_roundoff_tolerance']
                    con=float(np.max(abs(o['constraint']))); nd=float(np.max(abs(o['number_drift'])))
                    physics=mok and con<=th['scaled_constraint'] and nd<=th['comoving_number_relative_drift']
                    row=dict(case=cid,split=g['split'],rate=rate,model=label,
                             max_K_error=float(np.max(ke)),max_Pi_error_over_rho=float(np.max(pe)),
                             endpoint_K_error=float(ke[-1]),moment_margin=float(np.min(margins)),
                             polynomial_min=float(np.min(pmin)),constraint=con,number_drift=nd,
                             reference_valid=q['valid'],accuracy_pass=aok,moment_pass=bool(physics),
                             polynomial_pass=bool(pok),accepted=bool(q['valid'] and aok and physics and pok),
                             first_accuracy_failure=first_time(bad,times),
                             first_polynomial_failure=first_time(pmin < -th['realizability_roundoff_tolerance'],times),
                             maximum_epsilon=None if o['epsilon'] is None else float(max(o['epsilon'])))
                    rows.append(row)
                    case_curves[label]={'K':o['K'].tolist(),'Pi_over_reference_rho':(o['Pi']/r['rho']).tolist(),
                                        'polynomial_minimum':pmin.tolist(),'K_error':ke.tolist(),'Pi_error':pe.tolist()}
                if g['id']=='development' and distribution['id']=='P4_plus': examples[str(rate)]=case_curves
                print(cid, 'reference',q['valid'], 'accepted',','.join(r['model'] for r in rows[-3:] if r['accepted']) or 'none',flush=True)
    aggregate=[]
    for split in ['development','heldout']:
        for rate in protocol['rates']:
            for model in protocol['closures']:
                sub=[r for r in rows if r['split']==split and r['rate']==rate and r['model']==model]
                aggregate.append(dict(split=split,rate=rate,model=model,cases=len(sub),
                    accepted=sum(r['accepted'] for r in sub),accuracy_pass=sum(r['accuracy_pass'] for r in sub),
                    moment_pass=sum(r['moment_pass'] for r in sub),polynomial_pass=sum(r['polynomial_pass'] for r in sub),
                    worst_K_error=max(r['max_K_error'] for r in sub),worst_Pi_error=max(r['max_Pi_error_over_rho'] for r in sub)))
    grouped={q['case']:[r for r in rows if r['case']==q['case']] for q in quality}
    minima={}
    for split in ['development','heldout']:
        minima[split]={k:0 for k in ['A','B','C','none']}
        for case in grouped.values():
            if case[0]['split']==split:
                passing=[r['model'] for r in case if r['accepted']]
                minima[split][passing[0] if passing else 'none']+=1
    summary=dict(schema=1,scope='Fixed-interval synthetic classical closure benchmark; not real-world validation',
       protocol_sha256=sha(ROOT/'closure_protocol.json'),
       sources={str(p.relative_to(ROOT.parent)):sha(p) for p in [ROOT/'collision_model.py',ROOT/'closure_model.py',Path(__file__)]},
       environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
       cases=len(quality),reference_trajectories=2*len(quality),reduced_trajectories=len(rows),
       time_samples=len(times),horizon=end,reference_all_valid=all(q['valid'] for q in quality),
       reference_quality={key:(min(q[key] for q in quality) if key=='minimum_brightness' else max(q[key] for q in quality))
           for key in ['K_agreement','Pi_agreement','brightness_agreement','minimum_brightness','constraint','number_drift']},
       aggregates=aggregate,smallest_accepted_tested_closure=minima,
       negative_polynomial_but_feasible_moments=sum(r['moment_pass'] and not r['polynomial_pass'] for r in rows),
       higher_order_not_better_in_K=sum(c[2]['max_K_error']>c[1]['max_K_error'] for c in grouped.values()),
       elapsed_seconds=time.perf_counter()-started,reference_seconds=reference_time,reduced_seconds=reduced_time)
    buf=io.StringIO(newline=''); writer=csv.DictWriter(buf,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    text=buf.getvalue().encode(); (out/'cases.csv').write_bytes(text)
    (out/'cases.csv.gz').write_bytes(gzip.compress(text,mtime=0))
    summary['cases_csv_sha256']=hashlib.sha256(text).hexdigest()
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (out/'reference_quality.json').write_text(json.dumps(quality,indent=2,allow_nan=False)+'\n')
    (out/'examples.json').write_text(json.dumps(examples,separators=(',',':'),allow_nan=False))
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    if not summary['reference_all_valid']:
        raise RuntimeError('Reference validation failed. Raw outcomes saved, not accepted.')
    return summary


def compare_saved(out: Path, baseline: Path):
    """Regression checks; never replace computation with stored answers."""
    a=json.loads((out/'summary.json').read_text()); b=json.loads((baseline/'summary.json').read_text())
    for key in ['protocol_sha256','cases','reference_all_valid','smallest_accepted_tested_closure',
                'negative_polynomial_but_feasible_moments','higher_order_not_better_in_K']:
        if a[key]!=b[key]: raise AssertionError(f'Changed regression field: {key}')
    aa,bb=a['aggregates'],b['aggregates']
    if len(aa)!=len(bb): raise AssertionError('Changed aggregate schema.')
    for x,y in zip(aa,bb):
        if x.keys()!=y.keys(): raise AssertionError('Changed aggregate fields.')
        for key in x:
            if x[key]==y[key]: continue
            if key not in ('worst_K_error','worst_Pi_error') or not np.isclose(x[key],y[key],rtol=2e-5,atol=2e-9):
                raise AssertionError(f'Aggregate {x["split"]}/{x["rate"]}/{x["model"]}: {key} changed')
    print('Saved aggregate regression passed; all per-case outcomes are retained as generated CSV, not compared individually.')



def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'closure_results')
    parser.add_argument('--check-against',type=Path)
    args=parser.parse_args()
    result=execute(json.loads((ROOT/'closure_protocol.json').read_text()),args.out)
    if args.check_against: compare_saved(args.out,args.check_against)
    print(json.dumps({k:result[k] for k in ['cases','reference_all_valid','smallest_accepted_tested_closure']},indent=2))
if __name__=='__main__':main()
