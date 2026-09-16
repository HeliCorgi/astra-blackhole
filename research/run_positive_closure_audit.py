"""Compare positive angular closures at unchanged retained variable count.
Copyright 2026 HeliCorgi. Apache-2.0.

All previous cases are reused, not a fresh held-out population. Every reduced
model receives only its retained initial moments. Full reference tails enter
an explicitly labelled AFTER-THE-FACT diagnostic, never its evolution.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,platform,time
from pathlib import Path
import numpy as np
import scipy
from numpy.polynomial import legendre as leg
from collision_model import RelaxationModel
from closure_model import ClosureModel,initial_state,moment_diagnostics,observe_series
from positive_closure_model import PositiveClosureModel,AngularMaxEntropy
ROOT=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def write_csv(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def summary_of(rows,quality,protocol):
    methods=protocol['methods'];aggregates=[]
    for split in ('development','heldout'):
        for rate in protocol['rates']:
            for method in methods:
                sub=[r for r in rows if r['split']==split and r['rate']==rate and r['method']==method]
                good=[r for r in sub if r['integration_success']]
                aggregates.append(dict(split=split,rate=rate,method=method,cases=len(sub),
                    accepted=sum(r['accepted'] for r in sub),accuracy_pass=sum(r['accuracy_pass'] for r in sub),
                    distribution_pass=sum(r['distribution_pass'] for r in sub),
                    failures=sum(not r['integration_success'] for r in sub),
                    worst_K_error=max((r['max_K_error'] for r in good),default=None),
                    worst_Pi_error=max((r['max_Pi_error_over_rho'] for r in good),default=None)))
    paired=[]
    for label in ('B','C'):
        for split in ('development','heldout','all'):
            pairs=[]
            for q in quality:
                z=next(r for r in rows if r['case']==q['case'] and r['method']==label+'_zero')
                p=next(r for r in rows if r['case']==q['case'] and r['method']==label+'_positive')
                if split!='all' and z['split']!=split:continue
                pairs.append((z,p))
            paired.append(dict(label=label,split=split,cases=len(pairs),
                zero_accepted=sum(a['accepted'] for a,b in pairs),positive_accepted=sum(b['accepted'] for a,b in pairs),
                zero_accuracy_pass=sum(a['accuracy_pass'] for a,b in pairs),positive_accuracy_pass=sum(b['accuracy_pass'] for a,b in pairs),
                zero_distribution_fail=sum(not a['distribution_pass'] for a,b in pairs),
                positive_distribution_fail=sum(not b['distribution_pass'] for a,b in pairs),
                rescued=sum(not a['accepted'] and b['accepted'] for a,b in pairs),
                lost=sum(a['accepted'] and not b['accepted'] for a,b in pairs),
                both_fail=sum(not a['accepted'] and not b['accepted'] for a,b in pairs),
                K_improves=sum(b['integration_success'] and b['max_K_error']<a['max_K_error'] for a,b in pairs),
                K_worsens=sum(b['integration_success'] and b['max_K_error']>a['max_K_error'] for a,b in pairs),
                Pi_improves=sum(b['integration_success'] and b['max_Pi_error_over_rho']<a['max_Pi_error_over_rho'] for a,b in pairs),
                Pi_worsens=sum(b['integration_success'] and b['max_Pi_error_over_rho']>a['max_Pi_error_over_rho'] for a,b in pairs)))
    minima={}
    for family in ('zero','positive'):
        minima[family]={}
        for split in ('development','heldout'):
            counts={x:0 for x in ['A','B','C','none']}
            for q in quality:
                rr={r['method']:r for r in rows if r['case']==q['case']}
                if rr['A_zero']['split']!=split:continue
                names=['A_zero','B_'+family,'C_'+family]
                passing=[n[0] for n in names if rr[n]['accepted']]
                counts[passing[0] if passing else 'none']+=1
            minima[family][split]=counts
    return dict(aggregates=aggregates,paired_comparisons=paired,smallest_accepted=minima)


def execute(protocol,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();end=protocol['horizon'];th=protocol['thresholds']
    times=np.unique(np.r_[np.linspace(0,end,protocol['linear_time_samples']),
                         np.geomspace(1e-6,end,protocol['early_time_samples'])])
    mu=np.linspace(-1,1,protocol['reference_mu_samples'])
    rows=[];quality=[];curves={};diagnostics=[]
    for geometry in protocol['geometries']:
        for distribution in protocol['distributions']:
            for rate in protocol['rates']:
                cid=f"{geometry['id']}/{distribution['id']}/nu{rate}"
                config=(geometry['density'],geometry['Hb'],distribution['coefficients'])
                refs=[];br=[]
                for L,method,rtol,atol in [
                    (protocol['reference_lmax'],protocol['reference_method'],protocol['rtol'],protocol['atol']),
                    (protocol['verification_lmax'],protocol['verification_method'],protocol['verification_rtol'],protocol['verification_atol'])]:
                    m=RelaxationModel(L,rate)
                    s=m.evolve(end=end,y0=initial_state(L,*config),method=method,rtol=rtol,atol=atol)
                    o=observe_series(m,s,times);refs.append(o)
                    br.append(leg.legvander(mu,L)[:,::2]@o['c'])
                r,rh=refs
                q=dict(case=cid,K_agreement=float(max(abs(r['K']/rh['K']-1))),
                    Pi_agreement=float(max(abs(r['Pi']-rh['Pi'])/rh['rho'])),
                    brightness_agreement=float(np.max(abs(br[0]-br[1])/rh['rho'])),
                    minimum_brightness=float(min(np.min(br[0]/r['rho']),np.min(br[1]/rh['rho']))),
                    constraint=float(max(max(abs(r['constraint'])),max(abs(rh['constraint'])))),
                    number_drift=float(max(max(abs(r['number_drift'])),max(abs(rh['number_drift'])))))
                q['valid']=(q['K_agreement']<=th['reference_relative_curvature_agreement'] and
                    q['Pi_agreement']<=th['reference_pressure_agreement_over_density'] and
                    q['brightness_agreement']<=th['reference_brightness_agreement_over_density'] and
                    q['minimum_brightness']>=-th['realizability_roundoff_tolerance'] and
                    q['constraint']<=th['scaled_constraint'] and q['number_drift']<=th['comoving_number_relative_drift'])
                quality.append(q)
                cur={'time':times.tolist(),'reference_K':r['K'].tolist(),'reference_rho':r['rho'].tolist(),
                     'reference_Pi':r['Pi'].tolist()}
                for method,L in protocol['methods'].items():
                    row=dict(case=cid,split=geometry['split'],rate=rate,method=method,retained_order=L,
                        integration_success=False,failure=None,max_K_error=None,max_Pi_error_over_rho=None,
                        endpoint_K_error=None,moment_margin=None,distribution_minimum=None,
                        constraint=None,number_drift=None,moment_residual=None,tail_quadrature_residual=None,
                        evolution_K_agreement=None,evolution_Pi_agreement=None,max_hessian_condition=None,
                        verification_pass=False,reference_valid=q['valid'],accuracy_pass=False,
                        moment_pass=False,distribution_pass=False,accepted=False)
                    try:
                        positive=method.endswith('_positive')
                        m=(PositiveClosureModel(L,rate,protocol['positive_quadrature'],protocol['positive_tolerance'])
                           if positive else ClosureModel(L,rate))
                        y0=initial_state(L,*config)
                        s=m.evolve(y0,end,rtol=protocol['rtol'],atol=protocol['atol'])
                        o=observe_series(m,s,times)
                        dg=[moment_diagnostics(c) for c in o['c'].T]
                        margins=np.array([d['moment_margin'] for d in dg])
                        minimum=np.array([d['polynomial_minimum'] for d in dg])
                        if positive:
                            e=m.reconstruction;fits=[e.fit(c) for c in o['c'].T]
                            minlogs=np.array([e.minimum_log_density(f) for f in fits]);minimum=np.exp(minlogs)
                            checks=[e.independent_residual(c,f,protocol['positive_verification_quadrature']) for c,f in zip(o['c'].T,fits)]
                            mr=max(x['moment_residual'] for x in checks);tr=max(x['tail_residual'] for x in checks)
                            v=PositiveClosureModel(L,rate,protocol['positive_verification_quadrature'],protocol['positive_verification_tolerance'])
                            vs=v.evolve(y0,end,method='DOP853',rtol=protocol['verification_rtol'],atol=protocol['verification_atol'])
                            vo=observe_series(v,vs,times)
                            ka=float(max(abs(o['K']/vo['K']-1)));pa=float(max(abs(o['Pi']-vo['Pi'])/vo['rho']))
                            verified=mr<=th['positive_moment_residual'] and tr<=th['positive_quadrature_tail_agreement'] and max(ka,pa)<=th['positive_evolution_agreement']
                            row.update(moment_residual=mr,tail_quadrature_residual=tr,evolution_K_agreement=ka,
                                evolution_Pi_agreement=pa,max_hessian_condition=max(f['hessian_condition'] for f in fits))
                            cur[method+'_lambda']=[f['multipliers'].tolist() for f in fits]
                            # Reference-tail reconstruction is only an explanatory diagnostic.
                            errors=[];zeros=[]
                            for it in range(0,len(times),10):
                                c=r['c'][:L//2+1,it];f=e.fit(c);truth=r['c'][L//2+1,it]
                                errors.append(abs(f['tail']-truth)/c[0]);zeros.append(abs(truth)/c[0])
                            diagnostics.append(dict(case=cid,method=method,
                                kind='reference-state reconstruction ONLY; no feedback to evolution',
                                max_omitted_error_positive=float(max(errors)),max_omitted_error_zero=float(max(zeros))))
                        else:
                            verified=True
                        ke=abs(o['K']/r['K']-1);pe=abs(o['Pi']-r['Pi'])/r['rho']
                        con=float(max(abs(o['constraint'])));nd=float(max(abs(o['number_drift'])))
                        aok=bool(max(ke)<=th['curvature_relative_error'] and max(pe)<=th['pressure_gap_error_over_reference_density'])
                        mok=bool(min(margins)>=-th['realizability_roundoff_tolerance'] and con<=th['scaled_constraint'] and nd<=th['comoving_number_relative_drift'])
                        pok=bool(min(minimum)>0 if positive else min(minimum)>=-th['realizability_roundoff_tolerance'])
                        row.update(integration_success=True,max_K_error=float(max(ke)),max_Pi_error_over_rho=float(max(pe)),
                            endpoint_K_error=float(ke[-1]),moment_margin=float(min(margins)),distribution_minimum=float(min(minimum)),
                            constraint=con,number_drift=nd,verification_pass=bool(verified),accuracy_pass=aok,moment_pass=mok,
                            distribution_pass=pok,accepted=bool(q['valid'] and verified and aok and mok and pok))
                        cur[method]=dict(c=o['c'].tolist(),K=o['K'].tolist(),Pi=o['Pi'].tolist(),rho=o['rho'].tolist(),
                            distribution_minimum=minimum.tolist(),K_error=ke.tolist(),Pi_error=pe.tolist())
                    except (ValueError,RuntimeError,np.linalg.LinAlgError,FloatingPointError) as exc:
                        row['failure']=str(exc)
                    rows.append(row)
                curves[cid]=cur
                with (out/'progress.jsonl').open('a') as f:f.write(json.dumps({'reference':q,'rows':rows[-5:]},allow_nan=False)+'\n')
                print(cid,'valid reference:',q['valid'],'passes:',[x['method'] for x in rows[-5:] if x['accepted']],flush=True)
    summary=dict(schema=1,scope='Same-state positive angular-closure comparison; model-internal, not new physical observations',
        base_commit=protocol['base_commit'],protocol_sha256=sha(ROOT/'positive_closure_protocol.json'),
        source_sha256={str(p.relative_to(ROOT.parent)):sha(p) for p in [ROOT/'collision_model.py',ROOT/'closure_model.py',ROOT/'positive_closure_model.py',Path(__file__)]},
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        horizon=end,time_samples=len(times),cases=len(quality),reference_trajectories=2*len(quality),
        reduced_trajectories=len(rows),positive_verification_trajectories=sum(x['method'].endswith('positive') and x['integration_success'] for x in rows),
        all_references_valid=all(q['valid'] for q in quality),all_reduced_succeeded=all(x['integration_success'] for x in rows),
        all_positive_verifications_passed=all(x['verification_pass'] for x in rows if x['method'].endswith('positive')),
        elapsed_seconds=time.perf_counter()-start,**summary_of(rows,quality,protocol))
    write_csv(out/'cases.csv',rows);dump(out/'cases.json',rows);dump(out/'reference_quality.json',quality)
    dump(out/'tail_diagnostics.json',diagnostics);dump(out/'summary.json',summary);dump(out/'protocol.json',protocol)
    # Full curves are artifacts, never inputs to the predictive closures.
    (out/'curves.json').write_text(json.dumps(curves,separators=(',',':'),allow_nan=False)+'\n')
    if not summary['all_references_valid'] or not summary['all_reduced_succeeded'] or not summary['all_positive_verifications_passed']:
        raise RuntimeError('Numerical/domain failures saved; do not present them as passed.')
    return summary


def check_saved(out,baseline):
    a=json.loads((out/'cases.json').read_text());b=json.loads((baseline/'cases.json').read_text())
    if len(a)!=len(b):raise AssertionError('Different number of cases.')
    count=0
    for x,y in zip(a,b):
        if x.keys()!=y.keys():raise AssertionError('Changed fields.')
        for k in x:
            if isinstance(x[k],float):
                if not np.isclose(x[k],y[k],rtol=5e-5,atol=2e-8):raise AssertionError((x['case'],x['method'],k,x[k],y[k]))
            elif x[k]!=y[k]:raise AssertionError((x['case'],k))
            count+=1
    print(f'Per-row regression passed: {len(a)} rows, {count} fields.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=ROOT/'positive_closure_results')
    p.add_argument('--check-against',type=Path);args=p.parse_args()
    r=execute(json.loads((ROOT/'positive_closure_protocol.json').read_text()),args.out)
    if args.check_against:check_saved(args.out,args.check_against)
    print(json.dumps({k:r[k] for k in ['cases','all_references_valid','all_reduced_succeeded','smallest_accepted']},indent=2))
if __name__=='__main__':main()
