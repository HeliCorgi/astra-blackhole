"""Reproduce one published WDW sector and audit same-solution coordinates."""
from __future__ import annotations
import argparse,csv,gc,hashlib,json,platform,time
from pathlib import Path
from dataclasses import asdict
import numpy as np
import scipy
from chiba_model import Data,Diamond,green,exact_slice,same_slice_flux
ROOT=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def execute(out: Path):
    protocol=json.loads((ROOT/'chiba_protocol.json').read_text());out.mkdir(parents=True,exist_ok=True)
    t0=time.perf_counter();rows=[];quality=[];curves={};fielddata={}
    times=np.array(protocol['clock_times']);ct=np.linspace(0,times[-1],protocol['curve_times'])
    up,vp=np.meshgrid(protocol['field_probe_U'],protocol['field_probe_V']);up=up.ravel();vp=vp.ravel()
    for kappa in protocol['kappas']:
        data=Data(kappa,protocol['sigma'],protocol['rh'],protocol['ub']);allrows=[];allfields=[]
        exact=green(data,up,vp,protocol['green_order'])
        tighter=green(data,up,vp,protocol['green_verification_order'])
        gr=max(float(np.max(abs(a-b))) for a,b in zip(exact,tighter))
        for step in protocol['steps']:
            model=Diamond(data,step,protocol['umax'],protocol['vmax'])
            r=[model.slice(t) for t in times];allrows.append(r)
            for x in r:rows.append(dict(kappa=kappa,step=step,**x))
            allfields.append(model.evaluate(up,vp)[0])
            if step==protocol['steps'][-1]:
                curve=[model.slice(t) for t in ct];curves[str(kappa)]=curve
                fu=np.linspace(-.5,.3,161);fv=np.linspace(0,.7,141)
                # The field itself, including the boundary, needs no derivative margin.
                V,U=np.meshgrid(fv,fu,indexing='ij')
                env=model.interp(np.c_[V.ravel(),U.ravel()]).reshape(U.shape)
                fielddata[str(kappa)]=dict(U=fu.tolist(),V=fv.tolist(),modulus_squared=(abs(env)**2).tolist())
            del model;gc.collect()
        refs=[exact_slice(data,t,order=protocol['green_verification_order'],samples=protocol['green_slice_samples']) for t in times]
        fine=allrows[-1];med=allrows[-2]
        Q=2*np.sqrt(np.pi)*data.sigma/data.kappa
        q=dict(kappa=kappa,incoming_charge=Q,
               max_fine_medium_X=max(abs(a['mean_X']-b['mean_X']) for a,b in zip(fine,med)),
               max_coarse_medium_X=max(abs(a['mean_X']-b['mean_X']) for a,b in zip(allrows[0],med)),
               max_field_difference_from_green=float(max(abs(allfields[-1]-exact[0]))),
               max_X_difference_from_green=max(abs(a['mean_X']-b['mean_X']) for a,b in zip(fine,refs)),
               max_norm_relative=max(abs(a['norm']/Q-1) for a in fine),
               green_refinement_max_absolute=gr,
               max_negative_over_absolute=max(a['negative_over_absolute'] for a in fine),
               green_clock_slices=refs)
        # Crossing of a signed expectation, NOT a singularity-hitting probability.
        c=curves[str(kappa)];us=np.array([r['mean_U'] for r in c]);idx=np.flatnonzero(us>=0)
        q['mean_U_zero_time']=None if not len(idx) else float(np.interp(0,us[idx[0]-1:idx[0]+1],ct[idx[0]-1:idx[0]+1]))
        quality.append(q);print('kappa',kappa,'X error',q['max_X_difference_from_green'],'crossing',q['mean_U_zero_time'],flush=True)
    mapchecks=[same_slice_flux(Data(kappa=protocol['pullback_kappa']),t,
                              protocol['pullback_upper_U'],protocol['pullback_lower_U']) for t in protocol['pullback_times']]
    base=Data(kappa=protocol['boundary_check_kappa']);far=Data(kappa=base.kappa,ub=protocol['boundary_check_ub'])
    boundary=[]
    for t in (0.,.1,.3,.6):
        a=exact_slice(base,t,order=384,samples=2001);b=exact_slice(far,t,order=384,samples=2001)
        boundary.append(dict(t=t,X_difference=abs(a['mean_X']-b['mean_X']),norm_difference=abs(a['norm']-b['norm'])))
    plateau=[]
    for kappa in protocol['kappas']:
        d=Data(kappa=kappa);g0=d.boundary(0.)
        for x in protocol['plateau_x']:
            u=-np.exp(-x)/4;v=-u;p,du,dv=green(d,u,v,384)
            # t_M=0 <=> boost tau=0, and q_tau is a conserved-current density.
            dtau=-u*du+v*dv
            plateau.append(dict(kappa=kappa,x=x,modulus_squared=float(abs(p)**2),
                       plateau_squared=float(abs(g0)**2),ratio_to_plateau=float(abs(p/g0)**2),
                       signed_KG_density_tau=float(-2*np.imag(np.conj(p)*dtau))))
    th=protocol['thresholds'];checks={
        'field_green':all(q['max_field_difference_from_green']<th['finite_difference_green_field_absolute'] for q in quality),
        'expectation_green':all(q['max_X_difference_from_green']<th['finite_difference_green_mean_X_absolute'] for q in quality),
        'grid_refinement':all(q['max_fine_medium_X']<th['fine_medium_mean_X_absolute'] for q in quality),
        'charge':all(q['max_norm_relative']<th['fine_norm_relative_to_incoming'] for q in quality),
        'green_refinement':all(q['green_refinement_max_absolute']<th['green_refinement_absolute'] for q in quality),
        'coordinate_flux':all(q['relative_difference']<th['same_slice_flux_relative'] for q in mapchecks),
        'boundary':all(q['X_difference']<th['distant_boundary_mean_X_absolute'] for q in boundary)}
    result=dict(schema=1,protocol=protocol,sources={n:sha(ROOT/n) for n in ('chiba_model.py','chiba_protocol.json',Path(__file__).name)},
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        quality=quality,same_slice_coordinate_checks=mapchecks,boundary_checks=boundary,
        plateau=dict(amplitude=float(abs(base.boundary(0.))),squared=float(abs(base.boundary(0.))**2),
                     conclusion='Same scalar pullback has nonzero constant tail at finite boost time; not an L2(dx) state. No projection/reinitialization performed.'),
        checks=checks,all_checks_pass=all(checks.values()),elapsed_seconds=time.perf_counter()-t0,
        limitations='Published-equation reproduction, not raw-author-data match, peer review, positive probability measure or singularity-resolution proof.')
    for filename,obj in [('summary.json',result),('curves.json',curves),('fields.json',fielddata),('plateau.json',plateau)]:save(out/filename,obj)
    with (out/'clock_slices.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    if not result['all_checks_pass']:raise RuntimeError('Checks failed; full results preserved, not accepted.')
    return result


def compare(a,b):
    if a['protocol']!=b['protocol'] or a['sources']!=b['sources'] or a['checks']!=b['checks']:raise AssertionError('Protocol/source/check change')
    for x,y in zip(a['quality'],b['quality'],strict=True):
        if x['kappa']!=y['kappa']:raise AssertionError('Case change')
        for xx,yy in zip(x['green_clock_slices'],y['green_clock_slices'],strict=True):
            for field in ('norm','mean_X','mean_U','mean_V'):
                np.testing.assert_allclose(xx[field],yy[field],rtol=1e-7,atol=2e-9)
        for field in ('mean_U_zero_time',):
            if x[field] is None or y[field] is None:
                if x[field]!=y[field]:raise AssertionError('Crossing classification changed')
            else:np.testing.assert_allclose(x[field],y[field],rtol=1e-6,atol=2e-8)
    print('Selected per-clock observable regression passed; not equality of every field.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=ROOT/'chiba_results');p.add_argument('--check-against',type=Path);a=p.parse_args()
    r=execute(a.out)
    if a.check_against:compare(r,json.loads(a.check_against.read_text()))
    print(json.dumps(r['checks'],indent=2))
if __name__=='__main__':main()
