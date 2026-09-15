"""Reproduce a hierarchy of matched-state counterexamples.

No model is trained here. All results are synthetic Einstein--Vlasov solutions.
The first omitted moment is varied while lower moments and geometry are held
fixed. Paired futures score distinguishability, not an AI accuracy contest.
"""
from __future__ import annotations
import argparse,json,platform,sys,time
from pathlib import Path
import numpy as np
import scipy
from kinetic_model import MomentumBasis,matched_pair,massless_pair
from local_jets import CurvatureJets
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
from predictors import polynomial

def dump(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def sample(basis,solution,weights,times):
    states=solution.sol(times).T
    obs=[basis.observe(y,weights) for y in states]
    return {k:np.array([o[k] for o in obs]) for k in obs[0]}

def pair_run(basis,pair,times):
    observations=[];solutions=[]
    for key in ('plus','minus'):
        sol=basis.evolve(pair[key],0.,times[-1])
        solutions.append(sol)
        observations.append(sample(basis,sol,pair[key],times))
    op,om=observations
    gap=200*abs(op['K']-om['K'])/(op['K']+om['K'])
    return observations,solutions,gap

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/'results')
    ap.add_argument('--cohort-per-level',type=int,default=8)
    args=ap.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True)
    t0=time.perf_counter()
    centers=np.geomspace(.05,8.,10)
    basis=MomentumBasis(centers)
    times=np.linspace(0.,.45,451)
    protocol=dict(date='2026-09-16',status='Exploratory numerical counterexample audit, not preregistered.',
      model='Kantowski--Sachs Einstein--Vlasov, c=8*pi*G=1, rest mass=1, Lambda=0.',
      geometry='ds^2=-dt^2+a^2 dx^2+b^2 dOmega^2; a0=b0=1; Hb0=-1; rho0=.8; Ha0=.6.',
      distribution='10 smooth lognormal radial-number kernels; initially isotropic directions; no particle correlations.',
      centers=centers,log_width=.10,radial_quadrature=12,angular_quadrature=32,
      initial_same='particle density n, energy density rho, both directional pressures, metric, extrinsic curvature, K.',
      levels=[1,2,3],matching='match N and M0..M_k; bound |delta weight| <= .8 * baseline weight.',
      selection='LP maximizes first omitted INITIAL moment only, never future K.',
      integration='DOP853 rtol=2e-11 atol=2e-13; finite-radius safety stop b=.1.',
      evaluation_times=times,cohort_seeds=list(range(601,601+args.cohort_per_level)),
      curve_gap_definition='200*abs(KA-KB)/(KA+KB) percent; this is pair separation, NOT prediction error.',
      minimax_lower_bound='100*abs(KA-KB)/(KA+KB) percent for one prediction from identical summaries, over these two candidates.',
      limitations=['symmetry-reduced trapped homogeneous region; no exterior matching or event-horizon construction',
      'classical collisionless particles; no collisions, quantum fluctuations, two-particle correlations, rotation, evaporation',
      'new combinations of omitted moments, not a new physical law or claim of literature novelty',
      'each level constructs a DIFFERENT pair; the three gaps are not an accuracy comparison on one pair',
      'finite interval only; no inference of a universal breakdown curvature or of singularity avoidance'])
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2,default=dump))
    jets=CurvatureJets(export=out/'symbolic_curvature_jets.txt')
    primary=[];arrays={'time':times};all_solutions={};all_pairs={}
    for level in (1,2,3):
        pair=matched_pair(basis,level);all_pairs[level]=pair
        obs,sols,gap=pair_run(basis,pair,times);all_solutions[level]=sols
        jp=jets.evaluate(basis,pair['plus'],basis.initial_geometry(pair['plus']))
        jm=jets.evaluate(basis,pair['minus'],basis.initial_geometry(pair['minus']))
        preserved=np.max(abs(jp[:level]-jm[:level])/np.maximum(1,abs(jp[:level])))
        assert preserved<1e-11
        assert abs(jp[level]-jm[level])>1e-4
        row=dict(retained_through_M=level,initial_K=float(obs[0]['K'][0]),
                 initial_geometry=basis.initial_geometry(pair['plus']),
                 initial_moments_plus=pair['initial_moments_plus'],initial_moments_minus=pair['initial_moments_minus'],
                 weights_plus=pair['plus'],weights_minus=pair['minus'],
                 initial_jet_plus=jp,initial_jet_minus=jm,initial_jet_difference=jp-jm,
                 largest_relative_matching_residual=pair['relative_match_residual'],
                 matching_K_derivative_orders=list(range(level)),first_different_K_derivative_order=level,
                 endpoint_K_plus=float(obs[0]['K'][-1]),endpoint_K_minus=float(obs[1]['K'][-1]),
                 endpoint_pair_separation_percent=float(gap[-1]),
                 endpoint_two_candidate_minimax_lower_bound_percent=float(gap[-1]/2),
                 maximum_scaled_constraint_residual=max(float(np.max(abs(o['scaled_constraint']))) for o in obs),
                 all_spheres_future_trapped=bool(all(np.all(o['Hb']<0) for o in obs)),
                 minimum_b=min(float(o['b'].min()) for o in obs),
                 pressure_derivative_gap_at_start=float(obs[0]['dpr'][0]-obs[1]['dpr'][0]))
        cross=np.flatnonzero(gap/2>.1)
        row['first_sample_two_candidate_lower_bound_above_0p1percent']=None if len(cross)==0 else dict(time=float(times[cross[0]]),K_growth_plus=float(obs[0]['K'][cross[0]]/obs[0]['K'][0]))
        primary.append(row)
        arrays[f'level{level}_gap_percent']=gap
        for name,o in zip(('plus','minus'),obs):
            for key,value in o.items():arrays[f'level{level}_{name}_{key}']=value
        print('primary level',level,'pair separation',row['endpoint_pair_separation_percent'],flush=True)
    # Numerical controls on identical initial distributions, not retuned pairs.
    numerical=[]
    fine=MomentumBasis(centers,radial_order=24,angular_order=64)
    finer=MomentumBasis(centers,radial_order=32,angular_order=96)
    shorttimes=np.linspace(0.,.45,61)
    for level in (1,2,3):
        pair=all_pairs[level]
        for key,oldsol in zip(('plus','minus'),all_solutions[level]):
            w=pair[key];y0=basis.initial_geometry(w)
            tight=basis.evolve(w,0.,.45,y0,rtol=5e-13,atol=5e-15)
            alt=basis.evolve(w,0.,.45,y0,rtol=3e-12,atol=3e-14,method='RK45')
            sf=fine.evolve(w,0.,.45,y0,rtol=3e-12,atol=3e-14)
            sf2=finer.evolve(w,0.,.45,y0,rtol=3e-12,atol=3e-14)
            old=sample(basis,oldsol,w,shorttimes)['K']
            vals=[sample(basis,tight,w,shorttimes)['K'],sample(basis,alt,w,shorttimes)['K'],
                  sample(fine,sf,w,shorttimes)['K'],sample(finer,sf2,w,shorttimes)['K']]
            numerical.append(dict(level=level,side=key,
                  tighter_tolerance_max_relative_difference=float(np.max(abs(vals[0]/old-1))),
                  RK45_vs_DOP853_max_relative_difference=float(np.max(abs(vals[1]/old-1))),
                  fine_quadrature_max_relative_difference=float(np.max(abs(vals[2]/old-1))),
                  finer_vs_fine_quadrature_max_relative_difference=float(np.max(abs(vals[3]/vals[2]-1)))))
    # Verify moment derivatives with an independent complex-step differentiation.
    y=np.array([.1,-.2,.4,-1.5]);w=all_pairs[1]['plus'];dy=basis.rhs(0.,y,w)
    cs=basis.tensors(y+1j*1e-25*dy,w,4);c=basis.tensors(y,w,4)
    residual=[]
    for k in range(4):
        for i in range(k+1):
            j=k-i;expected=-((1+2*i)*y[2]+(2+2*j)*y[3])*c[i,j]+(2*k-1)*(y[2]*c[i+1,j]+y[3]*c[i,j+1])
            numeric=np.imag(cs[i,j])/1e-25
            residual.append(float(abs(numeric-expected)/max(1.,abs(expected))))
    assert max(residual)<1e-11
    # Vacuum control has the exact Schwarzschild interior relation K=48 M^2/b^6.
    zero=np.zeros(len(centers));vac=basis.evolve(zero,0.,.45)
    vo=sample(basis,vac,zero,shorttimes)
    exact_K=48/vo['b']**6
    mass=vo['b']/2*(1+(vo['Hb']*vo['b'])**2)
    vacuum=dict(max_relative_K_difference=float(np.max(abs(vo['K']/exact_K-1))),
                max_mass_invariant_error=float(np.max(abs(mass-1))),
                max_metric_relation_error=float(np.max(abs(vo['a']**2/(2/vo['b']-1)-1))))
    # Massless factorization control: different radial spectra, same energy -> same gravity
    # for this initially isotropic, collisionless, homogeneous setup.
    massless=MomentumBasis(centers,mass=0.)
    mp=massless_pair(massless)
    mobs,_,mgap=pair_run(massless,mp,shorttimes)
    massless_result=dict(max_pair_separation_percent=float(np.max(mgap)),
        weights_plus=mp['plus'],weights_minus=mp['minus'],
        weights_l1_difference=float(np.sum(abs(mp['plus']-mp['minus']))),
        explanation='For massless initially isotropic spectra, E(q,mu,a,b)=q*D(mu,a,b); all stress integrals factor through the same radial energy integral.')
    # Additional constructed pairs. No post-selection by future separation.
    cohort=[]
    for level in (1,2,3):
        for seed in protocol['cohort_seeds']:
            pp=matched_pair(basis,level,seed=seed)
            oo,_,gg=pair_run(basis,pp,shorttimes)
            cohort.append(dict(level=level,seed=seed,
                               pair_separation_at_0p45_percent=float(gg[-1]),
                               matching_residual=pp['relative_match_residual'],
                               max_scaled_constraint=max(float(np.max(abs(o['scaled_constraint']))) for o in oo)))
    cohort_summary=[]
    for level in (1,2,3):
        values=np.array([r['pair_separation_at_0p45_percent'] for r in cohort if r['level']==level])
        cohort_summary.append(dict(level=level,n_pairs=len(values),min=float(values.min()),median=float(np.median(values)),max=float(values.max())))
    # Reuse the ORIGINAL polynomial predictor, with no re-fitting of global physics parameters.
    # Local coefficients are, as before, fitted to each past-only window.
    forecasts=[]
    for level in (1,2,3):
        for key,sol in zip(('plus','minus'),all_solutions[level]):
            for now in (.15,.25,.35,.4):
                histtime=np.linspace(now-.1,now,41)
                hist=sample(basis,sol,all_pairs[level][key],histtime)['K']
                pred=polynomial(histtime,hist,.025)
                truth=basis.observe(sol.sol(now+.025),all_pairs[level][key])['K']
                err={k:None if v is None else float(100*abs(v/truth-1)) for k,v in pred['predicted_K_by_degree'].items()}
                forecasts.append(dict(level=level,side=key,origin=now,horizon=.025,truth=float(truth),
                                      predicted=pred['predicted_K_by_degree'],absolute_relative_error_percent=err))
    np.savez_compressed(out/'trajectories.npz',**arrays)
    summary=dict(protocol=protocol,primary=primary,numerical_controls=numerical,
                 moment_hierarchy_complex_step_max_scaled_residual=max(residual),vacuum_control=vacuum,
                 massless_spectrum_control=massless_result,cohort_summary=cohort_summary,
                 polynomial_transfer_summary={degree:dict(n=24,max_error_percent=max(r['absolute_relative_error_percent'][degree] for r in forecasts),mean_error_percent=float(np.mean([r['absolute_relative_error_percent'][degree] for r in forecasts]))) for degree in ('1','2','3')},
                 total_runtime_seconds=time.perf_counter()-t0,
                 environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    (out/'results.json').write_text(json.dumps(summary,indent=2,default=dump))
    (out/'cohort.json').write_text(json.dumps(cohort,indent=2,default=dump))
    (out/'legacy_predictor_on_new_model.json').write_text(json.dumps(forecasts,indent=2,default=dump))
    print(json.dumps({k:summary[k] for k in ['cohort_summary','polynomial_transfer_summary','vacuum_control','moment_hierarchy_complex_step_max_scaled_residual','total_runtime_seconds']},indent=2),flush=True)
    print('Results saved to',out,flush=True)
if __name__=='__main__':main()
