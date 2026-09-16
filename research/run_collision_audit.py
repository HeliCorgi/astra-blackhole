"""Run the recorded finite-time RTA audit; source and record travel together."""
from __future__ import annotations
from pathlib import Path
import argparse, json, platform, hashlib
import numpy as np
import scipy
import sympy as sp
from collision_model import RelaxationModel, fluid_evolve, symmetric_gap
from angular_model import AngularBasis, radial_weights

ROOT=Path(__file__).resolve().parent


def exact_checks():
    """Symbolic identities, not Lean proofs or analytic error bounds."""
    ha,hb,q,nu=sp.symbols('Ha Hb q nu'); c=sp.symbols('c0 c2 c4 c6 c8')
    rho,pr=c[0],c[0]/3+2*c[1]/15;pt=(rho-pr)/2
    db=-(3*hb**2+q+pr)/2;da=-pt-db-ha**2-hb**2-ha*hb
    K=4*((da+ha**2)**2+2*(db+hb**2)**2+2*(ha*hb)**2+(hb**2+q)**2)
    mu=sp.Symbol('mu');P=[sp.legendre(ell,mu) for ell in range(0,9,2)]
    A=sp.Matrix([[sp.integrate((mu*(1-mu**2)*sp.diff(p,mu)-4*mu**2*p)*v,(mu,-1,1))*(2*ell+1)/2 for p in P] for ell,v in zip(range(0,9,2),P)])
    dc=-4*hb*sp.Matrix(c)+(ha-hb)*A*sp.Matrix(c)-nu*sp.Matrix([0,*c[1:]])
    variables=[ha,hb,q,*c];rhs=[da,db,-2*hb*q,*dc]
    jets=[sp.expand(K)]
    for _ in range(3):jets.append(sp.expand(sum(sp.diff(jets[-1],v)*r for v,r in zip(variables,rhs))))
    differences={}
    for ell in (4,6,8):
        base={ha:sp.Rational(3,5),hb:-1,q:1,**{v:0 for v in c},c[0]:sp.Rational(4,5)}
        plus={**base,c[ell//2]:sp.Rational(16,25)};minus={**base,c[ell//2]:-sp.Rational(16,25)}
        degree=ell//2-1
        ds=[sp.factor(j.subs(plus)-j.subs(minus)) for j in jets[:degree+1]]
        assert all(d==0 for d in ds[:-1]) and ds[-1]!=0 and not ds[-1].has(nu)
        differences[str(ell)]={'first_distinct_derivative':degree,'jet_gaps':[str(d) for d in ds], 'leading_gap_independent_of_rate':True}
    energy=sp.simplify(dc[0]+(ha+3*hb)*rho+(ha-hb)*pr)
    constraint=2*ha*hb+hb**2+q-rho
    dconstraint=sp.expand(sum(sp.diff(constraint,v)*r for v,r in zip(variables,rhs)))
    # Calculate exact propagated constraint rather than assume its coefficient.
    constraint_factor=sp.factor(dconstraint/constraint)
    assert energy==0 and sp.simplify(constraint_factor+ha+2*hb)==0
    return {'angular_transport_matrix':[[str(x) for x in row] for row in A.tolist()],
            'energy_continuity_residual':str(energy),'constraint_dot_over_constraint':str(constraint_factor),
            'initial_curvature_jet_gaps':differences,'status':'SymPy exact algebra; not Lean formalization'}


def run(out: Path):
    p=json.loads((ROOT/'collision_protocol.json').read_text());out.mkdir(parents=True,exist_ok=True)
    end=p['end']; times=np.linspace(0,end,101);mu=np.linspace(-1,1,513)
    ref_model=RelaxationModel(); fluid=fluid_evolve(end,density=p['density'],number=p['number'])
    fluid_obs=[ref_model.observe(np.r_[y,np.zeros(len(ref_model.ells)-1)]) for y in fluid.sol(times).T]
    fluid_K=np.array([o['K'] for o in fluid_obs]); rows=[];curves={'time':times.tolist(),'fluid_K':fluid_K.tolist(),'pairs':[]}
    for nu in p['rates']:
      for ell in p['ells']:
        model=RelaxationModel(p['lmax'],nu);outputs=[];main_K=[];diagnostics=[];comparisons=[];brightness_comparisons=[]
        for sign in (1,-1):
            amp=sign*p['amplitude'];y0=model.initial(ell,amp,p['density'],p['number'])
            sol=model.evolve(ell,amp,end,rtol=p['rtol'],atol=p['atol'],method=p['method'],y0=y0)
            ys=sol.sol(times).T;obs=[model.observe(y) for y in ys];ks=np.array([o['K'] for o in obs]);main_K.append(ks);outputs.append(obs[-1])
            sampled_min=min(float(np.min(model.brightness(y,mu))/y[5]) for y in ys)
            maxC=max(abs(o['scaled_constraint']) for o in obs)
            number_error=max(abs(o['comoving_number']/p['number']-1) for o in obs)
            assert sampled_min > -1e-8 and maxC<1e-8 and number_error<1e-8
            diagnostics.append({'sampled_minimum_I_over_rho':sampled_min,'max_scaled_constraint':maxC,'max_comoving_number_relative_error':number_error})
            # Same observable/horizon, two independent time solvers at higher angular resolution.
            for method in ('Radau','DOP853'):
                fine=RelaxationModel(160,nu);ss=fine.evolve(ell,amp,end,rtol=2e-12,atol=2e-14,method=method,y0=fine.initial(ell,amp,p['density'],p['number']))
                kk=np.array([fine.observe(y)['K'] for y in ss.sol(times).T])
                comparisons.append(float(np.max(abs(kk/ks-1))))
                if method == 'Radau':
                    brightness_comparisons.append(max(float(np.max(abs(fine.brightness(ss.sol(t),mu)-model.brightness(sol.sol(t),mu)))/sol.sol(t)[5]) for t in times[::5]))
        gap=symmetric_gap(outputs[0]['K'],outputs[1]['K']);resolution=10*max(*comparisons,100*np.finfo(float).eps)
        row={'rate':nu,'ell':ell,'end':end,'endpoints':outputs,
             'pair_K_gap_percent':100*gap,'pair_K_gap_resolved':bool(gap>resolution),
             'empirical_resolution_relative':resolution,
             'max_fluid_K_error_percent':100*max(abs(o['K']/fluid_K[-1]-1) for o in outputs),
             'comparison_max_K_relative_difference':max(comparisons),
             'sampled_I_refinement_max_difference_over_rho':max(brightness_comparisons),'diagnostics':diagnostics,
             'collision_only_survival_factor':float(np.exp(-nu*end))}
        rows.append(row);curves['pairs'].append({'rate':nu,'ell':ell,'plus_K':main_K[0].tolist(),'minus_K':main_K[1].tolist()})
        print(f"rate={nu:4g} P{ell}: pair={100*gap:.8g}%, fluid error={row['max_fluid_K_error_percent']:.6g}%, resolved={row['pair_K_gap_resolved']}",flush=True)
    collisionless=[]
    for ell in p['ells']:
      for sign in (1,-1):
        old=AngularBasis(np.geomspace(.05,8,10),ell,sign*.8,mass=0,angular_order=128,radial_order=16)
        w=radial_weights(old);ss=old.evolve(w,0,end,rtol=2e-12,atol=2e-14)
        kk=np.array([old.observe(y,w)['K'] for y in ss.sol(times).T])
        model=RelaxationModel(160,0);spc=model.evolve(ell,sign*.8,end,rtol=2e-12,atol=2e-14)
        k2=np.array([model.observe(y)['K'] for y in spc.sol(times).T])
        error=float(np.max(abs(k2/kk-1)));assert error<1e-7
        collisionless.append({'ell':ell,'sign':sign,'max_K_relative_difference':error})
    result={'protocol':p,'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__},
            'exact_checks':exact_checks(),'fluid_endpoint':fluid_obs[-1],'cases':rows,
            'collisionless_reference':collisionless,
            'note':'Numerical disagreement is NOT a certified error bound. Unresolved tiny pair gaps are not physical detections.',
            'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('collision_model.py','run_collision_audit.py','collision_protocol.json')}}
    (out/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    (out/'curves.json').write_text(json.dumps(curves)+'\n')
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,default=ROOT/'collision_results');ap.add_argument('--check-against',type=Path);args=ap.parse_args()
    actual=run(args.out)
    if args.check_against:
        saved=json.loads(args.check_against.read_text())
        assert saved['protocol']==actual['protocol']
        assert saved['exact_checks']==actual['exact_checks']
        for new,old in zip(actual['cases'],saved['cases'],strict=True):
            assert (new['rate'],new['ell'])==(old['rate'],old['ell'])
            for i in (0,1):
                for field in ('K','rho','pr','pt','a','b'):
                    np.testing.assert_allclose(new['endpoints'][i][field],old['endpoints'][i][field],rtol=2e-8,atol=2e-10)
        print('Saved endpoints and symbolic identities agree within the stated numerical regression tolerance.')
if __name__=='__main__':main()
