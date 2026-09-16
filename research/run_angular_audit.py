"""Reproduce angular-state audit. Paired separations are not forecast errors."""
from __future__ import annotations
import argparse, json, platform, time
from pathlib import Path
import numpy as np
import scipy
import sympy as sp
from angular_model import AngularBasis, radial_weights, general_jets, angular_identities
from kinetic_model import MomentumBasis, massless_pair
from local_jets import CurvatureJets
ROOT = Path(__file__).resolve().parent
CENTERS = np.geomspace(.05, 8., 10)


def sample(basis, sol, weights, times):
    obs = [basis.observe(y, weights) for y in sol.sol(times).T]
    return {k: np.array([o[k] for o in obs]) for k in obs[0]}


def pair(mass, ell, alpha, times, jets=None, radial_order=12, angular_order=32):
    bases = [AngularBasis(CENTERS, ell, sign*alpha, mass=mass,
              radial_order=radial_order, angular_order=angular_order) for sign in (1,-1)]
    w = radial_weights(bases[0]); y0 = np.array([0.,0.,.6,-1.])
    sols = [b.evolve(w, 0., float(times[-1]), y0) for b in bases]
    obs = [sample(b,s,w,times) for b,s in zip(bases,sols)]
    p,m = obs
    gap = 200*abs(p['K']-m['K'])/(p['K']+m['K'])
    c = [b.tensors(y0,w,4) for b in bases]
    retained = [(i,j) for i,j in c[0] if 2*(i+j)<ell]
    row = dict(mass=mass,ell=ell,alpha=alpha,
      initial_K=float(p['K'][0]),initial_n=float(p['n'][0]),
      initial_pressure=float(p['pr'][0]),
      matched_moment_max_absolute_residual=max(abs(c[0][ij]-c[1][ij]) for ij in retained),
      initial_dpr_gap=float(p['dpr'][0]-m['dpr'][0]),
      endpoint_K_plus=float(p['K'][-1]),endpoint_K_minus=float(m['K'][-1]),
      endpoint_pair_separation_percent=float(gap[-1]),
      maximum_scaled_constraint=max(float(np.max(abs(o['scaled_constraint']))) for o in obs),
      minimum_b=min(float(o['b'].min()) for o in obs),
      all_spheres_future_trapped=all(bool(np.all(o['Hb']<0)) for o in obs))
    if jets is not None:
        jp,jm=[general_jets(jets,b,w,y0) for b in bases]
        first=ell//2-1
        row.update(initial_K_jets_plus=jp.tolist(),initial_K_jets_minus=jm.tolist(),
          initial_K_jet_gaps=(jp-jm).tolist(),first_different_K_derivative_order=first)
        assert np.max(abs(jp[:first]-jm[:first])/np.maximum(1,abs(jp[:first])))<1e-11
        assert abs(jp[first]-jm[first])>1e-5
    assert row['matched_moment_max_absolute_residual']<1e-12
    assert row['maximum_scaled_constraint']<1e-8
    assert row['all_spheres_future_trapped']
    return row,bases,w,y0,sols,obs,gap


def exact_checks(jets):
    """Finite polynomial identities and massless jets in exact rational arithmetic."""
    mu=sp.Symbol('mu'); alpha=sp.Rational(4,5); rho=sp.Rational(4,5)
    output=[]
    for ell in (4,6,8):
        identities=angular_identities(ell)
        values=[]
        for sign in (1,-1):
            d={jets.ha:sp.Rational(3,5),jets.hb:-1,jets.B:1}
            for (i,j),symbol in jets.C.items():
                d[symbol]=rho*sp.integrate(mu**(2*i)*(1-mu**2)**j*
                     (1+sign*alpha*sp.legendre(ell,mu)),(mu,-1,1))/2
            values.append([sp.factor(e.subs(d)) for e in jets.expressions])
        gap=[sp.factor(p-m) for p,m in zip(*values)]
        first=ell//2-1
        assert all(g==0 for g in gap[:first]) and gap[first]!=0
        output.append(dict(ell=ell,angular_integrals=identities,
           exact_massless_K_jet_gaps=[str(g) for g in gap],
           first_different_order=first))
    return output


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/'angular_results')
    ap.add_argument('--check-against',type=Path,help='Check new endpoints and exact identities against a saved baseline.')
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    reference=json.loads(args.check_against.read_text()) if args.check_against else None
    started=time.perf_counter();times=np.linspace(0.,.45,91)
    protocol=json.loads((ROOT/'angular_protocol.json').read_text())
    jets=CurvatureJets();exact=exact_checks(jets)
    primary=[];numerical=[];curves={'time':times.tolist()}
    for mass in (1.,0.):
        for ell in (4,6,8):
            row,bases,w,y0,sols,obs,gap=pair(mass,ell,.8,times,jets)
            primary.append(row);curves[f'm{int(mass)}_ell{ell}_separation_percent']=gap.tolist()
            # Compare the SAME distributions and initial geometry, not re-fitted pairs.
            for side,(b,s,o) in enumerate(zip(bases,sols,obs)):
                variants=[('tighter',b,dict(rtol=5e-13,atol=5e-15)),
                          ('RK45',b,dict(method='RK45',rtol=3e-12,atol=3e-14)),
                          ('quadrature',AngularBasis(CENTERS,ell,(1-2*side)*.8,
                             mass=mass,radial_order=24,angular_order=64),
                           dict(rtol=3e-12,atol=3e-14))]
                for name,bb,options in variants:
                    ss=bb.evolve(w,0.,.45,y0,**options)
                    kk=sample(bb,ss,w,times)['K']
                    error=float(np.max(abs(kk/o['K']-1)))
                    assert error<1e-7
                    numerical.append(dict(mass=mass,ell=ell,side=side,
                                          variant=name,max_relative_K_difference=error))
            print('mass',mass,'ell',ell,'separation %',row['endpoint_pair_separation_percent'],flush=True)
    scan=[]
    for mass in (1.,0.):
        for ell in (4,6,8):
            for alpha in (.2,.5):
                scan.append(pair(mass,ell,alpha,times)[0])
    # Kinematic control: imposed equal scale factors, NOT a coupled Einstein solution.
    isotropic=[];hierarchy=[]
    for mass in (1.,0.):
        bp,bm=[AngularBasis(CENTERS,4,a,mass=mass) for a in (.8,-.8)]
        w=radial_weights(bp);y=np.array([np.log(.8),np.log(.8),-.7,-.7])
        op,om=[b.observe(y,w) for b in (bp,bm)]
        isotropic.append(max(float(abs(op[k]-om[k])) for k in ('rho','pr','pt','dpr','dpt')))
        # Complex-step derivative of the explicit momentum integrals.
        y=np.array([.1,-.2,.4,-1.5]);dy=bp.rhs(0.,y,w)
        c=bp.tensors(y,w,4);cs=bp.tensors(y+1j*1e-25*dy,w,4)
        for k in range(4):
            for i in range(k+1):
                j=k-i
                rhs=-((1+2*i)*y[2]+(2+2*j)*y[3])*c[i,j]+(2*k-1)*(y[2]*c[i+1,j]+y[3]*c[i,j+1])
                hierarchy.append(float(abs(cs[i,j].imag/1e-25-rhs)/max(1.,abs(rhs))))
    zero=pair(1.,4,0.,times)[-1]
    mb=MomentumBasis(CENTERS,mass=0.);mp=massless_pair(mb)
    mo=[sample(mb,mb.evolve(mp[s],0.,.45),mp[s],times) for s in ('plus','minus')]
    radial_control=float(np.max(200*abs(mo[0]['K']-mo[1]['K'])/(mo[0]['K']+mo[1]['K'])))
    controls=dict(zero_amplitude_max_separation_percent=float(np.max(zero)),
      imposed_isotropic_scaling_max_stress_response_gap=max(isotropic),
      complex_step_moment_hierarchy_max_scaled_residual=max(hierarchy),
      massless_radial_spectrum_max_separation_percent=radial_control)
    assert controls['zero_amplitude_max_separation_percent']<1e-10
    assert max(isotropic)<1e-12 and max(hierarchy)<1e-11 and radial_control<1e-8
    results=dict(protocol=protocol,primary=primary,amplitude_scan=scan,
       exact_checks=exact,numerical_controls=numerical,negative_controls=controls,
       interpretation='Model-internal pair differences, not AI errors or observations; no singularity resolution.',
       environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),
       runtime_seconds=time.perf_counter()-started)
    if reference is not None:
        assert results['protocol']==reference['protocol'], 'Protocol changed; review the baseline explicitly.'
        assert results['exact_checks']==reference['exact_checks']
        assert len(results['primary'])==len(reference['primary'])
        for new,old in zip(results['primary'],reference['primary']):
            assert (new['mass'],new['ell'])==(old['mass'],old['ell'])
            for key in ('endpoint_K_plus','endpoint_K_minus','endpoint_pair_separation_percent'):
                np.testing.assert_allclose(new[key],old[key],rtol=1e-7,atol=1e-9)
        results['saved_baseline_comparison']='passed (rtol=1e-7, atol=1e-9; protocol and exact identities equal)'
    for name,obj in [('results.json',results),('curves.json',curves)]:
        (args.out/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
    print(json.dumps(controls,indent=2));print('Saved to',args.out)

if __name__=='__main__':main()
