"""Finite/model-internal checks. Copyright 2026 HeliCorgi. Apache-2.0."""
from pathlib import Path
import sys
import unittest
import numpy as np
from numpy.polynomial import legendre as leg
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'research'))
from positive_closure_model import AngularMaxEntropy,PositiveClosureModel,MomentDomainError
from closure_model import ClosureModel,initial_state,moment_diagnostics,observe_series
from collision_model import RelaxationModel

class PositiveClosureTests(unittest.TestCase):
    def test_isotropic(self):
        for L in (2,4):
            engine=AngularMaxEntropy(L)
            fit=engine.fit(np.r_[2.,np.zeros(L//2)])
            self.assertLess(np.max(abs(fit['multipliers'])),1e-11)
            self.assertLess(abs(fit['tail']),1e-11)
            self.assertAlmostEqual(engine.minimum_log_density(fit),0.,places=11)

    def test_manufactured_moment_inversion(self):
        for L,lam in [(2,[-3.]),(2,[5.]),(4,[-2.,1.]),(4,[4.,-3.]),(4,[-5.,4.])]:
            engine=AngularMaxEntropy(L)
            mean,_,_,_=engine._statistics(np.array(lam))
            c=np.r_[1.,(2*engine.degrees+1)*mean]
            fit=engine.fit(c)
            np.testing.assert_allclose(fit['multipliers'],lam,rtol=0,atol=1e-9)
            self.assertLess(engine.independent_residual(c,fit)['moment_residual'],1e-10)

    def test_independent_adaptive_integral(self):
        for L,c in [(2,[1.,-1.8]),(4,[1.,0.,3.])]:
            engine=AngularMaxEntropy(L);c=np.array(c);fit=engine.fit(c)
            for ell in range(0,L+3,2):
                integral=quad(lambda x: float(np.exp(engine.log_density(fit,np.array([x])))[0])*leg.legval(x,[0.]*ell+[1.])/2,
                              -1.,1.,epsabs=1e-12,epsrel=1e-12)[0]*(2*ell+1)
                target=c[ell//2]/c[0] if ell<=L else fit['tail']/c[0]
                self.assertAlmostEqual(integral,target,places=10)

    def test_positive_reconstruction_without_extra_moments(self):
        for L,c in [(2,[1.,-1.8]),(4,[1.,0.,3.])]:
            c=np.array(c);self.assertLess(moment_diagnostics(c)['polynomial_minimum'],0.)
            self.assertGreater(moment_diagnostics(c)['moment_margin'],0.)
            engine=AngularMaxEntropy(L);fit=engine.fit(c)
            self.assertGreater(np.exp(engine.minimum_log_density(fit)),0.)
            self.assertLess(engine.independent_residual(c,fit)['moment_residual'],1e-9)
            self.assertEqual(len(fit['multipliers']),L//2)

    def test_infeasible_and_boundary_rejected_not_clipped(self):
        for L,c in [(2,[1.,6.]),(2,[1.,-2.5]),(4,[1.,0.,20.]),(4,[1.,-2.,-3.]),(2,[-1.,0.])]:
            with self.assertRaises(MomentDomainError):AngularMaxEntropy(L).fit(c)

    def test_invalid_configuration_and_nonfinite(self):
        for kw in [dict(lmax=6),dict(lmax=2,quadrature_order=8),dict(lmax=2,tolerance=0),dict(lmax=2,max_iterations=0)]:
            with self.assertRaises(ValueError):AngularMaxEntropy(**kw)
        with self.assertRaises(ValueError):PositiveClosureModel(2,-1)
        with self.assertRaises(ValueError):AngularMaxEntropy(2).fit([1.,np.nan])

    def test_explicit_nonconvergence_is_failure(self):
        with self.assertRaises(RuntimeError):AngularMaxEntropy(4,max_iterations=1).fit([1.,0.,3.])

    def test_deterministic_and_scale_invariant(self):
        e=AngularMaxEntropy(4);c=np.array([.8,-.4,.2])
        a=e.fit(c);e.fit([1.,0.,3.]);b=e.fit(7*c)
        np.testing.assert_allclose(a['multipliers'],b['multipliers'],rtol=0,atol=1e-12)
        self.assertAlmostEqual(a['tail']*7,b['tail'],places=10)

    def test_algebraic_tail_changes_only_highest_retained_row(self):
        for L in (2,4):
            m=PositiveClosureModel(L,4);z=ClosureModel(L,4)
            y=initial_state(L,.8,-1,{'4':.6});y[6]=-.7
            fit=m.reconstruction.fit(y[5:])
            difference=m.rhs(0,y)-z.rhs(0,y)
            np.testing.assert_allclose(difference[:-1],0.,atol=1e-14)
            self.assertAlmostEqual(difference[-1],(y[2]-y[3])*m.tail_coupling[-1]*fit['tail'],places=12)
            self.assertEqual(len(difference),len(y))

    def test_independent_weak_transport_projection(self):
        # Integration by parts gives the adjoint operator
        # -mu(1-mu^2) P_l' - (1+mu^2)P_l.
        for L in (2,4):
            m=PositiveClosureModel(L,3);y=initial_state(L,.8,-1,{'4':.6});y[6]=-.3
            e=m.reconstruction;fit=e.fit(y[5:]);mu,w=leg.leggauss(160)
            I=y[5]*np.exp(e.log_density(fit,mu));dc=[]
            for j,ell in enumerate(range(0,L+1,2)):
                p=np.r_[np.zeros(ell),1.];P=leg.legval(mu,p);dP=leg.legval(mu,leg.legder(p))
                adj=-mu*(1-mu**2)*dP-(1+mu**2)*P
                v=-4*y[3]*y[5+j]+(y[2]-y[3])*(2*ell+1)*np.sum(w/2*I*adj)
                if ell:v-=m.rate*y[5+j]
                dc.append(v)
            np.testing.assert_allclose(m.rhs(0,y)[5:],dc,rtol=0,atol=2e-10)

    def test_energy_and_number_conservation_equations(self):
        for L in (2,4):
            m=PositiveClosureModel(L,16);y=initial_state(L,.8,-1,{'4':.6});y[6]=-.3
            dy=m.rhs(0,y);ha,hb=y[2:4];rho,c2=y[5:7]
            pr=rho/3+2*c2/15;pt=(rho-pr)/2
            self.assertAlmostEqual(dy[5],-(ha+2*hb)*rho-ha*pr-2*hb*pt,places=12)
            self.assertAlmostEqual(dy[4],-(ha+2*hb)*y[4],places=12)

    def test_collision_terms_preserve_energy(self):
        for L in (2,4):
            a=PositiveClosureModel(L,0);b=PositiveClosureModel(L,7)
            y=initial_state(L,.8,-1,{'4':.6});y[6]=-.3
            diff=b.rhs(0,y)-a.rhs(0,y)
            np.testing.assert_allclose(diff[:6],0.,atol=1e-14)
            np.testing.assert_allclose(diff[6:],-7*y[6:],atol=1e-14)

    def test_hessian_is_positive_on_examples(self):
        for L in (2,4):
            e=AngularMaxEntropy(L)
            for x in [np.zeros(L//2),np.full(L//2,-3),np.full(L//2,4)]:
                _,H,_,_=e._statistics(x)
                self.assertGreater(np.linalg.eigvalsh(H).min(),0.)

    def test_minimum_exponent_and_brightness(self):
        e=AngularMaxEntropy(4);fit=e.fit(np.array([1.,0.,3.]))
        sampled=min(e.log_density(fit,np.linspace(-1,1,10001)))
        exact=e.minimum_log_density(fit)
        self.assertLessEqual(exact,sampled+1e-12)
        self.assertLess(sampled-exact,1e-6)

    def test_short_evolution_constraints_and_solver_agreement(self):
        for L in (2,4):
            m=PositiveClosureModel(L,16);y=initial_state(L,.8,-1,{'4':.6})
            s=m.evolve(y,.05);o=observe_series(m,s,np.linspace(0,.05,11))
            v=PositiveClosureModel(L,16,160,2e-14)
            sv=v.evolve(y,.05,method='DOP853',rtol=2e-12,atol=2e-14)
            ov=observe_series(v,sv,np.linspace(0,.05,11))
            self.assertLess(max(abs(o['constraint'])),1e-8)
            self.assertLess(max(abs(o['number_drift'])),1e-8)
            self.assertLess(max(abs(o['K']/ov['K']-1)),1e-8)

if __name__=='__main__':unittest.main()
