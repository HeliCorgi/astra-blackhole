"""Independent checks for the three zero-tail closures.
Copyright 2026 HeliCorgi. Apache-2.0.
"""
from pathlib import Path
import sys,unittest
import numpy as np
import sympy as sp
from numpy.polynomial import legendre as leg
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'research'))
from closure_model import ClosureModel,initial_state,moment_diagnostics,observe_series
from collision_model import RelaxationModel,fluid_evolve

class TestClosureStates(unittest.TestCase):
    def test_transport_projection_independent(self):
        x=sp.Symbol('x'); A=ClosureModel(4,0).A
        for j,l in enumerate((0,2,4)):
            p=sp.legendre(l,x); ap=x*(1-x*x)*sp.diff(p,x)-4*x*x*p
            for i,k in enumerate((0,2,4)):
                value=(2*k+1)*sp.integrate(sp.legendre(k,x)*ap,(x,-1,1))/2
                self.assertAlmostEqual(A[i,j],float(value),places=14)

    def test_low_orders_match_full_rhs(self):
        for order in (0,2,4):
            m=ClosureModel(order,4); h=RelaxationModel(16,4)
            y=initial_state(order,.8,-1.,{'4':.4})
            if order>=2:y[6]=.03
            z=np.r_[y,np.zeros((5+len(h.ells))-len(y))]
            np.testing.assert_allclose(m.rhs(0,y),h.rhs(0,z)[:len(y)],rtol=1e-14,atol=1e-14)

    def test_initial_summary_is_matched(self):
        for density,hb in ((.8,-1.),(.5,-.85),(1.1,-1.15)):
            states=[initial_state(l,density,hb,{'4':.6}) for l in (0,2,4,128)]
            for y in states:
                np.testing.assert_allclose(y[:6],states[-1][:6],rtol=0,atol=0)
                self.assertLess(y[2]+2*y[3],0)
                self.assertAlmostEqual(2*y[2]*y[3]+y[3]**2+1,y[5])

    def test_energy_conservation_identity(self):
        for order in (0,2,4):
            m=ClosureModel(order,16); y=initial_state(order,.8,-1.,{'4':.3})
            if order>=2:y[6]=-.1
            rho=y[5];pr=rho/3+(2*y[6]/15 if order>=2 else 0);pt=(rho-pr)/2
            expected=-(y[2]+2*y[3])*rho-y[2]*pr-2*y[3]*pt
            self.assertAlmostEqual(m.rhs(0,y)[5],expected,places=14)

    def test_fluid_limit_is_independent_reference(self):
        m=ClosureModel(0,0); y=initial_state(0,.8,-1.,{})
        sol=m.evolve(y,.2,rtol=1e-12,atol=1e-14); ref=fluid_evolve(end=.2)
        np.testing.assert_allclose(sol.sol(.2),ref.sol(.2),rtol=1e-9,atol=1e-11)

    def test_negative_polynomial_can_have_feasible_moments(self):
        d=moment_diagnostics(np.array([1.,-1.2]))
        self.assertGreater(d['moment_margin'],0)
        self.assertLess(d['polynomial_minimum'],0)

    def test_infeasible_fourth_moment_is_rejected(self):
        d=moment_diagnostics(np.array([1.,0.,-4.]))
        self.assertLess(d['moment_margin'],0)

    def test_positive_measure_construction(self):
        m1,m2=.3,.2
        c2=15*(m1-1/3)/2; c4=315*(m2-1/5-4*c2/35)/8
        d=moment_diagnostics(np.array([1.,c2,c4]))
        self.assertGreater(d['moment_margin'],0)
        location=m2/m1;probability=m1*m1/m2
        self.assertAlmostEqual(probability*location,d['m1'])
        self.assertAlmostEqual(probability*location**2,d['m2'])
        self.assertTrue(0<probability<1 and 0<location<1)

    def test_polynomial_minimum_with_stationary_points(self):
        for c in [np.array([1.,-.6,.3]),np.array([1.,.8,-.4]),np.array([1.,-1.2])]:
            coef=np.zeros(2*len(c)-1);coef[::2]=c
            grid=float(leg.legval(np.linspace(-1,1,10001),coef).min())
            exact=moment_diagnostics(c)['polynomial_minimum']
            self.assertLessEqual(exact,grid+1e-12)
            self.assertLess(grid-exact,1e-6)

    def test_evolution_constraints_and_method_comparison(self):
        for order in (0,2,4):
            m=ClosureModel(order,16); y=initial_state(order,.8,-1.,{'4':.6})
            s=m.evolve(y,.3);r=m.evolve(y,.3,method='DOP853',rtol=2e-12,atol=2e-14)
            t=np.linspace(0,.3,21);a=observe_series(m,s,t);b=observe_series(m,r,t)
            self.assertLess(np.max(abs(a['K']/b['K']-1)),1e-8)
            self.assertLess(np.max(abs(a['constraint'])),1e-8)
            self.assertLess(np.max(abs(a['number_drift'])),1e-8)

    def test_constraint_propagation_symbolically(self):
        ha,hb,q,rho,pr,pt=sp.symbols('ha hb q rho pr pt')
        db=-(3*hb**2+q+pr)/2; da=-pt-db-ha**2-hb**2-ha*hb
        dr=-(ha+2*hb)*rho-ha*pr-2*hb*pt
        C=2*ha*hb+hb**2+q-rho
        dC=2*da*hb+2*ha*db+2*hb*db-2*hb*q-dr
        self.assertEqual(sp.factor(dC+(ha+2*hb)*C),0)

    def test_invalid_inputs(self):
        for order in (1,3,6):
            with self.assertRaises(ValueError):ClosureModel(order,1)
        with self.assertRaises(ValueError):ClosureModel(2,-1)
        with self.assertRaises(ValueError):initial_state(4,.8,-1.,{'4':1.})
        with self.assertRaises(ValueError):initial_state(4,.8,0.,{})
        with self.assertRaises(ValueError):moment_diagnostics(np.array([0.,1.]))
        with self.assertRaises(ValueError):ClosureModel(2,1).evolve(np.ones(3),.1)

if __name__=='__main__':unittest.main()
