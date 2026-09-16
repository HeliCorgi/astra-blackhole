"""Finite quantum geometry checks, not proofs of quantum gravity. Apache-2.0."""
from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'research'))
from wdw_model import Packet,WDWModel,classical_path,classical_ensemble

class TestWDWGeometry(unittest.TestCase):
    def test_classical_hamilton_equations(self):
        t=np.linspace(0,6,30);dt=1e-5;x,p=classical_path(2.,-1.,t);E=np.sqrt(1+np.exp(-4))
        xp,pp=classical_path(2.,-1.,t+dt);xm,pm=classical_path(2.,-1.,t-dt)
        np.testing.assert_allclose((xp-xm)/(2*dt),p/E,atol=1e-9)
        np.testing.assert_allclose((pp-pm)/(2*dt),np.exp(-2*x)/E,atol=1e-9)
    def test_classical_constraint_and_radius(self):
        t=np.linspace(0,10,50);x,p=classical_path(2.,-1.,t)
        np.testing.assert_allclose(p*p+np.exp(-2*x),1+np.exp(-4),atol=1e-13)
        self.assertTrue(np.all(np.diff(2-x-t)<0))
    def test_x_reflection_is_not_radius_bounce(self):
        turn=-np.arcsinh(-np.exp(2.));dt=1e-5
        x,p=classical_path(2.,-1.,np.array([turn-dt,turn,turn+dt]))
        self.assertLess(p[0],0);self.assertGreater(p[2],0)
        self.assertAlmostEqual((x[-1]-x[0])/(2*dt),0,places=8)
        self.assertAlmostEqual(-1-(x[-1]-x[0])/(2*dt),-1,places=8)
    def test_schwarzschild_radius_relation(self):
        t=np.linspace(0,7,21);turn=np.arcsinh(np.exp(2));E=np.sqrt(1+np.exp(-4))
        x,_=classical_path(2.,-1.,t);r=np.exp(-x-t)/4;M=E*np.exp(-turn)/4
        np.testing.assert_allclose(r/(2*M),1/(1+np.exp(2*(t-turn))),rtol=2e-14)
    def test_positive_symmetric_discrete_operator(self):
        m=WDWModel(.2,-2,4,.05);n=len(m.x)
        D=(2*np.eye(n)-np.eye(n,k=1)-np.eye(n,k=-1))/.05**2
        A=.2**2*(D+.05**2*D@D/12)+np.diag(np.exp(-2*m.x))
        np.testing.assert_allclose(m.operator(np.eye(n)),A,rtol=1e-14,atol=1e-12)
        self.assertGreater(m.eigenvalues.min(),0)
    def test_unitarity_energy_and_wdw_constraint(self):
        m=WDWModel(.2,dx=.04);o=m.evolve(Packet(.2,.4),np.linspace(0,6,13))
        self.assertLess(o['norm_error'],1e-11);self.assertLess(o['energy_relative_drift'],1e-11)
        self.assertLess(o['discrete_constraint_relative_residual'],1e-8)
    def test_reversibility(self):
        m=WDWModel(.2,dx=.04);p=Packet(.2,.4);v=p.initial(m.x);a=m.U.T@v
        v2=m.U@(a*np.exp(-1j*m.energy/.2));v3=m.U@((m.U.T@v2)*np.exp(1j*m.energy/.2))
        np.testing.assert_allclose(v,v3,atol=2e-12)
    def test_gaussian_initial_moments(self):
        m=WDWModel(.2,dx=.01);p=Packet(.2,.4,.5);o=m.evolve(p,[0])['initial'];target=p.continuum_moments()
        for k in target:self.assertAlmostEqual(o[k],target[k],delta=2e-6)
    def test_matched_marginals_but_opposite_covariance(self):
        a,b=Packet(.2,.4,-.5),Packet(.2,.4,.5);x=np.linspace(-4,16,2001)
        np.testing.assert_allclose(abs(a.initial(x))**2,abs(b.initial(x))**2,rtol=2e-14,atol=2e-17)
        aa,bb=a.continuum_moments(),b.continuum_moments()
        for k in ['mean_x','mean_p','var_x','var_p']:self.assertEqual(aa[k],bb[k])
        self.assertEqual(aa['covariance_xp'],-bb['covariance_xp'])
    def test_phase_changes_future_variance(self):
        m=WDWModel(.2,dx=.04)
        a=m.evolve(Packet(.2,.4,-.5),[0,6]);b=m.evolve(Packet(.2,.4,.5),[0,6])
        self.assertGreater(abs(a['var_x'][-1]-b['var_x'][-1]),.1)
    def test_ensemble_initial_moments(self):
        p=Packet(.2,.4,.5);e=classical_ensemble(p,np.array([0.,1.]),32)
        self.assertAlmostEqual(e['mean_x'][0],2.,places=13);self.assertAlmostEqual(e['var_x'][0],.16,places=12)
    def test_ensemble_resolution(self):
        p=Packet(.1,np.sqrt(.1));t=np.linspace(0,6,11)
        a=classical_ensemble(p,t,48);b=classical_ensemble(p,t,72)
        self.assertLess(max(abs(a['mean_x']-b['mean_x'])),1e-7)
    def test_short_box_is_not_certified_by_norm(self):
        p=Packet(.2,.4,.5);a=WDWModel(.2,dx=.04).evolve(p,[0,6]);b=WDWModel(.2,-4,4,.04).evolve(p,[0,6])
        self.assertLess(b['norm_error'],1e-11);self.assertGreater(abs(a['mean_x'][-1]-b['mean_x'][-1]),1e-3)
    def test_grid_refinement(self):
        p=Packet(.2,.4);o=[WDWModel(.2,dx=d).evolve(p,[0,3])['mean_x'][-1] for d in [.04,.02,.01]]
        self.assertLess(abs(o[1]-o[2]),abs(o[0]-o[1])/8)
    def test_invalid_inputs(self):
        for h in (0,-1,np.nan):
            with self.assertRaises(ValueError):Packet(h,.4)
        with self.assertRaises(ValueError):WDWModel(.2,right=-5)
        with self.assertRaises(ValueError):WDWModel(.2,dx=.013)
        with self.assertRaises(ValueError):WDWModel(.2,dx=.04).evolve(Packet(.1,.4),[0,1])
        with self.assertRaises(ValueError):classical_ensemble(Packet(.2,.4),[0],2)

if __name__=='__main__':unittest.main()
