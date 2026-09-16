"""Unbounded-moment audit tests. Copyright 2026 HeliCorgi. Apache-2.0."""
from pathlib import Path
import sys,unittest
import numpy as np
import sympy as sp
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from wdw_model import Packet,WDWModel,classical_path
from wdw_tail_model import states_at,log_moment,continuum_packet,threshold_overlap,tail_density,classical_log_moments

class TestWDWTails(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=WDWModel(.2,-4,24,.04);cls.packet=Packet(.2,np.sqrt(.2))

    def test_exact_gaussian_exponential_moments(self):
        m=self.model;p=self.packet;v=p.initial(m.x)
        for q in (2,4,6):self.assertAlmostEqual(log_moment(m.x,abs(v)**2,q,0),q*q*p.sigma**2/2,places=12)

    def test_initial_state_is_not_roundtrip_polluted(self):
        v,_=states_at(self.model,self.packet,[0.,1.]);np.testing.assert_array_equal(v[:,0],self.packet.initial(self.model.x))

    def test_same_spectral_propagator_and_norm(self):
        _,c=states_at(self.model,self.packet,[0.,.5,2.,6.]);self.assertLess(c['max_norm_error'],1e-11);self.assertLess(c['max_direct_state_l2_difference'],1e-11)

    def test_roundtrip_negative_control(self):
        _,c=states_at(self.model,self.packet,[0.]);real=log_moment(self.model.x,abs(self.packet.initial(self.model.x))**2,6,0)
        bad=log_moment(self.model.x,abs(c['initial_roundtrip'])**2,6,0)
        self.assertGreater(bad-real,20.) # norm accuracy cannot protect unbounded moments

    def test_partial_integral_not_conditional_mean(self):
        x=np.array([-1.,0.,1.]);p=np.array([.2,.3,.5])
        self.assertAlmostEqual(np.exp(log_moment(x,p,2,0,0,0)),.2*np.exp(-2)+.15)

    def test_cutoff_and_cap_monotonicity(self):
        v,_=states_at(self.model,self.packet,[2.]);p=abs(v[:,0])**2
        for cap in (False,True):
            a=[log_moment(self.model.x,p,6,2,cutoff=c,cap=cap) for c in (4.,8.,12.,16.)]
            self.assertTrue(np.all(np.diff(a)>0))

    def test_cap_is_bounded_and_dominates_partial(self):
        x=self.model.x;v,_=states_at(self.model,self.packet,[2.]);p=abs(v[:,0])**2
        for c in (4.,8.,12.):
            capped=log_moment(x,p,6,2,cutoff=c,cap=True)
            self.assertLessEqual(capped,6*c+1e-12)
            self.assertGreaterEqual(capped,log_moment(x,p,6,2,cutoff=c))

    def test_zero_moment_mass_returns_minus_infinity(self):
        self.assertEqual(log_moment(np.array([0.]),np.array([0.]),2,0),-np.inf)

    def test_invalid_inputs_rejected(self):
        for p in (np.array([-1.]),np.array([np.nan])):
            with self.assertRaises(ValueError):log_moment([0.],p,2,0)
        with self.assertRaises(ValueError):log_moment([0.],[1.],0,0)
        with self.assertRaises(ValueError):log_moment([0.],[1.],2,0,cap=True)
        with self.assertRaises(ValueError):states_at(self.model,Packet(.1,.4),[0.])

    def test_continuum_initial_normalization(self):
        p=self.packet;z=quad(lambda x:abs(continuum_packet(p,x))**2,-6,10,epsabs=1e-12)[0]
        self.assertAlmostEqual(z,1.,places=12)

    def test_threshold_overlap_refinement(self):
        a,err=threshold_overlap(self.packet);b,_=threshold_overlap(self.packet,14,1e-14)
        self.assertLess(abs(a-b),1e-11);self.assertGreater(abs(a),1000*err)

    def test_threshold_sine_integral_coefficient(self):
        s,x,t=sp.symbols('s x t',positive=True,real=True)
        f=2*s*x/(s*s+x*x)**2
        # Abel-regularized integral int_0^inf k exp(-s*k) sin(k*x) dk.
        self.assertEqual(sp.limit(x**3*f.subs(s,sp.I*t),x,sp.oo),2*sp.I*t)

    def test_tail_prediction_not_fitted(self):
        B,_=threshold_overlap(self.packet);x=np.array([20.,40.]);a=tail_density(x,2,.2,B);X=x+np.log(.4)-np.euler_gamma
        np.testing.assert_allclose(a*X**6,16*4*abs(B)**2/np.pi**2,rtol=1e-14)

    def test_classical_speed_bound_and_finite_moments(self):
        p=self.packet;ts=np.array([0.,.5,2.,6.]);z=classical_log_moments(p,ts,6,96)
        initial=18*p.sigma**2;self.assertTrue(np.all(z>=initial-1e-12));self.assertTrue(np.all(z<=initial+12*ts+1e-12))

    def test_small_probability_can_dominate_weighted_mean(self):
        x=np.array([0.,10.]);prob=np.array([1.-1e-12,1e-12]);total=np.exp(log_moment(x,prob,6,0,0))
        self.assertGreater(total,1e13);self.assertLess(prob[1],1e-10)

if __name__=='__main__':unittest.main()
