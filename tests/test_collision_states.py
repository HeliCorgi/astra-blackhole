"""Checks for the restricted massless relaxation model; no empirical validation."""
from pathlib import Path
import sys, unittest
import numpy as np
from numpy.polynomial import legendre as leg
from scipy.integrate import quad, solve_ivp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'research'))
from collision_model import RelaxationModel, transport_matrix, fluid_evolve

class TestCollisionStates(unittest.TestCase):
    def test_transport_independent_projection(self):
        L=16;A=transport_matrix(L);x,w=leg.leggauss(64)
        for j,ell in enumerate(range(0,L+1,2)):
            p=np.eye(1,L+1,ell).ravel()
            rhs=x*(1-x*x)*leg.legval(x,leg.legder(p))-4*x*x*leg.legval(x,p)
            for i,k in enumerate(range(0,L+1,2)):
                v=leg.legval(x,np.eye(1,L+1,k).ravel())
                self.assertAlmostEqual(A[i,j],float(w@(rhs*v)*(2*k+1)/2),places=11)

    def test_initial_matching(self):
        for ell in (4,6,8):
            m=RelaxationModel(16,4)
            a,b=[m.initial(ell,s*.8) for s in (1,-1)]
            np.testing.assert_array_equal(a[:7],b[:7])
            oa,ob=m.observe(a),m.observe(b)
            for key in ('rho','pr','pt','n','K','Ha','Hb'):
                self.assertEqual(oa[key],ob[key])
            self.assertLess(abs(oa['scaled_constraint']),1e-14)
            self.assertGreater(np.min(m.brightness(a,np.linspace(-1,1,1001))),0)

    def test_collision_conserves_number_and_energy(self):
        free,coll=RelaxationModel(16,0),RelaxationModel(16,7)
        y=free.initial(); y[6]=.1;y[7]=.2
        change=coll.rhs(0,y)-free.rhs(0,y)
        np.testing.assert_allclose(change[:6],0,atol=1e-14)
        np.testing.assert_allclose(change[6:],-7*y[6:],atol=1e-14)
        # Negative control: indiscriminate damping would destroy energy.
        self.assertNotEqual(-7*y[5],0)

    def test_exact_continuity(self):
        m=RelaxationModel(16,16);y=m.initial();y[:4]=[.1,-.1,.3,-.6];y[6]=.05
        o=m.observe(y);d=m.rhs(0,y);theta=y[2]+2*y[3]
        self.assertAlmostEqual(d[5],-theta*o['rho']-y[2]*o['pr']-2*y[3]*o['pt'],places=13)
        self.assertAlmostEqual(d[4],-theta*y[4],places=13)

    def test_local_equilibrium_matching(self):
        for n,rho in ((.2,.8),(1.,.4),(.03,2.)):
            T=rho/(3*n);A=n/(2*T**3)
            nn=quad(lambda q:q*q*A*np.exp(-q/T),0,np.inf)[0]
            ee=quad(lambda q:q**3*A*np.exp(-q/T),0,np.inf)[0]
            self.assertAlmostEqual(nn,n,places=10);self.assertAlmostEqual(ee,rho,places=10)

    def test_collision_only_relaxation_and_entropy(self):
        m=RelaxationModel(16,3);c=m.initial()[5:]
        def rhs(t,x):return np.r_[0.,-m.rate*x[1:]]
        sol=solve_ivp(rhs,(0,.4),c,rtol=1e-11,atol=1e-13)
        expected=c.copy();expected[1:]*=np.exp(-1.2)
        np.testing.assert_allclose(sol.y[:,-1],expected,rtol=1e-10,atol=1e-12)
        x,w=leg.leggauss(128);p=leg.legval(x,[0,0,0,0,1]);f=1+.8*p
        entropy_production=float(w@((f-1)*np.log(f))/2)
        self.assertGreater(entropy_production,0)

    def test_isotropic_prescribed_geometry_control(self):
        # Not a KS Einstein solution: an explicitly prescribed isotropic transport control.
        m=RelaxationModel(16,4);c=m.initial(8)[5:];H=-.3
        def rhs(t,x):return -4*H*x+np.r_[0.,-m.rate*x[1:]]
        sol=solve_ivp(rhs,(0,.3),c,rtol=1e-11,atol=1e-13)
        expected=c*np.exp(-4*H*.3);expected[1:]*=np.exp(-4*.3)
        np.testing.assert_allclose(sol.y[:,-1],expected,rtol=1e-9,atol=1e-12)
        self.assertEqual(sol.y[1,-1],0) # no P2 pressure source from P8 without shear.

    def test_fluid_limit_and_number(self):
        f=fluid_evolve(.2);m=RelaxationModel(16);target=m.observe(np.r_[f.y[:,-1],np.zeros(8)])['K']
        errors=[]
        for rate in (16,256):
            mm=RelaxationModel(32,rate);ss=mm.evolve(end=.2)
            o=mm.observe(ss.y[:,-1]);errors.append(abs(o['K']/target-1))
            self.assertAlmostEqual(o['comoving_number'],.8/3,places=9)
        self.assertLess(errors[1],errors[0])

    def test_reject_low_angular_resolution_even_when_K_stable(self):
        coarse=RelaxationModel(32,0);fine=RelaxationModel(128,0)
        a=coarse.evolve(4,.8);b=fine.evolve(4,.8)
        self.assertLess(np.min(coarse.brightness(a.y[:,-1],np.linspace(-1,1,1025))),0)
        self.assertGreater(np.min(fine.brightness(b.y[:,-1],np.linspace(-1,1,1025))),0)
        self.assertLess(abs(coarse.observe(a.y[:,-1])['K']/fine.observe(b.y[:,-1])['K']-1),1e-8)

    def test_invalid_inputs(self):
        for rate in (-1,np.inf,np.nan):
            with self.assertRaises(ValueError):RelaxationModel(rate=rate)
        with self.assertRaises(ValueError):RelaxationModel(7)
        with self.assertRaises(ValueError):RelaxationModel().initial(amplitude=1)
        with self.assertRaises(ValueError):RelaxationModel().evolve(end=-.1)

if __name__=='__main__':unittest.main()
