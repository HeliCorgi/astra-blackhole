"""Independent numerical/algebraic checks; not proof of physical equivalence."""
import sys,unittest
from pathlib import Path
import numpy as np
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from chiba_model import Data,Diamond,green,currents,pullback,same_slice_flux,exact_slice

class TestChibaReproduction(unittest.TestCase):
    def test_initial_data(self):
        d=Data();u=np.linspace(-.5,.2,21)
        np.testing.assert_allclose(green(d,u,np.zeros_like(u))[0],d.boundary(u),rtol=0,atol=0)
        self.assertAlmostEqual(abs(d.boundary(d.ub)),np.exp(-8),places=16)

    def test_initial_derivative(self):
        d=Data();u=np.array([-.4,-.3,-.2]);eps=1e-6
        numerical=(d.boundary(u+eps)-d.boundary(u-eps))/(2*eps)
        np.testing.assert_allclose(d.derivative(u),numerical,rtol=2e-9,atol=2e-8)

    def test_green_derivatives(self):
        d=Data();u=np.array([-.3,-.1,.2]);v=np.array([.1,.3,.2]);eps=1e-5
        p,du,dv=green(d,u,v)
        fd_u=(green(d,u+eps,v)[0]-green(d,u-eps,v)[0])/(2*eps)
        fd_v=(green(d,u,v+eps)[0]-green(d,u,v-eps)[0])/(2*eps)
        np.testing.assert_allclose(du,fd_u,rtol=3e-7,atol=1e-7)
        np.testing.assert_allclose(dv,fd_v,rtol=3e-7,atol=1e-7)

    def test_green_mixed_wdw_residual(self):
        d=Data();u=np.array([-.3,-.1,.2]);v=np.array([.1,.3,.2]);eps=1e-5
        p,_,_=green(d,u,v)
        derivative=(green(d,u+eps,v)[2]-green(d,u-eps,v)[2])/(2*eps)
        np.testing.assert_allclose(derivative,-p/d.kappa**2,rtol=3e-7,atol=1e-6)

    def test_exact_plane_wave_stencil(self):
        for delta in (.01,.001):
            z=1j*delta/.1
            self.assertAlmostEqual(abs((2-(1+z))/(1-z)-1),0,places=15)

    def test_diamond_cell_residual(self):
        m=Diamond(Data(),step=.002,umax=.1,vmax=.3);f=m.envelope;z=1j*m.step/.1
        residual=(1-z)*f[1:,1:]-f[1:,:-1]-f[:-1,1:]+(1+z)*f[:-1,:-1]
        self.assertLess(np.max(abs(residual)),2e-14)
        np.testing.assert_array_equal(f[:,0],0)

    def test_field_against_green(self):
        d=Data(ub=-.65);m=Diamond(d,step=.0005,umax=.1,vmax=.4)
        u=np.array([-.3,-.2,-.1]);v=np.array([.05,.2,.3]);exact=green(d,u,v)[0]
        self.assertLess(np.max(abs(m.evaluate(u,v)[0]-exact)),2e-5)

    def test_refinement(self):
        d=Data(ub=-.65);u=np.array([-.3,-.2,-.1]);v=np.array([.05,.2,.3]);target=green(d,u,v)[0]
        errors=[]
        for h in (.002,.001):
            m=Diamond(d,h,umax=.1,vmax=.4);errors.append(np.max(abs(m.evaluate(u,v)[0]-target)))
        self.assertLess(errors[1],errors[0]/2)

    def test_signed_current_not_density(self):
        omega=3.;p=np.array([1+0j]);q,_=currents(p,-.5j*omega*p,-.5j*omega*p)
        qneg,_=currents(p,.5j*omega*p,.5j*omega*p)
        self.assertAlmostEqual(q[0],2*omega);self.assertAlmostEqual(qneg[0],-2*omega)

    def test_charge_conservation(self):
        d=Data(ub=-.65);expected=2*np.sqrt(np.pi)*d.sigma/d.kappa
        for t in (0,.1,.3):
            self.assertLess(abs(exact_slice(d,t,samples=1201)['norm']/expected-1),2e-7)

    def test_coordinate_inverse(self):
        u=np.array([-.4,-.2,-.01]);v=np.array([.1,.2,.4]);x,tau=pullback(u,v)
        np.testing.assert_allclose(-np.exp(-x-tau)/4,u,rtol=1e-14)
        np.testing.assert_allclose(np.exp(-x+tau)/4,v,rtol=1e-14)

    def test_same_segment_flux_and_moment(self):
        for t in (0,.1,.3):
            d=same_slice_flux(Data(),t,order=256)
            self.assertLess(d['relative_difference'],1e-9)
            self.assertLess(d['mean_X_absolute_difference'],1e-10)

    def test_plateau_is_not_L2_gaussian(self):
        d=Data();x=40.;u=-np.exp(-x)/4;v=-u;p,_,_=green(d,u,v,384)
        self.assertGreater(abs(d.boundary(0)),0)
        self.assertLess(abs(p/d.boundary(0)-1),1e-7)

    def test_map_residual_algebra(self):
        x,t=sp.symbols('x t',real=True);u=-sp.exp(-x-t)/4;v=sp.exp(-x+t)/4
        f=sp.Function('f');U,V=sp.symbols('U V');phi=f(U,V)
        Dx=lambda q:-U*sp.diff(q,U)-V*sp.diff(q,V)
        Dt=lambda q:-U*sp.diff(q,U)+V*sp.diff(q,V)
        self.assertEqual(sp.simplify(Dt(Dt(phi))-Dx(Dx(phi))+4*U*V*sp.diff(phi,U,V)),0)
        self.assertEqual(sp.simplify(-4*u*v-sp.exp(-2*x)/4),0)

    def test_invalid_inputs(self):
        for k in (0,-1,float('nan')):
            with self.assertRaises(ValueError):Data(kappa=k)
        with self.assertRaises(ValueError):green(Data(),-.6,.1)
        with self.assertRaises(ValueError):pullback(.1,.1)
        with self.assertRaises(ValueError):Diamond(Data(),step=-.1)

if __name__=='__main__':unittest.main()
