import sys, unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'research'))
from trajectory_chaos_model import *

class TestTrajectoryChaos(unittest.TestCase):
    def test_ks_exact_energy(self):
        t=np.linspace(0,100,201);x,p=ks_path(2,-1,t)
        np.testing.assert_allclose(ks_energy(x,p),ks_energy(2,-1),rtol=2e-14,atol=2e-14)
    def test_ks_long_time_ftle_tends_zero(self):
        for pert in [(1e-8,0),(0,1e-8),(1e-8,1e-8)]:
            lam,_=finite_time_lyapunov_ks(2,-1,pert,np.array([1000.,10000.]))
            self.assertLess(abs(lam[-1]),2e-4)
            self.assertLess(abs(lam[-1]),abs(lam[0])+1e-12)
    def test_gauss_exact_average(self):
        est=gauss_ensemble_lyapunov(300000,12345)
        self.assertLess(abs(est/gauss_exact_lyapunov()-1),.005)
    def test_gauss_single_orbit_positive(self):
        self.assertGreater(gauss_finite_lyapunov(np.pi-3,10000),2.0)
    def test_deterministic_repeat(self):
        a=gauss_orbit(np.pi-3,100);b=gauss_orbit(np.pi-3,100)
        self.assertTrue(np.array_equal(a,b))
    def test_nearby_gauss_orbits_separate(self):
        a=gauss_orbit(np.pi-3,40);b=gauss_orbit(np.pi-3+1e-12,40)
        self.assertTrue(np.any(abs(a-b)>1e-3))
    def test_invalid(self):
        with self.assertRaises(ValueError):gauss_map(0)
        with self.assertRaises(ValueError):finite_time_lyapunov_ks(perturb=(0,0))
if __name__=='__main__':unittest.main()
