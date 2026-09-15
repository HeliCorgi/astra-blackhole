from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'research'))
from predictors import polynomial,neural,power_law
from kinetic_model import MomentumBasis,matched_pair

class TestInterfaces(unittest.TestCase):
    def test_exact_power_law(self):
        t=np.linspace(-.1,0,41);k=(1-t)**(-4)
        p=polynomial(t,k,.1)
        for v in p['predicted_K_by_degree'].values():
            self.assertAlmostEqual(v,.9**(-4),places=11)
        self.assertAlmostEqual(power_law(1.,4.,.1)['predicted_K'],.9**(-4),places=12)
    def test_invalid_samples(self):
        with self.assertRaises(ValueError):polynomial([0.]*8,[1.]*8,.1)
        with self.assertRaises(ValueError):polynomial(range(8),[-1.]*8,.1)
        with self.assertRaises(ValueError):power_law(1.,4.,1.)
        with self.assertRaises(ValueError):neural(np.linspace(-.1,0,41),np.ones(41),1.)

class TestMatchedPhysicalStates(unittest.TestCase):
    def test_same_stress_and_curvature(self):
        b=MomentumBasis(np.geomspace(.05,8.,10))
        for level in (1,2,3):
            p=matched_pair(b,level)
            y1=b.initial_geometry(p['plus']);y2=b.initial_geometry(p['minus'])
            np.testing.assert_allclose(y1,y2,rtol=0,atol=1e-13)
            o1=b.observe(y1,p['plus']);o2=b.observe(y2,p['minus'])
            for key in ['rho','pr','pt','n','K']:
                self.assertAlmostEqual(o1[key],o2[key],places=12)
            self.assertGreater(p['first_omitted_gap'],1e-5)
            self.assertGreater(p['plus'].min(),0)
            self.assertGreater(p['minus'].min(),0)
            self.assertLess(o1['outgoing_expansion'],0)
            self.assertLess(abs(o1['constraint']),1e-12)

if __name__=='__main__':unittest.main()
