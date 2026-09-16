"""Angular-state checks. Exact algebra is not a proof of physical validity."""
from pathlib import Path
import sys, unittest
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'research'))
from angular_model import AngularBasis,radial_weights,general_jets,angular_identities
from local_jets import CurvatureJets
from run_angular_audit import exact_checks


class TestAngularStates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.centers=np.geomspace(.05,8.,10)
        cls.jets=CurvatureJets()

    def bases(self,ell=4,mass=1.,alpha=.8):
        return [AngularBasis(self.centers,ell,sign*alpha,mass=mass) for sign in (1,-1)]

    def test_exact_angular_cancellation(self):
        for ell in (4,6,8):
            ids=angular_identities(ell)
            for key,value in ids.items():
                i,j=map(int,key.split(','))
                if 2*(i+j)<ell:self.assertEqual(sp.Rational(value),0)
            self.assertNotEqual(sp.Rational(ids[f'{ell//2},0']),0)

    def test_same_initial_stress_and_delayed_curvature_response(self):
        for mass in (1.,0.):
            for ell in (4,6,8):
                bp,bm=self.bases(ell,mass);w=radial_weights(bp)
                y=np.array([0.,0.,.6,-1.])
                cp,cm=[b.tensors(y,w,4) for b in (bp,bm)]
                for i,j in cp:
                    if 2*(i+j)<ell:self.assertAlmostEqual(cp[i,j],cm[i,j],places=12)
                jp,jm=[general_jets(self.jets,b,w,y) for b in (bp,bm)]
                first=ell//2-1
                np.testing.assert_allclose(jp[:first],jm[:first],rtol=1e-12,atol=1e-11)
                self.assertGreater(abs(jp[first]-jm[first]),1e-5)
                self.assertAlmostEqual(jp[0],bp.observe(y,w)['K'],places=11)

    def test_exact_massless_jet_delays(self):
        result=exact_checks(self.jets)
        for row in result:
            gap=[sp.Rational(v) for v in row['exact_massless_K_jet_gaps']]
            first=row['first_different_order']
            self.assertTrue(all(v==0 for v in gap[:first]))
            self.assertNotEqual(gap[first],0)

    def test_pressure_response_is_shear_times_omitted_moment(self):
        for mass in (1.,0.):
            bp,bm=self.bases(mass=mass);w=radial_weights(bp);y=np.array([0.,0.,.6,-1.])
            cp,cm=[b.tensors(y,w,2) for b in (bp,bm)]
            op,om=[b.observe(y,w) for b in (bp,bm)]
            expected=(y[2]-y[3])*(cp[2,0]-cm[2,0])
            self.assertAlmostEqual(op['dpr']-om['dpr'],expected,places=13)
            self.assertAlmostEqual(op['dpt']-om['dpt'],-expected/2,places=13)

    def test_massless_angular_difference_survives(self):
        bp,bm=self.bases(mass=0.);w=radial_weights(bp);y=np.array([0.,0.,.6,-1.])
        k=[b.observe(b.evolve(w,0.,.2,y).sol(.2),w)['K'] for b in (bp,bm)]
        self.assertGreater(abs(k[0]-k[1]),.01)

    def test_imposed_isotropic_scaling_control(self):
        bp,bm=self.bases(mass=0.);w=radial_weights(bp)
        y=np.array([np.log(.8),np.log(.8),-.7,-.7])
        op,om=[b.observe(y,w) for b in (bp,bm)]
        for key in ('rho','pr','pt','dpr','dpt'):
            self.assertAlmostEqual(op[key],om[key],places=12)

    def test_invalid_angular_inputs(self):
        for ell in (2,5,4.5):
            with self.assertRaises(ValueError):AngularBasis(self.centers,ell=ell)
        for alpha in (1.,-1.,float('nan')):
            with self.assertRaises(ValueError):AngularBasis(self.centers,amplitude=alpha)
        with self.assertRaises(ValueError):AngularBasis(self.centers,ell=8,angular_order=8)

    def test_zero_amplitude_is_identical(self):
        bp,bm=self.bases(alpha=0.);w=radial_weights(bp)
        y=np.array([.1,-.2,.4,-1.5])
        np.testing.assert_array_equal(bp.rhs(0.,y,w),bm.rhs(0.,y,w))

if __name__=='__main__':unittest.main()
