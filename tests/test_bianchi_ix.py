import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'research'))
from bianchi_ix_model import potential,wall_term_and_gradient,evolve,kasner_from_state,bkl_step,first_bkl_divergence

class TestBianchiIX(unittest.TestCase):
    def test_threefold_potential_symmetry(self):
        th=2*np.pi/3;R=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
        b=np.array([.17,-.09]);self.assertAlmostEqual(potential(*b),potential(*(R@b)),places=11)
    def test_gradient_finite_difference(self):
        s,bp,bm=3.1,.13,-.07;W,gp,gm=wall_term_and_gradient(s,bp,bm);h=1e-6
        fp=(wall_term_and_gradient(s,bp+h,bm)[0]-wall_term_and_gradient(s,bp-h,bm)[0])/(2*h)
        fm=(wall_term_and_gradient(s,bp,bm+h)[0]-wall_term_and_gradient(s,bp,bm-h)[0])/(2*h)
        self.assertAlmostEqual(gp,fp,places=7);self.assertAlmostEqual(gm,fm,places=7)
    def test_continuous_bounces_match_bkl(self):
        ts=np.array([2.,6.,15.]);sol=evolve([0,0,1,.31],2,15,t_eval=ts)
        k=[kasner_from_state(t,y) for t,y in zip(ts,sol.y.T)];u=[x['u'] for x in k]
        self.assertLess(abs(u[1]/bkl_step(u[0])-1),3e-4)
        self.assertLess(abs(u[2]/bkl_step(u[1])-1),3e-5)
        for x in k:
            self.assertAlmostEqual(x['sum'],1,places=12)
            self.assertLess(x['wall_ratio'],4e-4)
    def test_bkl_sensitive_but_deterministic(self):
        rows=first_bkl_divergence(1.7112429032939538,1e-12,.1,100)
        self.assertEqual(rows[-1][0],35);self.assertGreater(rows[-1][3],.1)
        self.assertEqual(first_bkl_divergence(1.7112429032939538,0,1e-9,20)[-1][3],0)
    def test_solver_refinement(self):
        t=np.array([2.,6.,15.]);a=evolve([0,0,1,.31],2,15,t_eval=t,max_step=.01)
        b=evolve([0,0,1,.31],2,15,t_eval=t,max_step=.005,rtol=1e-11,atol=1e-13)
        np.testing.assert_allclose(a.y,b.y,rtol=2e-8,atol=2e-10)
if __name__=='__main__':unittest.main()
