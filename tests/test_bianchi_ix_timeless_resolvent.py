from __future__ import annotations
import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))
from bianchi_ix_timeless_smatrix import WDWGridSpec,build_constraint,near_zero_modes
from bianchi_ix_timeless_resolvent import resolvent_action,dyson_relative,t_action,shell_fraction

class TimelessResolventTests(unittest.TestCase):
    def small(self):
        return build_constraint(WDWGridSpec(
            ns=5,nx=6,ny=6,s_bounds=(2,5),
            bp_bounds=(-2.5,2.5),bm_bounds=(-2.5,2.5),
            potential_abs_cap=5))

    def test_eta_positive(self):
        g=self.small()
        with self.assertRaises(ValueError):
            resolvent_action(g["C"],np.ones(g["C"].shape[0]),0.0,0.0)

    def test_dyson_identity(self):
        g=self.small()
        n=g["C"].shape[0]
        b=np.arange(1,n+1,dtype=float)
        b=b/np.linalg.norm(b)
        vd=.03*g["F_B"].ravel()
        r=dyson_relative(g["C"],b,0.0,.2,vd)
        self.assertLess(r,2e-11)

    def test_zero_potential_has_zero_t(self):
        g=self.small()
        n=g["C"].shape[0]
        v=np.ones(n,complex)/np.sqrt(n)
        t,x=t_action(g["C"],v,0.0,.2,np.zeros(n))
        self.assertLess(np.linalg.norm(t),1e-15)
        self.assertLess(np.linalg.norm(x),1e-15)

    def test_shell_fraction_bounds(self):
        g=self.small()
        vals,U=near_zero_modes(g["C"],3)
        v=U[:,0]+.3j*U[:,1]
        f=shell_fraction(U,v)
        self.assertGreaterEqual(f,0.0)
        self.assertLessEqual(f,1.0+1e-12)

if __name__=="__main__":
    unittest.main()
