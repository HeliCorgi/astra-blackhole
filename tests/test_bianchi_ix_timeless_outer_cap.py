from __future__ import annotations
import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))

from bianchi_ix_timeless_smatrix import WDWGridSpec,build_constraint,near_zero_modes
from bianchi_ix_timeless_resolvent import t_action as plain_t_action
from bianchi_ix_timeless_outer_cap import outer_cap_profile,dyson_relative,t_action

class TimelessOuterCAPTests(unittest.TestCase):
    def small(self):
        return build_constraint(WDWGridSpec(
            ns=5,nx=6,ny=6,s_bounds=(2,5),
            bp_bounds=(-2.5,2.5),bm_bounds=(-2.5,2.5),
            potential_abs_cap=5))

    def test_outer_profile(self):
        p=outer_cap_profile((7,8,9),3)
        self.assertEqual(p.shape,(7,8,9))
        self.assertGreaterEqual(float(p.min()),0.0)
        self.assertLessEqual(float(p.max()),1.0)
        self.assertAlmostEqual(float(p[0,4,4]),1.0,places=12)
        self.assertEqual(float(p[3,4,4]),0.0)

    def test_dyson_identity_with_outer_absorber(self):
        g=self.small();n=g["C"].shape[0]
        b=np.arange(1,n+1,dtype=float);b/=np.linalg.norm(b)
        outer=.1*outer_cap_profile(g["shape"],2).ravel()
        vd=.03*g["F_B"].ravel()
        r=dyson_relative(g["C"],b,0.0,.2,outer,vd)
        self.assertLess(r,2e-11)

    def test_zero_outer_matches_plain_resolvent(self):
        g=self.small();vals,U=near_zero_modes(g["C"],3);v=U[:,0]
        vd=.03*g["F_B"].ravel();zero=np.zeros(g["C"].shape[0])
        a,_=t_action(g["C"],v,0.0,.2,zero,vd)
        b,_=plain_t_action(g["C"],v,0.0,.2,vd)
        self.assertLess(np.linalg.norm(a-b),2e-12)

    def test_zero_physical_absorber_has_zero_t(self):
        g=self.small();vals,U=near_zero_modes(g["C"],3);v=U[:,0]
        outer=.1*outer_cap_profile(g["shape"],2).ravel()
        t,x=t_action(g["C"],v,0.0,.2,outer,np.zeros(g["C"].shape[0]))
        self.assertLess(np.linalg.norm(t),1e-15)
        self.assertLess(np.linalg.norm(x),1e-15)

if __name__=="__main__":
    unittest.main()
