from __future__ import annotations
import sys,unittest
from pathlib import Path
import numpy as np
from scipy import sparse
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"research"))
from bianchi_ix_timeless_smatrix import WDWGridSpec,build_constraint,near_zero_modes,low_energy_edge_min_packet,smatrix_action,commutator_metrics
class TimelessSMatrixTests(unittest.TestCase):
    def small(self):return build_constraint(WDWGridSpec(ns=6,nx=8,ny=8,s_bounds=(2,5),bp_bounds=(-2.5,2.5),bm_bounds=(-2.5,2.5),potential_abs_cap=5))
    def test_constraint_hermitian(self):
        g=self.small();self.assertLess(sparse.linalg.norm(g["C"]-g["C"].getH()),1e-13)
    def test_absorber_bounds(self):
        g=self.small();self.assertGreaterEqual(g["F_B"].min(),0);self.assertLessEqual(g["F_B"].max(),1)
    def test_near_zero_mode(self):
        g=self.small();vals,U=near_zero_modes(g["C"],3);self.assertEqual(U.shape[1],3);self.assertTrue(np.all(np.diff(np.abs(vals))>=-1e-12))
    def test_edge_min_packet(self):
        g=self.small();vals,p,m=low_energy_edge_min_packet(g,6);self.assertAlmostEqual(np.linalg.norm(p),1,places=12);self.assertLessEqual(m["edge_mass_one_cell"],1)
    def test_v0_zero_identity_and_commutator(self):
        g=self.small();vals,U=near_zero_modes(g["C"],3);v=U[:,0]+.2j*U[:,1];v/=np.linalg.norm(v)
        y=smatrix_action(g["C"],g["F_B"],v,.25,0,g["spec"].hbar);self.assertLess(np.linalg.norm(y-v),2e-8)
        r=commutator_metrics(g["C"],g["F_B"],v,.25,0,g["spec"].hbar);self.assertLess(r["relative"],2e-8)
if __name__=="__main__":unittest.main()
