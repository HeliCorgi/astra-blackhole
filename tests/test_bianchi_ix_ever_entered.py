from __future__ import annotations
import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"research"))
from quantum_bianchi_ix_model import QuantumBianchiIX,PacketSpec
from bianchi_ix_histories import propagate_linear
from bianchi_ix_ever_entered import signed_B_distance,absorber_profile,binary_sector,pair_step,CAPSpec

class EverEnteredTests(unittest.TestCase):
 def test_binary_boundary_is_time_independent_geometry(self):
  bp=np.array([-1.,0.,1.]);bm=np.array([.2,-.2,.2])
  d=signed_B_distance(bp,bm);self.assertEqual(d.shape,(3,));self.assertLess(d[0],0);self.assertGreater(d[2],0)
 def test_smooth_absorber_bounds_and_symmetry(self):
  g=QuantumBianchiIX(bp_bounds=(-2,2),bm_bounds=(-2,2),nx=30,ny=30);f=absorber_profile(g,.3)
  self.assertGreaterEqual(f.min(),0);self.assertLessEqual(f.max(),1)
  self.assertLess(np.max(abs(f-f[:,::-1])),1e-12)
 def test_hard_binary_partitions_grid(self):
  g=QuantumBianchiIX(bp_bounds=(-2,2),bm_bounds=(-2,2),nx=30,ny=30);b=binary_sector(g)
  self.assertEqual(b.dtype,bool);self.assertTrue(np.any(b));self.assertTrue(np.any(~b))
 def test_v0_zero_negative_control(self):
  g=QuantumBianchiIX(bp_bounds=(-3,3),bm_bounds=(-3,3),nx=30,ny=30)
  v=g.initial(PacketSpec());full,no,st=pair_step(g,v,v,2.05,.1,CAPSpec(0,.3,30))
  self.assertLess(np.linalg.norm(full-no),1e-12);self.assertLess(st["gram_step_drift"],1e-11)
 def test_positive_cap_is_contractive_on_no_branch(self):
  g=QuantumBianchiIX(bp_bounds=(-3,3),bm_bounds=(-3,3),nx=30,ny=30)
  v=propagate_linear(g,g.initial(PacketSpec()),2,2.1,.05,20);full,no,_=pair_step(g,v,v,2.15,.1,CAPSpec(.1,.3,30))
  self.assertLessEqual(np.linalg.norm(no),np.linalg.norm(full)+1e-12)
if __name__=="__main__":unittest.main()
