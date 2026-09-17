import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'research'))
from quantum_bianchi_ix_model import QuantumBianchiIX,PacketSpec,wall_grid,propagate_quantum,initial_wigner_ensemble
class TestQBIX(unittest.TestCase):
 def test_initial_normalization(self):
  g=QuantumBianchiIX(nx=40,ny=44);v=g.initial(PacketSpec());self.assertAlmostEqual(np.linalg.norm(v),1.,12)
 def test_threefold_wall_symmetry(self):
  b1,b2=.21,-.37;w=wall_grid(2,np.array([[b1]]),np.array([[b2]]))[0,0]
  th=2*np.pi/3;c,s=np.cos(th),np.sin(th);x,y=c*b1-s*b2,s*b1+c*b2
  wr=wall_grid(2,np.array([[x]]),np.array([[y]]))[0,0];self.assertAlmostEqual(w,wr,11)
 def test_short_unitarity(self):
  g=QuantumBianchiIX(nx=48,ny=52,potential_cap=30);r=propagate_quantum(g,PacketSpec(),2,2.1,.05,24,(2,2.1));self.assertLess(abs(r['rows'][-1]['norm']-1),1e-11)
 def test_cap_region_small_initial(self):
  g=QuantumBianchiIX();o=g.observables(g.initial(PacketSpec()),2);self.assertLess(o['capped_region_mass'],2e-6)
 def test_wigner_is_classically_allowed(self):
  p=PacketSpec();Y=initial_wigner_ensemble(p,1024,7);W=wall_grid(2,Y[0],Y[1]);self.assertTrue(np.all(Y[2]**2+Y[3]**2+W>0))
if __name__=='__main__':unittest.main()
