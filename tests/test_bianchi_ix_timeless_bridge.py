from __future__ import annotations
import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))
from quantum_bianchi_ix_model import QuantumBianchiIX,PacketSpec
from bianchi_ix_timeless_bridge import apply_sqrt_A,dotG_action,mismatch_row

class TimelessBridgeTests(unittest.TestCase):
    def test_sqrt_action_on_dense_eigenvector(self):
        g=QuantumBianchiIX(hbar=.2,bp_bounds=(-2,2),bm_bounds=(-2,2),nx=20,ny=20,potential_cap=30)
        n=g.nx*g.ny
        A=np.column_stack([g.apply_A(np.eye(n)[:,j],3.0) for j in range(n)])
        A=(A+A.T.conj())/2
        ev,U=np.linalg.eigh(A)
        idx=n//3
        v=U[:,idx]
        got,_=apply_sqrt_A(g,v,3.0,24)
        want=np.sqrt(max(ev[idx],0))*v
        self.assertLess(np.linalg.norm(got-want)/np.linalg.norm(want),1e-9)

    def test_dotG_nonzero_for_bianchi_ix(self):
        g=QuantumBianchiIX(hbar=.2,bp_bounds=(-3,3),bm_bounds=(-3,3),nx=24,ny=24)
        v=g.initial(PacketSpec())
        x=dotG_action(g,v,3.0,.02,36)
        self.assertGreater(np.linalg.norm(x),1e-5)

    def test_mismatch_is_positive(self):
        g=QuantumBianchiIX(hbar=.2,bp_bounds=(-3,3),bm_bounds=(-3,3),nx=24,ny=24)
        v=g.initial(PacketSpec())
        r=mismatch_row(g,3.0,.02,36,v)
        self.assertGreater(r["relative_naive_constraint_mismatch"],1e-6)

    def test_delta_refinement_stable_small_grid(self):
        g=QuantumBianchiIX(hbar=.2,bp_bounds=(-3,3),bm_bounds=(-3,3),nx=24,ny=24)
        v=g.initial(PacketSpec())
        a=mismatch_row(g,3.0,.02,42,v)["relative_naive_constraint_mismatch"]
        b=mismatch_row(g,3.0,.01,42,v)["relative_naive_constraint_mismatch"]
        self.assertLess(abs(a-b)/max(abs(b),1e-12),.08)

    def test_zero_vector_sqrt(self):
        g=QuantumBianchiIX(hbar=.2,bp_bounds=(-3,3),bm_bounds=(-3,3),nx=24,ny=24)
        z=np.zeros(g.nx*g.ny,complex)
        x,m=apply_sqrt_A(g,z,3.0,20)
        self.assertEqual(np.linalg.norm(x),0)
        self.assertEqual(m["dimension"],0)

if __name__=="__main__":unittest.main()
