from __future__ import annotations
import sys,unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))
from quantum_bianchi_ix_model import QuantumBianchiIX,PacketSpec
from quantum_bianchi_ix_long_model import wall_masks,fidelity
from bianchi_ix_histories import project_sector,linear_regrid,propagate_linear,decoherence_matrix,decoherence_measures


class HistoryTests(unittest.TestCase):
    def test_projectors_resolve_identity_and_are_orthogonal(self):
        g=QuantumBianchiIX(bp_bounds=(-3,3),bm_bounds=(-3,3),nx=30,ny=30)
        v=g.initial(PacketSpec())
        parts=[project_sector(g,v,3.0,k) for k in range(3)]
        self.assertLess(np.linalg.norm(sum(parts)-v),1e-14)
        for i in range(3):
            for j in range(i):
                self.assertLess(abs(np.vdot(parts[i],parts[j])),1e-14)

    def test_linear_regrid_is_linear_and_not_branch_normalized(self):
        a=QuantumBianchiIX(bp_bounds=(-3,3),bm_bounds=(-3,3),nx=30,ny=30)
        b=QuantumBianchiIX(bp_bounds=(-5,5),bm_bounds=(-5,5),nx=50,ny=50)
        v=a.initial(PacketSpec())
        x=project_sector(a,v,2.5,0); y=project_sector(a,v,2.5,1)
        lhs=linear_regrid(2*x-0.3j*y,a,b)
        rhs=2*linear_regrid(x,a,b)-0.3j*linear_regrid(y,a,b)
        self.assertLess(np.linalg.norm(lhs-rhs),1e-12)
        # A branch is intentionally not forced back to unit norm.
        self.assertLess(np.linalg.norm(linear_regrid(x,a,b)),1.0)

    def test_small_two_time_history_closure(self):
        g=QuantumBianchiIX(bp_bounds=(-3,3),bm_bounds=(-3,3),nx=30,ny=30)
        v=g.initial(PacketSpec())
        v=propagate_linear(g,v,2.0,2.1,.05,20)
        first=[project_sector(g,v,2.1,j) for j in range(3)]
        branches={}
        full=propagate_linear(g,v,2.1,2.2,.05,20)
        for j,x in enumerate(first):
            y=propagate_linear(g,x,2.1,2.2,.05,20)
            for k in range(3):
                branches[(j,k)]=project_sector(g,y,2.2,k)
        recombined=sum(branches.values())
        self.assertGreater(fidelity(recombined,full),.999999)
        keys,D=decoherence_matrix(branches)
        self.assertLess(np.max(np.abs(D-D.conj().T)),1e-12)
        # Different final sectors are exactly orthogonal projectors.
        for a,(j,k) in enumerate(keys):
            for b,(jj,kk) in enumerate(keys):
                if k!=kk:
                    self.assertLess(abs(D[a,b]),1e-12)

    def test_decoherence_matrix_positive_semidefinite(self):
        rng=np.random.default_rng(7)
        branches={i:rng.normal(size=20)+1j*rng.normal(size=20) for i in range(5)}
        _,D=decoherence_matrix(branches)
        m=decoherence_measures(D)
        self.assertGreater(m["minimum_eigenvalue_hermitian_part"],-1e-10)
        self.assertLess(m["hermiticity_max_abs"],1e-12)

    def test_no_zero_norm_failure(self):
        g=QuantumBianchiIX(bp_bounds=(-3,3),bm_bounds=(-3,3),nx=30,ny=30)
        z=np.zeros(g.nx*g.ny,complex)
        out=propagate_linear(g,z,2.0,2.2,.1,20)
        self.assertEqual(np.linalg.norm(out),0.0)


if __name__=="__main__":
    unittest.main()
