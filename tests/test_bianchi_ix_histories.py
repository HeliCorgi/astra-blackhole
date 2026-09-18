from __future__ import annotations
import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"research"))
from quantum_bianchi_ix_model import QuantumBianchiIX,PacketSpec
from bianchi_ix_histories import (project_sector,linear_regrid,propagate_linear,propagate_block,
                                  decoherence_matrix,decoherence_measures)

class HistoryTests(unittest.TestCase):
    def test_projectors_resolve_identity_and_are_orthogonal(self):
        g=QuantumBianchiIX(bp_bounds=(-3,3),bm_bounds=(-3,3),nx=30,ny=30);v=g.initial(PacketSpec())
        parts=[project_sector(g,v,3.0,k) for k in range(3)]
        self.assertLess(np.linalg.norm(sum(parts)-v),1e-14)
        for i in range(3):
            for j in range(i): self.assertLess(abs(np.vdot(parts[i],parts[j])),1e-14)

    def test_linear_regrid_is_linear_and_not_branch_normalized(self):
        a=QuantumBianchiIX(bp_bounds=(-3,3),bm_bounds=(-3,3),nx=30,ny=30)
        b=QuantumBianchiIX(bp_bounds=(-5,5),bm_bounds=(-5,5),nx=50,ny=50)
        v=a.initial(PacketSpec());x=project_sector(a,v,2.5,0);y=project_sector(a,v,2.5,1)
        self.assertLess(np.linalg.norm(linear_regrid(2*x-.3j*y,a,b)-2*linear_regrid(x,a,b)+.3j*linear_regrid(y,a,b)),1e-12)
        self.assertLess(np.linalg.norm(linear_regrid(x,a,b)),1.0)

    def test_common_block_propagator_preserves_recombination(self):
        g=QuantumBianchiIX(bp_bounds=(-3,3),bm_bounds=(-3,3),nx=30,ny=30)
        v=propagate_linear(g,g.initial(PacketSpec()),2,2.1,.05,20)
        V=np.column_stack([project_sector(g,v,2.1,j) for j in range(3)])
        out,st=propagate_block(g,V,2.1,2.2,.05,30)
        # One common projected unitary preserves the branch Gram matrix.
        self.assertLess(st["max_gram_step_drift"],1e-11)
        branches={(j,k):project_sector(g,out[:,j],2.2,k) for j in range(3) for k in range(3)}
        self.assertLess(np.linalg.norm(sum(branches.values())-out.sum(axis=1))/np.linalg.norm(out.sum(axis=1)),1e-12)
        keys,D=decoherence_matrix(branches);self.assertLess(np.max(np.abs(D-D.conj().T)),1e-12)
        for a,(j,k) in enumerate(keys):
            for b,(jj,kk) in enumerate(keys):
                if k!=kk:self.assertLess(abs(D[a,b]),1e-12)

    def test_decoherence_matrix_positive_semidefinite(self):
        rng=np.random.default_rng(7);branches={i:rng.normal(size=20)+1j*rng.normal(size=20) for i in range(5)}
        _,D=decoherence_matrix(branches);m=decoherence_measures(D)
        self.assertGreater(m["minimum_eigenvalue_hermitian_part"],-1e-10);self.assertLess(m["hermiticity_max_abs"],1e-12)

    def test_zero_block(self):
        g=QuantumBianchiIX(bp_bounds=(-3,3),bm_bounds=(-3,3),nx=30,ny=30)
        z=np.zeros((g.nx*g.ny,3),complex);out,_=propagate_block(g,z,2,2.2,.1,24)
        self.assertEqual(np.linalg.norm(out),0.0)
if __name__=="__main__":unittest.main()
