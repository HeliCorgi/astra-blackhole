from __future__ import annotations
import json,sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'research'))
from quantum_bianchi_ix_model import QuantumBianchiIX,PacketSpec
from quantum_bianchi_ix_long_model import positive_wall_terms,regrid_state,fidelity,branch_interference
from bianchi_ix_model import evolve,kasner_from_state,bkl_step

class LongAuditTests(unittest.TestCase):
    def test_wall_sectors_partition(self):
        g=QuantumBianchiIX(bp_bounds=(-2,2),bm_bounds=(-2,2),nx=30,ny=30)
        lab=np.argmax(positive_wall_terms(10,g.BP,g.BM),axis=0)
        self.assertEqual(lab.shape,(30,30));self.assertTrue(np.all((lab>=0)&(lab<3)))
    def test_regrid_roundtrip_gaussian(self):
        a=QuantumBianchiIX(bp_bounds=(-4,4),bm_bounds=(-4,4),nx=48,ny=48)
        b=QuantumBianchiIX(bp_bounds=(-6,6),bm_bounds=(-6,6),nx=72,ny=72)
        v=a.initial(PacketSpec());w,_=regrid_state(v,a,b);back,_=regrid_state(w,b,a)
        self.assertGreater(fidelity(v,back),.9999)
    def test_short_branch_recombination(self):
        g=QuantumBianchiIX(bp_bounds=(-4,4),bm_bounds=(-4,4),nx=40,ny=40,potential_cap=30)
        v=g.initial(PacketSpec());r=branch_interference(g,v,2,2.1,.1,20,None)
        self.assertGreater(r['recombination_fidelity'],.999)
        self.assertAlmostEqual(sum(r['split_weights']),1,places=10)
    def test_third_classical_bkl_transition(self):
        ss=np.array([15.,36.]);sol=evolve(np.array([0.,0.,1.,.31]),2,36,t_eval=ss,rtol=2e-9,atol=2e-11,max_step=.02)
        u0=kasner_from_state(15,sol.y[:,0])['u'];u1=kasner_from_state(36,sol.y[:,1])['u']
        self.assertLess(abs(bkl_step(u0)-u1)/u1,2e-6)
    def test_published_baseline_scope(self):
        d=json.loads((ROOT/'research/quantum_bianchi_ix_long_results/summary.json').read_text())
        self.assertEqual([x['central_wall'] for x in d['central_bounces']],['B-','B+'])
        self.assertGreater(d['regrid_cumulative_roundtrip_fidelity_product'],.9997)
        self.assertLess(d['numerical_controls']['max_saved_edge_mass_w1'],3e-4)
        self.assertLess(d['numerical_controls']['max_saved_capped_region_mass'],1e-5)
        self.assertGreater(d['interference']['second_wall_split_9p5_to_11']['interference_density_L1_half'],.08)
        self.assertGreater(d['interference']['third_wall_split_30_to_34']['interference_density_L1_half'],.12)
if __name__=='__main__':unittest.main()
