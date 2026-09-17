import json,math,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'research'))
from bianchi_ix_branch_info import entropy_bits,effective_branches,total_variation,js_bits,audit

class TestBranchInfo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((ROOT/'research/quantum_bianchi_ix_long_results/summary.json').read_text())
        cls.out=audit(cls.data)
    def test_entropy_limits(self):
        self.assertAlmostEqual(entropy_bits([1,0,0]),0.0,places=14)
        self.assertAlmostEqual(entropy_bits([1/3]*3),math.log2(3),places=14)
        self.assertAlmostEqual(effective_branches([1/3]*3),3.0,places=13)
    def test_distances(self):
        self.assertAlmostEqual(total_variation([1,0,0],[0,1,0]),1.0,places=14)
        self.assertAlmostEqual(js_bits([1,0,0],[1,0,0]),0.0,places=14)
    def test_third_wall_nearly_three_sector(self):
        row=next(x for x in self.out['selected'] if x['s']==32.46)
        self.assertGreater(row['quantum']['maximum_entropy_fraction'],0.99)
        self.assertGreater(row['quantum']['effective_branch_count'],2.98)
        self.assertLess(row['quantum']['dominant_sector_probability'],0.38)
    def test_third_wall_difference_exceeds_numerical_control_scale(self):
        n=self.out['numerical_scale_comparison']
        self.assertGreater(n['third_wall_ratio_to_max_local_control'],50)
        self.assertGreater(n['third_wall_interference_ratio_to_max_local_control'],10)
    def test_second_wall_quantum_and_classical_close(self):
        row=next(x for x in self.out['selected'] if x['s']==10.385)
        self.assertLess(row['quantum_classical_tv'],0.005)
        self.assertLess(row['quantum_classical_js_bits'],5e-5)

if __name__=='__main__':unittest.main()
