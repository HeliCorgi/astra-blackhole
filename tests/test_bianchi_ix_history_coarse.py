from __future__ import annotations
import json,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'research'))
from bianchi_ix_history_coarse import PARTITIONS,audit
class CoarseHistoryTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.out=audit(json.loads((ROOT/'research/bianchi_ix_history_coarse_input.json').read_text()))
 def row(self,n):return next(x for x in self.out['partitions'] if x['name']==n)
 def test_all_five_partitions(self):self.assertEqual(len(PARTITIONS),5)
 def test_full_coarse_is_diagonal(self):self.assertEqual(self.row('ignore_second')['max_offdiagonal_abs'],0);self.assertAlmostEqual(self.row('ignore_second')['trace_diagonal'],1,places=12)
 def test_fine_parent_value(self):self.assertAlmostEqual(self.row('fine_3')['max_offdiagonal_abs'],0.019418109072940764,places=14)
 def test_A_vs_Bs_suppresses_interference(self):self.assertLess(self.row('A_vs_Bs')['max_offdiagonal_abs'],self.row('fine_3')['max_offdiagonal_abs']/3)
 def test_other_binary_splits_above_control(self):
  for n in ('Bplus_vs_rest','Bminus_vs_rest'):self.assertGreater(self.row(n)['strict_signal_to_propagated_control_bound'],5)
if __name__=='__main__':unittest.main()
