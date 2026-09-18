from __future__ import annotations
import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"research"))
from bianchi_ix_classical_histories import summarize,history_labels

class ClassicalHistoryTests(unittest.TestCase):
    def test_summary_counts(self):
        x=np.array([[0,0,1,2],[2,2,1,0],[1,1,2,0]])
        s=summarize(x);self.assertEqual(sum(r["count"] for r in s["full_27"]),4)
        self.assertEqual(s["first_A_count"],2)
        self.assertAlmostEqual(sum(r["conditional_frequency"] for r in s["A_conditioned_9"]),1.0)
    def test_baseline_small_sample_labels(self):
        x=history_labels(n=256,ds=.02);self.assertEqual(x.shape,(3,256))
        self.assertTrue(np.all((x>=0)&(x<3)))
    def test_step_refinement_small_sample(self):
        a=history_labels(n=256,ds=.02);b=history_labels(n=256,ds=.01)
        self.assertLess(np.mean(np.any(a!=b,axis=0)),.02)
if __name__=="__main__":unittest.main()
