import importlib.util
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"))
spec=importlib.util.spec_from_file_location(
    "ks_affine_mass_obstruction_audit",
    ROOT/"research"/"ks_affine_mass_obstruction_audit.py",
)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class TestKSAffineMassObstruction(unittest.TestCase):
    def test_classical_affine(self):
        r=mod.classical_affine_symbolic()
        self.assertEqual(r["poisson_u_E"],"1")
        self.assertEqual(r["poisson_E_mu_plus_mu"],"0")
        self.assertEqual(r["poisson_mu_E_minus_mu"],"0")

    def test_null_bridge(self):
        r=mod.null_canonical_symbolic()
        self.assertIn("p_U*p_W - 4",r["constraint_factor"])
        self.assertEqual(r["mass_pT_pW_over_8"],"p_W*(-U*p_U + W*p_W)/8")

    def test_finite_box_negative_control(self):
        r=mod.finite_box_affine_residual(right=8.,dx=.08)
        self.assertGreater(r["relative_residual_for_[H,M0]=-ihM0"],0.1)
        self.assertLess(r["commutator_antihermiticity_residual"],2e-12)

    def test_no_go_record(self):
        r=mod.affine_weyl_no_go()
        self.assertEqual(r["status"],"PASS")
        self.assertIn("can obey",r["conclusion"])


if __name__=="__main__":
    unittest.main()
