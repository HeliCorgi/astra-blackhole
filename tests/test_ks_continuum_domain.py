import importlib.util
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("ks_continuum_domain_audit", ROOT/"research"/"ks_continuum_domain_audit.py")
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

class TestKSContinuumDomain(unittest.TestCase):
    def test_threshold_domain_boundary(self):
        self.assertTrue(mod.threshold_domain_is_finite(.5))
        self.assertTrue(mod.threshold_domain_is_finite(1.0))
        self.assertFalse(mod.threshold_domain_is_finite(1.5))
        self.assertFalse(mod.threshold_domain_is_finite(2.0))

    def test_reference_threshold_coefficient_stable_and_nonzero(self):
        a=mod.threshold_coefficient(h=.2,left=-4,right=8,dx=.01)
        b=mod.threshold_coefficient(h=.2,left=-6,right=10,dx=.005)
        self.assertGreater(abs(a),1e-6)
        self.assertLess(abs(a-b)/abs(b),1e-9)

    def test_finite_box_velocity_sylvester_and_bound(self):
        row=mod.finite_box_velocity_check(right=8,dx=.08)
        self.assertLess(row["sylvester_relative_residual"],1e-10)
        self.assertGreaterEqual(row["min_velocity_eigenvalue"],-1.000001)
        self.assertLessEqual(row["max_velocity_eigenvalue"],1.000001)
        self.assertGreater(row["min_discrete_H2_minus_P2_eigenvalue"],-1e-9)

if __name__=="__main__":
    unittest.main()
