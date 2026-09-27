from pathlib import Path
import sys
import unittest
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.steps.step_003_closure_delays_final import v_test_circular


class DirectedVTests(unittest.TestCase):
    def test_opposite_concentration_is_not_evidence_for_reference_direction(self):
        angles = np.linspace(-.1, .1, 50)
        forward, p_forward = v_test_circular(angles)
        opposite, p_opposite = v_test_circular(angles+np.pi)
        self.assertGreater(forward, 9)
        self.assertLess(opposite, -9)
        self.assertLess(p_forward, 1e-15)
        self.assertGreater(p_opposite, 1-1e-15)

    def test_weighted_projection_and_rotation(self):
        angles = np.linspace(-.3, 1.8, 30)
        weights = np.linspace(1, 2, 30)
        expected = np.sqrt(2)*np.sum(weights*np.cos(angles))/np.sqrt(np.sum(weights**2))
        v, p = v_test_circular(angles, weights=weights)
        np.testing.assert_allclose([v, p], [expected, stats.norm.sf(expected)])
        np.testing.assert_allclose(v_test_circular(angles+2., mu0=2., weights=weights), [v,p])


if __name__ == '__main__':
    unittest.main()
