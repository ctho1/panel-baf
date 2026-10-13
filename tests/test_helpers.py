import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from baf_model import lowess  # noqa: E402


class LowessTest(unittest.TestCase):
    def test_linear_data_is_reproduced(self):
        x = np.linspace(0, 1, 200)
        fit = lowess(2 * x + 1, x, frac=0.3)
        self.assertTrue(np.allclose(fit[:, 1], 2 * fit[:, 0] + 1, atol=1e-9))

    def test_robust_to_outliers(self):
        rng = np.random.default_rng(0)
        x = rng.uniform(0, 1, 500)
        y = np.sin(3 * x) + rng.normal(0, 0.05, 500)
        y[:25] += 5                                          # gross outliers
        fit = lowess(y, x, frac=0.3, it=3)
        self.assertLess(np.abs(fit[:, 1] - np.sin(3 * fit[:, 0])).max(), 0.1)


if __name__ == "__main__":
    unittest.main()
