import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from figure_genome import _repel, expected_levels  # noqa: E402


class FigureHelpersTest(unittest.TestCase):
    def test_expected_levels_at_full_purity(self):
        lev = {c: (lr, b) for c, _, lr, _, b in expected_levels(1.0)}
        self.assertAlmostEqual(lev[1][0], -1.0)              # 1 copy: log2(1/2)
        self.assertAlmostEqual(lev[3][0], np.log2(1.5))
        self.assertAlmostEqual(lev[4][0], 1.0)
        self.assertAlmostEqual(lev[1][1], 0.0)               # LOH: minor allele absent
        self.assertAlmostEqual(lev[3][1], 1 / 3)             # 2:1
        self.assertAlmostEqual(lev[4][1], 0.25)              # 3:1

    def test_expected_levels_at_zero_purity_are_diploid(self):
        for _, _, lr, _, b in expected_levels(0.0):
            self.assertAlmostEqual(lr, 0.0)
            self.assertAlmostEqual(b, 0.5)

    def test_repel_enforces_gap(self):
        out = _repel([10.0, 10.5, 30.0, 11.0], min_gap=2.0)
        s = np.sort(out)
        self.assertTrue(np.all(np.diff(s) >= 2.0 - 1e-9))
        self.assertAlmostEqual(out[2], 30.0)                 # isolated label stays put


if __name__ == "__main__":
    unittest.main()
