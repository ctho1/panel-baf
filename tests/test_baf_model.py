import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from baf_common import BGRID  # noqa: E402
from baf_model import robust_rate, snp_loglik, viterbi_segment  # noqa: E402


def simulate(b, n_snps=20000, mean_depth=0.15, seed=1):
    """Low-depth SNPs (HWE genotypes) whose heterozygous minor-allele fraction is b."""
    rng = np.random.default_rng(seed)
    p = rng.uniform(0.05, 0.5, n_snps)
    g = rng.binomial(2, p)                                # alt-allele count of the genotype
    n = rng.poisson(mean_depth * 20, n_snps)              # enriched for covered SNPs
    f = np.where(g == 0, 0.005, np.where(g == 2, 0.995, np.where(rng.random(n_snps) < 0.5, b, 1 - b)))
    k = rng.binomial(n, f)
    keep = n >= 2
    return k[keep], n[keep], p[keep]


class SnpLikelihoodTest(unittest.TestCase):
    def mle(self, b):
        k, n, p = simulate(b)
        ll = snp_loglik(k, n, p, h=1.0, w_ref=1.0).sum(axis=0)
        return BGRID[int(np.argmax(ll))]

    def test_balanced_region(self):
        self.assertGreaterEqual(self.mle(0.5), 0.46)

    def test_single_copy_loss_at_60pct_purity(self):
        self.assertAlmostEqual(self.mle(0.29), 0.29, delta=0.04)

    def test_high_purity_loh(self):
        self.assertLessEqual(self.mle(0.05), 0.10)

    def test_depth_one_snps_are_uninformative(self):
        # without reference bias a single read cannot distinguish b from 1-b -> flat likelihood
        ll = snp_loglik(np.array([0, 1]), np.array([1, 1]), np.array([0.3, 0.3]), w_ref=1.0)
        self.assertTrue(np.allclose(ll - ll[:, :1], 0))


class RobustRateTest(unittest.TestCase):
    def test_hotspots_ignored(self):
        rng = np.random.default_rng(0)
        c = rng.poisson(5.0, 1000).astype(float)
        c[:20] = 500                                       # capture hotspots
        lam, used = robust_rate(c, np.ones_like(c))
        self.assertAlmostEqual(lam, 5.0, delta=0.3)
        self.assertLessEqual(used, 980)


class ViterbiTest(unittest.TestCase):
    def test_detects_step(self):
        rng = np.random.default_rng(0)
        lr = np.r_[rng.normal(0, 0.15, 60), rng.normal(-0.5, 0.15, 60)]
        bll = np.zeros((len(lr), len(BGRID)))
        mu = np.round(np.arange(-1, 1.01, 0.05), 2)
        path = viterbi_segment(lr, bll, 0.15, mu, np.arange(0, len(BGRID), 2), -15.0)
        breaks = np.where(np.diff(path) != 0)[0] + 1
        self.assertEqual(len(breaks), 1)
        self.assertLessEqual(abs(breaks[0] - 60), 2)


if __name__ == "__main__":
    unittest.main()
