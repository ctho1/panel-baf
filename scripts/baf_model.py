"""Statistical core: low-depth SNP likelihood over the het minor-allele fraction b, robust Poisson rate for
off-target windows, and the joint (logR, b) Viterbi segmentation."""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import betaln, gammaln, logsumexp

from baf_common import BGRID

EPS = 0.005      # per-base error towards the specific other allele
RHO = 0.002      # beta-binomial overdispersion of het allele counts
DELTA = 0.01     # outlier floor (paralogous loci, artefacts)


def snp_loglik(k, n, p, h=1.0, w_ref=0.96, bgrid=BGRID):
    """log L(b) for each SNP (rows) over the b grid (columns).

        L(b) = (1-d) [ P(AA) Bin(k|n,e) + P(BB) Bin(k|n,1-e) + P(AB) 1/2 (BB(k|n,b') + BB(k|n,1-b')) ] + d/(n+1)

    k: alt reads, n: ref+alt reads, p: population alt-allele frequency, h: het-rate scale (HWE * h),
    w_ref: alt/ref sampling ratio (reference bias), b' = b with reference bias applied.
    """
    k = np.asarray(k, float)[:, None]
    n = np.asarray(n, float)[:, None]
    p = np.clip(np.asarray(p, float), 1e-3, 1 - 1e-3)[:, None]
    lchoose = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)

    def binom(f):
        return lchoose + k * np.log(f) + (n - k) * np.log1p(-f)

    def bbin(f):
        f = np.clip(f, EPS, 1 - EPS)
        a, b = f * (1 - RHO) / RHO, (1 - f) * (1 - RHO) / RHO
        return lchoose + betaln(k + a, n - k + b) - betaln(a, b)

    def biased(f):
        return f * w_ref / (f * w_ref + (1 - f))

    # HWE with het rate scaled by h; the het excess/deficit is taken equally from both homozygotes
    dh = (h - 1) * 2 * p * (1 - p)
    phr = np.clip((1 - p) ** 2 - dh / 2, 1e-6, None)
    pha = np.clip(p ** 2 - dh / 2, 1e-6, None)
    phet = np.clip(2 * p * (1 - p) + dh, 1e-6, None)
    tot = phr + pha + phet
    phr, pha, phet = phr / tot, pha / tot, phet / tot
    bg = bgrid[None, :]
    lhet = np.logaddexp(bbin(biased(bg)), bbin(biased(1 - bg))) - np.log(2)
    comps = np.stack([np.broadcast_to(np.log(phr) + binom(np.full_like(p, EPS)), lhet.shape),
                      np.broadcast_to(np.log(pha) + binom(np.full_like(p, 1 - EPS)), lhet.shape),
                      np.log(phet) + lhet])
    return np.logaddexp(np.log1p(-DELTA) + logsumexp(comps, axis=0), np.log(DELTA) - np.log(n + 1))


def thin(df, min_gap):
    """Keep one SNP per `min_gap` bp (avoid counting the same fragment twice). df sorted by chrom, pos."""
    keep = np.zeros(len(df), bool)
    last_c, last_p = None, -1e18
    for i, (c, p) in enumerate(zip(df.chrom.values, df.pos.values)):
        if c != last_c or p - last_p > min_gap:
            keep[i] = True
            last_c, last_p = c, p
    return df[keep]


def profiled_ll(snps, afcol, h, groups, w_ref):
    """Sum over groups (chromosome arms) of max_b sum_snps log L(b)."""
    ll = snp_loglik(snps.alt_n.values, snps.depth.values, snps[afcol].values, h, w_ref)
    g = pd.factorize(groups)[0]
    agg = np.zeros((g.max() + 1, len(BGRID)))
    np.add.at(agg, g, ll)
    return agg.max(axis=1).sum()


def robust_rate(c, e, n_iter=20, q=0.999):
    """Poisson rate multiplier for windows with counts c and expectations e; windows above the q-quantile
    (allowing 2x overdispersion) are treated as capture hotspots and discarded."""
    c = np.asarray(c, float)
    e = np.asarray(e, float)
    lam = max(np.median(c / e), 0.05)
    keep = np.ones(len(c), bool)
    for _ in range(n_iter):
        keep = c <= stats.poisson.ppf(q, 2 * lam * e) + 1
        new = c[keep].sum() / e[keep].sum()
        if abs(new - lam) < 1e-4 * lam:
            break
        lam = new
    return lam, int(keep.sum())


def viterbi_segment(lr, bll, sigma, mu_grid, b_idx, jump_logp):
    """Joint (logR, b) HMM over bins; Gaussian logR emission, BAF emission = per-bin log L(b);
    uniform jump to any state with log-prob jump_logp. Returns the state path."""
    T = len(lr)
    em_lr = np.zeros((T, len(mu_grid)))
    has = ~np.isnan(lr)
    em_lr[has] = stats.norm.logpdf((lr[has, None] - mu_grid[None, :]) / sigma) - np.log(sigma)
    em = (em_lr[:, :, None] + bll[:, b_idx][:, None, :]).reshape(T, -1)
    stay = np.log1p(-np.exp(jump_logp))
    delta = em[0].copy()
    bp = np.zeros((T, em.shape[1]), np.int32)
    for t in range(1, T):
        best = int(np.argmax(delta))
        st = delta + stay
        jp = delta[best] + jump_logp
        use_jump = jp > st
        bp[t] = np.where(use_jump, best, np.arange(em.shape[1]))
        delta = np.where(use_jump, jp, st) + em[t]
    path = np.zeros(T, np.int64)
    path[-1] = int(np.argmax(delta))
    for t in range(T - 1, 0, -1):
        path[t - 1] = bp[t, path[t]]
    return path
