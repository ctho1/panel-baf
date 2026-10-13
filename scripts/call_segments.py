#!/usr/bin/env python
"""
Segment / arm calls (amplification, gain, loss/LOH, CN-LOH) and tumour-cell-content estimate.

Purity/state fit ("ASCAT-lite"; the dominant state is assumed to be 1+1):
  expected b      = (1-p + p*nB) / (2(1-p) + p*(nA+nB))
  expected logR   = gamma * log2( (2(1-p) + p*(nA+nB)) / 2 )
  score(p, gamma) = sum_seg max_state [ LL_b,seg(b_exp) + LL_logR,seg ] - penalty(state != 1+1)
LL_b is the exact profile log-likelihood of each segment; LL_logR is Gaussian with the segment SE plus a
systematic floor. gamma < 1 allows for compression of the off-target log2 ratio.
Outputs (tmp/<sample>/): segments_called.tsv, arms_called.tsv, purity.json
"""
import json

import numpy as np
import pandas as pd

from baf_common import BGRID, load_arms, parse_args

LOGR_FLOOR = 0.06      # systematic logR uncertainty per segment
STATE_PENALTY = 3.0    # log-lik penalty for any state other than 1+1
GAIN_T, LOSS_T, AMP_T = 0.15, -0.15, 1.0
STATES = [(t - b, b) for t in range(1, 6) for b in range(0, t // 2 + 1)]


def state_ll(ll_curve, logR, se, p, gamma, nA, nB):
    nt = nA + nB
    b_exp = (1 - p + p * nB) / (2 * (1 - p) + p * nt)
    lb = np.interp(b_exp, BGRID, ll_curve)
    ll_r = 0.0
    if not np.isnan(logR):
        lr_exp = gamma * np.log2((2 * (1 - p) + p * nt) / 2)
        ll_r = -0.5 * ((logR - lr_exp) / np.sqrt(se ** 2 + LOGR_FLOOR ** 2)) ** 2
    return lb + ll_r - (0 if (nA, nB) == (1, 1) else STATE_PENALTY)


def _se(r):
    return r.logR_se if not np.isnan(r.logR_se) else 0.1


def fit(table, curves, mask):
    best, surface = (-np.inf, None, None), {}
    rows = [(r, c - c.max()) for (_, r), c in zip(table[mask].iterrows(), curves[mask])]
    for p in np.round(np.arange(0.20, 1.001, 0.01), 2):
        for g in np.round(np.arange(0.5, 1.21, 0.05), 2):
            tot = sum(max(state_ll(c, r.logR, _se(r), p, g, a, b) for a, b in STATES) for r, c in rows)
            surface[(p, g)] = tot
            if tot > best[0]:
                best = (tot, p, g)
    return best, surface


def assign(table, curves, p, g):
    out = []
    for (_, r), c in zip(table.iterrows(), curves):
        cc = c - c.max()
        out.append(max(STATES, key=lambda s: state_ll(cc, r.logR, _se(r), p, g, *s)))
    return out


def call(r, state, cen):
    """Final call from the fitted allelic state, checked against the raw evidence.
    Imbalance p-values ignore LD between SNPs -> strict thresholds."""
    nA, nB = state
    nt = nA + nB
    se = np.sqrt(_se(r) ** 2 + LOGR_FLOOR ** 2)
    z = r.logR / se if not np.isnan(r.logR) else 0.0
    pimb = r.p_imbalance
    if not np.isnan(r.logR) and r.logR >= AMP_T:
        return "Amplification", "high" if pimb < 1e-8 else "medium"
    if r.n_bins < 3:
        return "Not assessed (<3 Mb)", ""
    cs, ce = cen[r.chrom]
    if r.n_bins < 12 and r.start < ce + 6e6 and r.end > cs - 6e6 and pimb > 1e-4:
        return "Pericentromeric, not assessed", ""     # segmental-duplication-rich DNA
    if (nA, nB) == (1, 1):
        if abs(z) > 4 and abs(r.logR) > 0.25:
            return ("Loss (BAF n.s.)" if r.logR < 0 else "Gain (BAF n.s.)"), "low"
        if pimb < 0.01 and r.b_hi <= 0.45:
            return "Possible imbalance", "low"
        return ("Neutral" if r.b_lo >= 0.44 else "Neutral (BAF uninformative)"), ""
    conf = "high" if (pimb < 1e-4 and abs(z) > 3) or pimb < 1e-8 else ("medium" if pimb < 1e-3 or abs(z) > 3 else "low")
    if nt == 1:
        return ("Loss/LOH" if pimb < 1e-3 else "Loss (BAF n.s.)"), conf
    if nt == 2 and nB == 0:
        return "CN-LOH", conf
    if nt >= 3 and nA == nB:
        return "Gain (balanced)", conf
    if nt >= 3:
        return ("Gain (imbalanced)" if pimb < 1e-3 else "Gain (BAF n.s.)"), conf
    return f"State {nA}+{nB}", conf


def main():
    P, _ = parse_args(__doc__)
    cen = load_arms(P.ref)
    segs = pd.read_csv(P.t("segments.tsv"), sep="\t")
    arms = pd.read_csv(P.t("arms.tsv"), sep="\t")
    seg_ll, arm_ll = np.load(P.t("seg_ll.npy")), np.load(P.t("arm_ll.npy"))

    m = ((segs.n_bins >= 5) & ~(segs.logR >= AMP_T)).values     # fit on segments >= 5 Mb, no focal amps
    (score, p_hat, g_hat), surface = fit(segs, seg_ll, m)
    prof = {}
    for (p, g), v in surface.items():
        prof[p] = max(prof.get(p, -np.inf), v)
    ps = np.array(sorted(prof)); pv = np.array([prof[p] for p in ps])
    p_ci = [float(ps[pv >= pv.max() - 2].min()), float(ps[pv >= pv.max() - 2].max())]
    informative = any(s != (1, 1) for s in assign(segs[m], seg_ll[m], p_hat, g_hat))

    for tab, curves, name in [(segs, seg_ll, "segments"), (arms, arm_ll, "arms")]:
        st = assign(tab, curves, p_hat, g_hat)
        cc = [call(r, s, cen) for (_, r), s in zip(tab.iterrows(), st)]
        tab["call"] = [c[0] for c in cc]
        tab["confidence"] = [c[1] for c in cc]
        tab["fit_state"] = [f"{a}+{b}" for a, b in st]
        tab["fit_total_cn"] = [a + b for a, b in st]
        imp = []                       # purity implied by b alone for the fitted state
        for (_, r), (a, b) in zip(tab.iterrows(), st):
            den = r.b * (a + b - 2) + 1 - b
            imp.append(float(np.clip((1 - 2 * r.b) / den, 0, 1)) if a != b and den != 0 else np.nan)
        tab["purity_from_b"] = imp
        tab.to_csv(P.t(f"{name}_called.tsv"), sep="\t", index=False, float_format="%.4g")

    res = dict(purity=p_hat if informative else None, purity_CI_2LL=p_ci if informative else None,
               logR_gamma=g_hat if informative else None, score=score, fit_segments=int(m.sum()),
               state_penalty=STATE_PENALTY,
               note=None if informative else "no imbalanced/aneuploid segment -> purity not estimable")
    json.dump(res, open(P.t("purity.json"), "w"), indent=1)
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
    print(f"purity={res['purity']} (CI {res['purity_CI_2LL']}) gamma={res['logR_gamma']}")
    cols = ["arm", "start", "end", "n_bins", "logR", "b", "b_lo", "b_hi", "p_imbalance", "call", "confidence", "fit_state"]
    print(segs[["seg"] + cols].to_string(float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
