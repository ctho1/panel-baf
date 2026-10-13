#!/usr/bin/env python
"""
BAF + read-depth copy-number profile from a hybrid-capture panel (tumour only, no matched normal).

BAF   Common SNVs (1000G, MAF >= 5 %), allele counts from on-, near- and off-target reads. Off-target depth
      is typically ~0.1x, so genotypes cannot be called per SNP: every SNP with >= 2 reads contributes a
      likelihood over b, the minor-allele fraction of heterozygous SNPs (see baf_model.snp_loglik).
      Genotype priors use the 1000G superpopulation chosen by likelihood; the het rate h is fitted.
logR  Off-target fragments (targets +-1 kb masked) in 10 kb windows (mappability >= 0.9), GC-expected,
      capture hotspots removed, aggregated to 1 Mb bins, corrected for Alu density + GC (off-target capture
      is strongly Alu-driven), centred on the autosomal mode.
Segmentation  Joint HMM (Viterbi) over (logR, b) per chromosome arm.
Outputs (tmp/<sample>/): bins_1Mb.tsv, segments.tsv, arms.tsv, *_ll.npy, informative_snps.tsv.gz, params.json
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

from baf_common import AUTOSOMES, BGRID, BIN, CHROMS, REF_FASTA, bci, load_arms, load_wig, parse_args, read_bed, read_fai
from baf_model import DELTA, EPS, RHO, lowess, profiled_ll, robust_rate, snp_loglik, thin, viterbi_segment

WIN = 10_000
MASK_PAD = 1000
THIN_BP = 150
OFF_CAP = 15          # off-target SNPs: cap effective depth (capture hotspots / paralog mismapping)
MIN_UNIQUE_FRAC = 0.5      # min. fraction of reads with MAPQ >= 20 per 10-kb window
MIN_UNIQUE_FRAC_BIN = 0.8  # ... and per 1-Mb bin (edges of alt-contig / segmental-duplication regions)
POPS = ["AF", "EUR", "AFR", "AMR", "EAS", "SAS"]
JUMP_LOGP = -15.0


def build_bins(P, chroms):
    fai = read_fai(P.r(REF_FASTA))["len"]
    gcw, mpw = load_wig(P.r("gc_hg19_10kb.wig")), load_wig(P.r("map_hg19_10kb.wig"))
    bed = read_bed(P.bed)
    alu = pd.read_csv(P.r("alu_hg19.bed.gz"), sep="\t", header=None, names=["chrom", "s", "e"])
    frag = np.load(P.t("fragment_midpoints.npz"))
    win, alu_bin = [], {}
    for c in chroms:
        L, n = int(fai[c]), len(gcw[c])
        cnt = np.bincount(frag[c] // WIN, minlength=n)[:n]
        cnt_all = np.bincount(frag[c + "_all"] // WIN, minlength=n)[:n]
        masked = np.zeros(n, bool)
        sub = bed[bed.chrom == c]
        for s, e in zip(sub.s - MASK_PAD, sub.e + MASK_PAD):
            masked[max(s, 0) // WIN:(e - 1) // WIN + 1] = True
        am = np.zeros(L, bool)
        a = alu[alu.chrom == c]
        for s, e in zip(a.s.values, a.e.values):
            am[s:e] = True
        starts = np.arange(0, L, BIN)
        alu_bin[c] = np.add.reduceat(am, starts) / np.diff(np.r_[starts, L])
        win.append(pd.DataFrame(dict(chrom=c, start=np.arange(n) * WIN, count=cnt, count_all=cnt_all, gc=gcw[c],
                                     map=mpw[c], masked=masked)))
    win = pd.concat(win, ignore_index=True)
    # reads that cannot be placed uniquely (MAPQ < 20 for most reads: alt haplotype contigs such as the MHC,
    # segmental duplications) make a window unusable even if its reference mappability looks fine
    unique = (win.count_all == 0) | (win["count"] >= MIN_UNIQUE_FRAC * win.count_all)
    win["ok"] = (win["map"] >= 0.9) & (win.gc > 0.25) & (win.gc < 0.7) & ~win.masked & unique
    au = win[win.ok & (win.chrom != "chrX")]
    strata = [(g, robust_rate(grp["count"].values, np.ones(len(grp)))[0])
              for g, grp in au.groupby(np.round(au.gc / 0.005) * 0.005) if len(grp) >= 200]
    strata = np.array(strata)
    fit = lowess(np.log(strata[:, 1]), strata[:, 0], frac=0.4)
    win["expected"] = np.exp(np.interp(win.gc, fit[:, 0], fit[:, 1]))
    win["bin"] = win.start // BIN
    rows = []
    for (c, bi), grp in win.groupby(["chrom", "bin"], sort=False):
        g = grp[grp.ok]
        lam, nkeep = robust_rate(g["count"].values, g.expected.values) if len(g) >= 30 else (np.nan, 0)
        n_all = int(grp.count_all.sum())
        rows.append((c, bi * BIN, min((bi + 1) * BIN, int(fai[c])), int(grp["count"].sum()), len(g), nkeep, lam,
                     float(g.gc.mean()) if len(g) else np.nan, float(grp["map"].mean()),
                     float(alu_bin[c][bi]) if bi < len(alu_bin[c]) else np.nan,
                     float(grp["count"].sum() / n_all) if n_all else np.nan))
    bins = pd.DataFrame(rows, columns=["chrom", "start", "end", "count_raw", "n_win_ok", "n_win_used", "rate",
                                       "gc", "map", "alu_frac", "unique_frac"])
    win.drop(columns=["bin"]).to_csv(P.t("windows_10kb.tsv.gz"), sep="\t", index=False, float_format="%.4g")
    return bins


def correct_logR(bins):
    """Backfitting LOWESS of log2 rate on Alu fraction and GC (autosomes); centre on the autosomal mode."""
    ok = (bins.rate.notna() & (bins.rate > 0) & (bins.n_win_used >= 30)
          & (bins.unique_frac.fillna(0) >= MIN_UNIQUE_FRAC_BIN))   # e.g. MHC with alt haplotype contigs
    bins["usable"] = ok
    auto = (ok & bins.chrom.isin(AUTOSOMES)).values
    y = np.full(len(bins), np.nan)
    y[ok.values] = np.log2(bins.rate.values[ok.values].astype(float))
    mu = np.nanmedian(y[auto])
    f_alu = np.zeros(len(bins)); f_gc = np.zeros(len(bins))
    for _ in range(6):
        fa = lowess(y[auto] - mu - f_gc[auto], bins.alu_frac.values[auto], frac=0.4, it=3)
        f_alu = np.interp(bins.alu_frac.values, fa[:, 0], fa[:, 1]); f_alu -= np.median(f_alu[auto])
        fg = lowess(y[auto] - mu - f_alu[auto], bins.gc.values[auto], frac=0.4, it=3)
        f_gc = np.interp(bins.gc.values, fg[:, 0], fg[:, 1]); f_gc -= np.median(f_gc[auto])
    lr = y - mu - f_alu - f_gc
    kde = stats.gaussian_kde(lr[auto], bw_method=0.1)
    xs = np.linspace(-1, 1, 2001)
    bins["logR"] = np.where(ok, lr - xs[np.argmax(kde(xs))], np.nan)
    return bins


def infer_sex(P):
    """Male if fragment density in the male-specific chrY region reaches >= 15 % of the autosomal median."""
    frag = np.load(P.t("fragment_midpoints.npz"))
    auto = np.concatenate([np.bincount(frag[c] // BIN) for c in AUTOSOMES])
    y = frag["chrY"]
    y_density = ((y > 2_800_000) & (y < 28_800_000)).sum() / 26.0
    ratio = float(y_density / max(np.median(auto), 1))
    return ("male" if ratio >= 0.15 else "female"), ratio


def estimate_ref_bias(snps):
    """alt/ref sampling ratio from deep on/near-target heterozygous SNPs (median VAF)."""
    d = snps[(snps.region != "off") & (snps.depth >= 100)]
    vaf = d.alt_n / d.depth
    het = vaf[(vaf > 0.3) & (vaf < 0.7)]
    if len(het) < 20:
        return 0.96, len(het)
    m = float(np.median(het))
    return float(np.clip(m / (1 - m), 0.9, 1.05)), len(het)


def summarise(bins, bll, idx, **extra):
    lr = bins.logR.values[idx]; lr = lr[~np.isnan(lr)]
    sll = bll[idx].sum(axis=0)
    bm, blo, bhi = bci(sll)
    lrt = max(2 * (sll.max() - sll[-1]), 0.0)
    return dict(**extra, start=int(bins.start.values[idx[0]]), end=int(bins.end.values[idx[-1]]),
                n_bins=len(idx), n_bins_logR=len(lr),
                logR=float(np.median(lr)) if len(lr) else np.nan,
                logR_se=float(1.2533 * np.std(lr) / np.sqrt(len(lr))) if len(lr) > 2 else np.nan,
                n_snp_inf=int(bins.n_snp_inf.values[idx].sum()),
                n_snp_deep=int(bins.n_snp_deep.values[idx].sum()),
                b=bm, b_lo=blo, b_hi=bhi, LRT_imbalance=lrt,
                p_imbalance=float(0.5 * stats.chi2.sf(lrt, 1)) if lrt > 0 else 1.0), sll


def main():
    P, _ = parse_args(__doc__)
    arms = load_arms(P.ref)
    sex, y_ratio = infer_sex(P)
    seg_chroms = CHROMS if sex == "female" else AUTOSOMES
    bins = correct_logR(build_bins(P, CHROMS))
    bins["arm"] = [c + ("p" if s + BIN / 2 < arms[c][0] else "q") for c, s in zip(bins.chrom, bins.start)]
    auto = bins.usable & bins.chrom.isin(AUTOSOMES)
    d = np.diff(bins.logR[auto].values)
    sigma = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)

    # ---- SNPs, reference bias, population AF, het-rate calibration ----
    snps = pd.read_csv(P.t("snps.tsv.gz"), sep="\t")
    snps = snps[snps.chrom.isin(seg_chroms)]
    snps["arm"] = [c + ("p" if p < arms[c][0] else "q") for c, p in zip(snps.chrom, snps.pos)]
    w_ref, n_w = estimate_ref_bias(snps)
    inf = thin(snps[snps.depth >= 2], THIN_BP).reset_index(drop=True)   # depth-1 SNPs carry no information on b
    cap = (inf.region == "off") & (inf.depth > OFF_CAP)
    inf.loc[cap, "alt_n"] = np.round(inf.alt_n[cap] * OFF_CAP / inf.depth[cap]).astype(int)
    inf.loc[cap, "ref_n"] = OFF_CAP - inf.alt_n[cap]
    inf.loc[cap, "depth"] = OFF_CAP
    auto_inf = inf[inf.chrom.isin(AUTOSOMES)]
    pop_ll = {pop: profiled_ll(auto_inf, pop, 1.0, auto_inf.arm.values, w_ref) for pop in POPS}
    pop = max(pop_ll, key=pop_ll.get)
    hgrid = np.round(np.arange(0.70, 1.151, 0.025), 3)
    h_ll = np.array([profiled_ll(auto_inf, pop, h, auto_inf.arm.values, w_ref) for h in hgrid])
    h = float(hgrid[np.argmax(h_ll)])
    ll = snp_loglik(inf.alt_n.values, inf.depth.values, inf[pop].values, h, w_ref)
    ll -= ll.max(axis=1, keepdims=True)
    inf["bin"] = inf.pos.values // BIN
    key_b = pd.MultiIndex.from_arrays([bins.chrom, bins.start // BIN])
    row = key_b.get_indexer(pd.MultiIndex.from_arrays([inf.chrom, inf.bin]))
    bll = np.zeros((len(bins), len(BGRID)))
    np.add.at(bll, row, ll)
    bins["n_snp_inf"] = np.bincount(row, minlength=len(bins))
    bins["n_snp_deep"] = np.bincount(row, weights=(inf.depth.values >= 30), minlength=len(bins)).astype(int)
    est = np.array([bci(x) if n > 0 else (np.nan,) * 3 for x, n in zip(bll, bins.n_snp_inf)])
    bins["b_mle"], bins["b_lo"], bins["b_hi"] = est[:, 0], est[:, 1], est[:, 2]
    bins["b_lrt"] = [max(2 * (x.max() - x[-1]), 0) for x in bll]

    # ---- joint segmentation per chromosome arm ----
    mu_grid = np.round(np.arange(-2.5, 7.01, 0.05), 2)
    b_idx = np.arange(0, len(BGRID), 2)
    bins["seg"] = -1
    seg_rows, seg_ll, arm_rows, arm_ll = [], [], [], []
    for c in seg_chroms:
        for a in ("p", "q"):
            ix = np.where((bins.arm == c + a).values & (bins.usable.values | (bins.n_snp_inf.values > 0)))[0]
            if len(ix) < 3:
                continue
            r, l = summarise(bins, bll, ix, chrom=c, arm=c + a)
            arm_rows.append(r); arm_ll.append(l)
            path = viterbi_segment(bins.logR.values[ix], bll[ix], sigma, mu_grid, b_idx, JUMP_LOGP)
            brk = np.r_[0, np.where(np.diff(path) != 0)[0] + 1, len(ix)]
            for s0, s1 in zip(brk[:-1], brk[1:]):
                sid = len(seg_rows)
                bins.loc[bins.index[ix[s0:s1]], "seg"] = sid
                r, l = summarise(bins, bll, ix[s0:s1], seg=sid, chrom=c, arm=c + a)
                seg_rows.append(r); seg_ll.append(l)
    segs, armt = pd.DataFrame(seg_rows), pd.DataFrame(arm_rows)

    bins.to_csv(P.t("bins_1Mb.tsv"), sep="\t", index=False, float_format="%.4g")
    segs.to_csv(P.t("segments.tsv"), sep="\t", index=False, float_format="%.4g")
    armt.to_csv(P.t("arms.tsv"), sep="\t", index=False, float_format="%.4g")
    np.save(P.t("bin_ll.npy"), bll)
    np.save(P.t("seg_ll.npy"), np.array(seg_ll))
    np.save(P.t("arm_ll.npy"), np.array(arm_ll))
    inf.to_csv(P.t("informative_snps.tsv.gz"), sep="\t", index=False)
    params = dict(sex=sex, chrY_density_ratio=round(y_ratio, 4), sigma_logR_bin=round(float(sigma), 4),
                  population=pop, pop_loglik={k: round(v, 1) for k, v in pop_ll.items()},
                  het_rate_scale=h, h_profile=dict(zip(hgrid.tolist(), np.round(h_ll - h_ll.max(), 2).tolist())),
                  w_ref=round(w_ref, 4), n_hets_for_w_ref=n_w,
                  n_snps_covered=len(snps), n_inf_snps=len(inf), n_inf_snps_offtarget=int((inf.region == "off").sum()),
                  rho=RHO, eps=EPS, delta=DELTA, jump_logp=JUMP_LOGP, off_cap=OFF_CAP)
    json.dump(params, open(P.t("params.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in params.items() if k != "h_profile"}, indent=1))


if __name__ == "__main__":
    main()
