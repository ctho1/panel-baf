#!/usr/bin/env python
"""On-target depth per gene (GC-corrected, relative to the panel median) plus heterozygous on-target SNPs per
gene. Without a panel of normals, target-to-target capture efficiency is uncorrected -> only large deviations
(amplifications, homozygous deletions) are interpretable.
Outputs: output/<sample>/<sample>.gene_coverage.tsv, tmp/<sample>/<sample>.target_coverage.tsv"""
import numpy as np
import pandas as pd
import pysam
from statsmodels.nonparametric.smoothers_lowess import lowess

from baf_common import parse_args

P, _ = parse_args(__doc__)
t = pd.read_csv(P.t("target.regions.bed.gz"), sep="\t", header=None)
t = t.iloc[:, [0, 1, 2, t.shape[1] - 1]]
t.columns = ["chrom", "start", "end", "depth"]
fa = pysam.FastaFile(P.r("hg19_primary.fa"))
seqs = [fa.fetch(c, s, e).upper() for c, s, e in zip(t.chrom, t.start, t.end)]
t["gc"] = [(s.count("G") + s.count("C")) / max(len(s), 1) for s in seqs]
t["len"] = t.end - t.start

ref = pd.read_csv(P.r("refGene_hg19.txt.gz"), sep="\t", header=None, usecols=[2, 4, 5, 12],
                  names=["chrom", "txs", "txe", "gene"])


def pick(names):
    """Prefer a proper gene symbol over LOC/MIR/LINC/antisense entries, then the shortest name."""
    return sorted(names, key=lambda x: (x.startswith(("LOC", "MIR", "LINC")) or "-AS" in x, len(x)))[0]


genes = []
for c, s, e in zip(t.chrom, t.start, t.end):
    r = ref[(ref.chrom == c) & (ref.txs <= e + 2000) & (ref.txe >= s - 2000)]
    genes.append(pick(sorted(set(r.gene))) if len(r) else "intergenic")
t["gene"] = genes

lr = np.log2(t.depth.clip(lower=1) / t.depth.median())
fit = lowess(lr, t.gc, frac=0.5, it=3, return_sorted=True)
t["log2_rel"] = lr - np.interp(t.gc, fit[:, 0], fit[:, 1])
sd = 1.4826 * np.median(np.abs(t.log2_rel - t.log2_rel.median()))

snps = pd.read_csv(P.t("snps.tsv.gz"), sep="\t")
deep = snps[(snps.region != "off") & (snps.depth >= 50)].copy()
deep["vaf"] = deep.alt_n / deep.depth
het = deep[(deep.vaf > 0.05) & (deep.vaf < 0.95)]

rows = []
for (c, g), sub in t.groupby(["chrom", "gene"], sort=False):
    s0, s1 = sub.start.min() - 1000, sub.end.max() + 1000
    h = het[(het.chrom == c) & (het.pos >= s0) & (het.pos <= s1)]
    mb = np.minimum(h.vaf, 1 - h.vaf)
    rows.append(dict(chrom=c, gene=g, start=int(sub.start.min()), end=int(sub.end.max()), n_targets=len(sub),
                     target_bp=int(sub.len.sum()), mean_depth=float(np.average(sub.depth, weights=sub.len)),
                     log2_rel=float(np.median(sub.log2_rel)), z=float(np.median(sub.log2_rel) / (sd / np.sqrt(len(sub)))),
                     n_het_snps=len(h), median_minor_af=float(np.median(mb)) if len(h) else np.nan))
pd.DataFrame(rows).to_csv(P.o("gene_coverage.tsv"), sep="\t", index=False, float_format="%.4g")
t.to_csv(P.t("target_coverage.tsv"), sep="\t", index=False, float_format="%.4g")
print(f"target-level robust SD of GC-corrected log2 depth: {sd:.3f}")
