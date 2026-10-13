#!/usr/bin/env python
"""Genome-wide BAF / copy-number figure (the pipeline's only deliverable).
Inputs: tmp/<sample>/ (bins, bin log-likelihoods, SNPs, called segments, purity)
Outputs: output/<sample>/<sample>.genome_BAF_logR.png (600 dpi) and .pdf (vector, editable text)"""
import json

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")

from baf_common import BGRID, CHROMS, REF_FASTA, parse_args, read_fai  # noqa: E402
import figure_genome  # noqa: E402


def _extra_args(ap):
    ap.add_argument("--width-mm", type=float, default=figure_genome.STYLE["width_mm"],
                    help="final print width (174 Springer double column, 180 Nature)")
    ap.add_argument("--height-mm", type=float, default=figure_genome.STYLE["height_mm"])


P, ARGS = parse_args(__doc__, _extra_args)
bins = pd.read_csv(P.t("bins_1Mb.tsv"), sep="\t")
segs = pd.read_csv(P.t("segments_called.tsv"), sep="\t")
snps = pd.read_csv(P.t("snps.tsv.gz"), sep="\t")                 # uncapped depths for display
snps = snps[snps.chrom.isin(CHROMS)]
BLL = np.load(P.t("bin_ll.npy"))
params = json.load(open(P.t("params.json")))
purity = json.load(open(P.t("purity.json")))["purity"]
fai = read_fai(P.r(REF_FASTA))["len"]
cyto = pd.read_csv(P.r("cytoBand_hg19.txt.gz"), sep="\t", header=None, names=["chrom", "s", "e", "band", "stain"])
CEN = cyto[cyto.stain == "acen"].groupby("chrom").apply(lambda g: (g.s.min() + g.e.max()) / 2).to_dict()
plot_chroms = CHROMS if params["sex"] == "female" else CHROMS[:-1]
off, acc = {}, 0
for c in plot_chroms:
    off[c] = acc
    acc += int(fai[c])


def xf_genome(c, p):
    return np.asarray(p, float) + np.array([off[x] for x in np.atleast_1d(c)])


def window_baf(b, size=5, step=1):
    """Sliding windows of `size` 1-Mb bins (step `step`) within each arm: summed log-likelihood -> MLE + 95 % CI."""
    rows = []
    for _, g in b[b.n_snp_inf > 0].groupby("arm", sort=False):
        idx = g.index.values
        for k in range(0, max(len(idx) - size, 0) + 1, step):
            ii = idx[k:k + size]
            ll = BLL[ii].sum(axis=0)
            i = int(np.argmax(ll))
            ok = ll >= ll[i] - 1.92
            rows.append(dict(chrom=g.chrom.iloc[0], x=(b.start[ii].min() + b.end[ii].max()) / 2,
                             b=BGRID[i], lo=BGRID[ok].min(), hi=BGRID[ok].max()))
    return pd.DataFrame(rows, columns=["chrom", "x", "b", "lo", "hi"])


gb = bins[bins.chrom.isin(plot_chroms)]
w = window_baf(gb, 5)
w = w[(w.hi - w.lo) <= 0.35]                                   # drop uninformative windows
het = snps[snps.chrom.isin(plot_chroms) & (snps.depth >= 8)].assign(vaf=lambda x: x.alt_n / x.depth)
het = het[(het.vaf > 0.03) & (het.vaf < 0.97)]                 # apparent hets only
baf_points = pd.DataFrame(dict(
    x=np.r_[xf_genome(w.chrom, w.x), xf_genome(w.chrom, w.x), xf_genome(het.chrom, het.pos)],
    y=np.r_[w.b, 1 - w.b, het.vaf]))
ok = gb.usable
log2_points = pd.DataFrame(dict(x=xf_genome(gb.chrom[ok], (gb.start[ok] + gb.end[ok]) / 2), y=gb.logR[ok]))

seg_cols = ["chrom", "start", "end", "b", "logR", "call"]
seg_before = segs[seg_cols].copy(deep=True)
drawn, paths = figure_genome.draw_genome_figure(
    chroms=plot_chroms, offsets=off, chrom_len={c: int(fai[c]) for c in plot_chroms}, centromeres=CEN,
    baf_points=baf_points, log2_points=log2_points, segments=segs[seg_cols], purity=purity,
    out_prefix=P.o("genome_BAF_logR"), width_mm=ARGS.width_mm, height_mm=ARGS.height_mm)
# presentation only: drawn segment values, boundaries and calls must equal the segmentation output
pd.testing.assert_frame_equal(drawn[seg_cols].reset_index(drop=True), seg_before.reset_index(drop=True))
print("figure:", ", ".join(paths.values()))
