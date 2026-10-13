#!/usr/bin/env python
"""Genome-wide, per-chromosome and focal BAF / logR plots.
Outputs (output/<sample>/): genome_BAF_logR.png, focal.png, chromosomes/<chrom>.png"""
import json
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from baf_common import BGRID, BIN, CHROMS, parse_args  # noqa: E402
from baf_model import robust_rate  # noqa: E402

# reference palette (light): surface / ink / categorical slots
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
CALL_COL = {"Amplification": "#8f1d1d", "Gain (imbalanced)": "#e34948", "Gain (BAF n.s.)": "#f0a3a2",
            "Gain (balanced)": "#f0a3a2", "Loss/LOH": "#256abf", "Loss (BAF n.s.)": "#9ec5f4",
            "CN-LOH": "#6250d6", "Possible imbalance": "#c9c3f2"}

plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.facecolor": SURF, "figure.facecolor": SURF,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": SURF})

P, ARGS = parse_args(__doc__, lambda ap: ap.add_argument("--focal", default="EGFR;CDKN2A,CDKN2B"))
S = P.sample
bins = pd.read_csv(P.t("bins_1Mb.tsv"), sep="\t")
segs = pd.read_csv(P.o("segments.tsv"), sep="\t")
snps = pd.read_csv(P.t("snps.tsv.gz"), sep="\t")                 # uncapped depths for display
snps = snps[snps.chrom.isin(CHROMS)]
BLL = np.load(P.t("bin_ll.npy"))
params = json.load(open(P.t("params.json")))
purity = json.load(open(P.o("purity.json")))["purity"]
fai = pd.read_csv(P.r("hg19_primary.fa.fai"), sep="\t", header=None, usecols=[0, 1], names=["c", "L"]).set_index("c").L
cyto = pd.read_csv(P.r("cytoBand_hg19.txt.gz"), sep="\t", header=None, names=["chrom", "s", "e", "band", "stain"])
plot_chroms = CHROMS if params["sex"] == "female" else CHROMS[:-1]
off, acc = {}, 0
for c in plot_chroms:
    off[c] = acc
    acc += int(fai[c])
GLEN = acc


def xf_genome(c, p):
    if np.ndim(p):
        return np.asarray(p, float) + np.array([off[x] for x in np.atleast_1d(c)])
    return p + off[c]


def xf_chrom(c, p):
    return np.asarray(p, float) / 1e6


def chrom_axes(ax, ymin, ymax, labels=True):
    for i, c in enumerate(plot_chroms):
        if i % 2 == 1:
            ax.add_patch(Rectangle((off[c], ymin), int(fai[c]), ymax - ymin, color="#f0efec", lw=0, zorder=0))
    ax.set_xlim(0, GLEN)
    ax.set_xticks([off[c] + fai[c] / 2 for c in plot_chroms])
    ax.set_xticklabels([c.replace("chr", "") for c in plot_chroms] if labels else [], fontsize=8)
    ax.tick_params(axis="x", length=0)


def window_baf(b, size=3, step=1):
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


def draw_baf(ax, b, sg, sn, xf, wsize=3, min_depth=8):
    d = sn[sn.depth >= min_depth].assign(vaf=lambda x: x.alt_n / x.depth)
    d = d[(d.vaf > 0.03) & (d.vaf < 0.97)]      # apparent hets only; homozygous SNPs carry no BAF information
    o = d.region == "off"
    ax.scatter(xf(d.chrom[o], d.pos[o]), d.vaf[o], s=4, c=MUTED, alpha=0.35, lw=0, zorder=2)
    ax.scatter(xf(d.chrom[~o], d.pos[~o]), d.vaf[~o], s=9, c=INK2, alpha=0.75, lw=0, zorder=3, marker="D")
    w = window_baf(b, wsize)
    w = w[(w.hi - w.lo) <= 0.35]
    alpha = np.clip(1 - (w.hi - w.lo) / 0.4, 0.15, 0.95).values
    xb = xf(w.chrom, w.x)
    for y in (w.b, 1 - w.b):
        ax.scatter(xb, y, s=10, c=BLUE, alpha=alpha, lw=0, zorder=3)
    for _, r in sg.iterrows():
        x0, x1 = xf(r.chrom, r.start), xf(r.chrom, r.end)
        called = r.call in CALL_COL
        for y in (r.b, 1 - r.b):
            ax.plot([x0, x1], [y, y], color=ORANGE if called else MUTED, lw=2.2 if called else 1.2,
                    solid_capstyle="round", zorder=4)
    ax.axhline(0.5, color=AXIS, lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.set_ylim(-0.02, 1.02)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_ylabel("B-allele frequency")


def draw_logr(ax, b, sg, xf, ylim=(-1.5, 1.6)):
    ok = b.usable
    ax.scatter(xf(b.chrom[ok], (b.start[ok] + b.end[ok]) / 2), b.logR[ok].clip(*ylim), s=6, c=BLUE, alpha=0.45,
               lw=0, zorder=2)
    hi = b[ok & (b.logR > ylim[1])]
    ax.scatter(xf(hi.chrom, (hi.start + hi.end) / 2), np.full(len(hi), ylim[1] - 0.05), marker="^", s=30,
               c=CALL_COL["Amplification"], zorder=5)
    for _, r in sg.iterrows():
        if not np.isnan(r.logR):
            ax.plot([xf(r.chrom, r.start), xf(r.chrom, r.end)], [min(max(r.logR, ylim[0]), ylim[1])] * 2,
                    color=ORANGE, lw=2, solid_capstyle="round", zorder=4)
    ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
    ax.set_ylim(*ylim)
    ax.set_ylabel("log2 ratio (off-target)")


def draw_calls(ax, sg, xf):
    for _, r in sg.iterrows():
        if r.call in CALL_COL:
            x0 = xf(r.chrom, r.start)
            ax.add_patch(Rectangle((x0, 0), xf(r.chrom, r.end) - x0, 1, color=CALL_COL[r.call], lw=0))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_ylabel("Call", rotation=0, ha="right", va="center")


# ---------------- genome-wide ----------------
pur = f"tumour cell content ~{purity:.0%}" if purity is not None else "tumour cell content n/a"
fig, axs = plt.subplots(3, 1, figsize=(16, 8.2), sharex=True, gridspec_kw=dict(height_ratios=[4, 3, 0.45], hspace=0.08))
gb = bins[bins.chrom.isin(plot_chroms)]
draw_baf(axs[0], gb, segs, snps[snps.chrom.isin(plot_chroms)], xf_genome, wsize=5)
chrom_axes(axs[0], -0.02, 1.02, labels=False)
draw_logr(axs[1], gb, segs, xf_genome)
chrom_axes(axs[1], -1.5, 1.6, labels=False)
draw_calls(axs[2], segs, xf_genome)
chrom_axes(axs[2], 0, 1)
fig.suptitle(f"{S} – BAF and copy-number profile (panel, on- + off-target SNPs; hg19; {pur}; "
             f"{params['sex']}, AF: {params['population']})", x=0.125, y=0.965, ha="left", color=INK, fontsize=11)
h = [plt.Line2D([], [], ls="", marker="o", color=BLUE, label="BAF estimate b / 1−b, sliding 5 Mb (step 1 Mb; opacity ~ precision)"),
     plt.Line2D([], [], ls="", marker="D", color=INK2, ms=4, label="Het. SNP on/near-target"),
     plt.Line2D([], [], ls="", marker="o", color=MUTED, ms=3, label="Het. SNP off-target (≥8 reads)"),
     plt.Line2D([], [], color=ORANGE, lw=2, label="Segment with call"),
     plt.Line2D([], [], color=MUTED, lw=1.2, label="Segment without call")]
axs[0].legend(handles=h, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=5, frameon=False, fontsize=8, labelcolor=INK2)
ch = [Rectangle((0, 0), 1, 1, color=v, label=k) for k, v in CALL_COL.items() if k in set(segs.call)]
if ch:
    axs[2].legend(handles=ch, loc="upper left", bbox_to_anchor=(0, -0.9), ncol=len(ch), frameon=False, fontsize=8,
                  labelcolor=INK2)
fig.savefig(P.o("genome_BAF_logR.png"), dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------------- per chromosome ----------------
os.makedirs(os.path.join(P.out, "chromosomes"), exist_ok=True)
for c in plot_chroms:
    b, sg = bins[bins.chrom == c], segs[segs.chrom == c]
    fig, axs = plt.subplots(2, 1, figsize=(9, 5.2), sharex=True, gridspec_kw=dict(height_ratios=[1.2, 1], hspace=0.08))
    draw_baf(axs[0], b, sg, snps[snps.chrom == c], xf_chrom, wsize=3)
    draw_logr(axs[1], b, sg, xf_chrom, ylim=(-1.5, max(1.6, min(np.nanmax(b.logR.values) + 0.3, 6))))
    axs[1].set_xlabel(f"{c} position (Mb)")
    axs[0].set_title(f"{S} – {c}", loc="left", color=INK, fontsize=11)
    for ax in axs:
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_xlim(0, fai[c] / 1e6)
    fig.savefig(os.path.join(P.out, "chromosomes", f"{S}.{c}.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

# ---------------- focal zooms (genes from --focal, if covered by the panel) ----------------
win = pd.read_csv(P.t("windows_10kb.tsv.gz"), sep="\t")
tc = pd.read_csv(P.t("target_coverage.tsv"), sep="\t")
bins["corr"] = np.log2(bins.rate) - bins.logR       # per-bin Alu/GC/centring correction, reused for 250-kb blocks
groups = [[g for g in grp.split(",") if g] for grp in ARGS.focal.split(";")]
groups = [g for g in groups if tc.gene.isin(g).any()]
if groups:
    fig, axs = plt.subplots(2, len(groups), figsize=(4.4 * len(groups), 6.6), squeeze=False,
                            gridspec_kw=dict(height_ratios=[1.2, 1], hspace=0.1, wspace=0.25))
    for j, genes in enumerate(groups):
        gt = tc[tc.gene.isin(genes)]
        c = gt.chrom.iloc[0]
        mid = (gt.start.min() + gt.end.max()) / 2
        x0, x1 = max(mid - 2.5e6, 0), min(mid + 2.5e6, fai[c])
        band = cyto[(cyto.chrom == c) & (cyto.s <= mid) & (cyto.e > mid)].band
        label = "/".join(genes) if len(genes) < 3 else ",".join(genes)
        ww = win[(win.chrom == c) & (win.start >= x0) & (win.start < x1) & win.ok]
        bc = bins[bins.chrom == c].set_index("start")["corr"].interpolate(limit_direction="both")
        rows = []
        for blk, g in ww.groupby(ww.start // 250_000):
            if len(g) >= 10:
                lam, _ = robust_rate(g["count"].values, g.expected.values)
                rows.append(((blk * 250_000 + 125_000) / 1e6,
                             np.log2(max(lam, 1e-3)) - bc.get((blk * 250_000) // BIN * BIN, np.nan)))
        agg = pd.DataFrame(rows, columns=["x", "lr"])
        ax = axs[1, j]
        ax.scatter(agg.x, agg.lr, s=16, c=BLUE, alpha=0.85, lw=0)
        t = tc[(tc.chrom == c) & (tc.start >= x0) & (tc.end <= x1)]
        ax.scatter((t.start + t.end) / 2e6, t.log2_rel, s=24, c=ORANGE, marker="s", lw=0)
        ax.axhline(0, color=AXIS, lw=0.8)
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_xlabel(f"{c} position (Mb)")
        if j == 0:
            ax.set_ylabel("log2 relative depth")
        ax = axs[0, j]
        d = snps[(snps.chrom == c) & (snps.pos >= x0) & (snps.pos < x1) & (snps.depth >= 10)].copy()
        d["vaf"] = d.alt_n / d.depth
        d = d[(d.vaf > 0.02) & (d.vaf < 0.98)]
        ax.scatter(d.pos / 1e6, d.vaf, s=10 + 2 * np.sqrt(d.depth), c=np.where(d.region == "off", MUTED, INK2),
                   alpha=0.8, lw=0)
        ax.axhline(0.5, color=AXIS, lw=0.8, ls=(0, (3, 3)))
        ax.set_ylim(-0.02, 1.02)
        if j == 0:
            ax.set_ylabel("VAF of heterozygous SNPs (≥10 reads)")
        for a_ in axs[:, j]:
            a_.axvspan(gt.start.min() / 1e6 - 0.02, gt.end.max() / 1e6 + 0.02, color="#f0efec", zorder=0)
            a_.set_xlim(x0 / 1e6, x1 / 1e6)
        ax.set_title(f"{label} ({c.replace('chr', '')}{band.iloc[0] if len(band) else ''})", loc="left",
                     color=INK, fontsize=10)
    hl = [plt.Line2D([], [], ls="", marker="o", color=INK2, label="SNP on/near-target (size ~ depth)"),
          plt.Line2D([], [], ls="", marker="o", color=MUTED, label="SNP off-target"),
          plt.Line2D([], [], ls="", marker="o", color=BLUE, label="Off-target depth, 250-kb blocks (GC/Alu-corrected)"),
          plt.Line2D([], [], ls="", marker="s", color=ORANGE, label="On-target depth per target (GC-corrected, rel. panel median)")]
    fig.legend(handles=hl, loc="lower left", bbox_to_anchor=(0.08, -0.05), ncol=2, frameon=False, fontsize=8,
               labelcolor=INK2)
    fig.savefig(P.o("focal.png"), dpi=170, bbox_inches="tight")
    plt.close(fig)
print("plots written")
