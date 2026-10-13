"""Publication-ready genome-wide BAF / copy-number figure (presentation only).

Drawn at final print size (default 174 mm x 95 mm, Springer double column; 180 mm for Nature) without any
downscaling or tight-bbox cropping. Three aligned panels share the genomic x axis:
B-allele frequency, copy number (log2 ratio) and the call track.
rcParams live in styles/publication.mplstyle; figure-specific settings (colours, widths, layout) in STYLE.
The function never modifies its inputs and returns the table of drawn segments so callers can assert that
values, boundaries and calls are exactly those of the segmentation output.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle
from matplotlib.transforms import blended_transform_factory, offset_copy

STYLE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "styles", "publication.mplstyle")

MM = 1 / 25.4
STYLE = dict(
    width_mm=174.0, height_mm=95.0,
    margin_mm=dict(left=11.0, right=12.0, top=2.5, bottom=9.6),
    gap_ab_mm=2.2, gap_bc_mm=1.6, height_ratios=(1.0, 1.0, 0.15),
    point_color="#6e6e6e", point_size=1.2, point_alpha=0.35,
    # BAF panel: called segments dark, uncalled mid grey + thinner (mirrored b and 1 - b)
    seg_called="#1b1b1b", seg_called_lw=1.3, seg_uncalled="#9a9a9a", seg_uncalled_lw=0.7,
    # log2 panel: muted blue
    log2_called="#3b6ea5", log2_uncalled="#a9bfd9",
    ref_color="#7a7a7a", ref_lw=0.4,
    chrom_line="#b4b4b4", chrom_lw=0.35,
    centromere="#d2d2d2", centromere_lw=0.35, centromere_ls=(0, (0.6, 1.4)),
    expected_color="#a8a8a8", expected_lw=0.4, expected_ls=(0, (2.2, 1.8)),
    offscale_color="#1b1b1b", offscale_ms=3.2, offscale_fontsize=5.5, call_min_width_mm=0.4,
    margin_label_fontsize=5.5,
    baf_ylim=(-0.01, 1.01), log2_ylim=(-1.5, 1.5),
)
# (call value in segments table, figure label, face colour, hatch) -- legend order: gains, losses, LOH, other.
# Red family = gains, blue family = losses, distinct purple = CN-LOH (checked for deuteranopia/protanopia:
# smallest OKLab dE between any two fills 12.2 deut / 9.7 prot); possible imbalance = light hatching.
CALLS = [
    ("Amplification", "Amplification", "#67001F", None),
    ("Gain (imbalanced)", "Gain", "#D6604D", None),
    ("Gain (BAF n.s.)", "Gain (log2 only)", "#F4A582", None),
    ("Gain (balanced)", "Gain (balanced alleles)", "#FDDBC7", None),
    ("Loss/LOH", "Loss with LOH", "#2166AC", None),
    ("Loss (BAF n.s.)", "Loss (log2 only)", "#92C5DE", None),
    ("CN-LOH", "CN-LOH", "#AA3377", None),
    ("Possible imbalance", "Possible imbalance", "white", "//////"),
]
CALL_MAP = {c[0]: c for c in CALLS}


def expected_levels(purity):
    """Expected log2 ratio and minor-allele BAF for 1, 3 and 4 copies (3:1) at the given tumour cell content,
    diploid baseline, no compression."""
    p = purity
    return [  # (copies, label log2, log2, label BAF, minor BAF)
        (1, "CN 1", np.log2(1 - p / 2), "1:0", (1 - p) / (2 - p)),
        (3, "CN 3", np.log2(1 + p / 2), "2:1", 1 / (2 + p)),
        (4, "CN 4", np.log2(1 + p), "3:1", 1 / (2 + 2 * p)),
    ]


def _repel(values, min_gap):
    """Spread 1-D label positions so neighbours are at least min_gap apart (keeps the mean)."""
    order = np.argsort(values)
    v = np.asarray(values, float)[order]
    for _ in range(50):
        moved = False
        for i in range(1, len(v)):
            if v[i] - v[i - 1] < min_gap:
                d = (min_gap - (v[i] - v[i - 1])) / 2
                v[i - 1] -= d
                v[i] += d
                moved = True
        if not moved:
            break
    out = np.empty_like(v)
    out[order] = v
    return out


def _margin_labels(ax, ys, labels, st, panel_h_mm):
    """Right-margin labels for horizontal reference lines, de-overlapped, with short leader lines."""
    lo, hi = ax.get_ylim()
    to_mm = lambda y: (y - lo) / (hi - lo) * panel_h_mm
    pos_mm = _repel([to_mm(y) for y in ys], min_gap=st["margin_label_fontsize"] * 0.3528 * 1.25)
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    for y, ymm, lab in zip(ys, pos_mm, labels):
        y_lab = lo + ymm / panel_h_mm * (hi - lo)
        ax.plot([1.0, 1.006, 1.016], [y, y_lab, y_lab], transform=tr, color=st["expected_color"], lw=0.35,
                clip_on=False)
        ax.text(1.019, y_lab, lab, transform=tr, va="center", ha="left", fontsize=st["margin_label_fontsize"],
                color="#3c3c3c", clip_on=False)


def draw_genome_figure(*, chroms, offsets, chrom_len, centromeres, baf_points, log2_points, segments, purity,
                       out_prefix, width_mm=None, height_mm=None, style=None):
    """Draw and save <out_prefix>.png (600 dpi) and <out_prefix>.pdf. Inputs (not modified):
      baf_points  DataFrame x (genome coordinate), y (BAF)        -- windows (b, 1-b) and single het SNPs
      log2_points DataFrame x, y (log2 ratio per 1-Mb bin)
      segments    DataFrame chrom, start, end, b, logR, call
    Returns (drawn segment table, {format: path})."""
    st = dict(STYLE, **(style or {}))
    W = width_mm or st["width_mm"]
    H = height_mm or st["height_mm"]
    m = st["margin_mm"]
    glen = offsets[chroms[-1]] + chrom_len[chroms[-1]]
    gx = lambda c, p: offsets[c] + p
    seg_in = segments.copy(deep=True)

    plt.style.use(STYLE_FILE)
    fig = plt.figure(figsize=(W * MM, H * MM))
    usable_h = H - m["top"] - m["bottom"] - st["gap_ab_mm"] - st["gap_bc_mm"]
    hr = np.array(st["height_ratios"]) / sum(st["height_ratios"]) * usable_h
    pw = W - m["left"] - m["right"]
    rect = lambda y0_mm, h_mm: [m["left"] / W, y0_mm / H, pw / W, h_mm / H]
    y_c = m["bottom"]
    y_b = y_c + hr[2] + st["gap_bc_mm"]
    y_a = y_b + hr[1] + st["gap_ab_mm"]
    ax_a = fig.add_axes(rect(y_a, hr[0]))
    ax_b = fig.add_axes(rect(y_b, hr[1]), sharex=ax_a)
    ax_c = fig.add_axes(rect(y_c, hr[2]), sharex=ax_a)

    # chromosome boundaries + faint centromeres
    for ax in (ax_a, ax_b, ax_c):
        for c in chroms[1:]:
            ax.axvline(offsets[c], color=st["chrom_line"], lw=st["chrom_lw"], zorder=0.5)
    for ax in (ax_a, ax_b):
        for c in chroms:
            if c in centromeres:
                ax.axvline(gx(c, centromeres[c]), color=st["centromere"], lw=st["centromere_lw"],
                           ls=st["centromere_ls"], zorder=0.4)

    # points (rasterized: the PDF stays small, text and segments stay vector)
    ax_a.scatter(baf_points.x, baf_points.y, s=st["point_size"], c=st["point_color"], alpha=st["point_alpha"],
                 lw=0, rasterized=True, zorder=1)
    lo, hi = st["log2_ylim"]
    inr = log2_points[(log2_points.y >= lo) & (log2_points.y <= hi)]
    ax_b.scatter(inr.x, inr.y, s=st["point_size"], c=st["point_color"], alpha=st["point_alpha"], lw=0,
                 rasterized=True, zorder=1)

    # reference lines
    ax_a.axhline(0.5, color=st["ref_color"], lw=st["ref_lw"], zorder=0.6)
    ax_b.axhline(0.0, color=st["ref_color"], lw=st["ref_lw"], zorder=0.6)

    # purity-adjusted expected levels (faint dashed) with right-margin labels
    if purity is not None:
        lev = expected_levels(purity)
        for _, _, lr, _, b in lev:
            ax_b.axhline(lr, color=st["expected_color"], lw=st["expected_lw"], ls=st["expected_ls"], zorder=0.7)
            for y in (b, 1 - b):
                ax_a.axhline(y, color=st["expected_color"], lw=st["expected_lw"], ls=st["expected_ls"], zorder=0.7)
        _margin_labels(ax_b, [l[2] for l in lev], [l[1] for l in lev], st, hr[1])
        _margin_labels(ax_a, [1 - l[4] for l in lev], [l[3] for l in lev], st, hr[0])

    # segments: BAF mirrored (b, 1 - b) dark/grey; log2 muted blue (uncalled lighter + thinner)
    drawn = []
    for r in seg_in.itertuples(index=False):
        called = r.call in CALL_MAP
        lw = st["seg_called_lw"] if called else st["seg_uncalled_lw"]
        z = 3 if called else 2
        x0, x1 = gx(r.chrom, r.start), gx(r.chrom, r.end)
        if not np.isnan(r.b):
            col = st["seg_called"] if called else st["seg_uncalled"]
            ax_a.plot([x0, x1], [r.b, r.b], color=col, lw=lw, zorder=z)
            ax_a.plot([x0, x1], [1 - r.b, 1 - r.b], color=col, lw=lw, zorder=z)
        if not np.isnan(r.logR) and lo <= r.logR <= hi:
            ax_b.plot([x0, x1], [r.logR, r.logR], color=st["log2_called"] if called else st["log2_uncalled"],
                      lw=lw, zorder=z)
        drawn.append(dict(chrom=r.chrom, start=r.start, end=r.end, b=r.b, logR=r.logR, call=r.call, called=called))

    # off-scale log2 values: triangle at the panel edge + actual value (one marker per cluster of bins)
    for side, sel in (("hi", log2_points.y > hi), ("lo", log2_points.y < lo)):
        pts = log2_points[sel].sort_values("x")
        if pts.empty:
            continue
        cluster = np.cumsum(np.r_[True, np.diff(pts.x.values) > 5e6])
        for _, g in pts.groupby(cluster):
            k = g.y.idxmax() if side == "hi" else g.y.idxmin()
            x, v = g.x[k], g.y[k]
            yv = hi if side == "hi" else lo
            ax_b.plot([x], [yv], marker="^" if side == "hi" else "v", ms=st["offscale_ms"], mew=0,
                      color=st["offscale_color"], clip_on=False, zorder=4)
            ax_b.annotate(f"{v:.1f}".replace("-", "−"), (x, yv), xytext=(2.2, -1.2 if side == "hi" else 3.0),
                          textcoords="offset points", fontsize=st["offscale_fontsize"], ha="left", va="center",
                          color=st["offscale_color"])

    # call track (focal events get a minimum display width so they stay visible in print)
    min_w = st["call_min_width_mm"] / pw * glen
    for r in seg_in.itertuples(index=False):
        if r.call in CALL_MAP:
            _, _, fc, hatch = CALL_MAP[r.call]
            x0, x1 = gx(r.chrom, r.start), gx(r.chrom, r.end)
            focal = x1 - x0 < min_w
            if focal:
                x0, x1 = (x0 + x1 - min_w) / 2, (x0 + x1 + min_w) / 2
            ax_c.add_patch(Rectangle((x0, 0), x1 - x0, 1, facecolor=fc, hatch=hatch,
                                     edgecolor="#8c8c8c" if hatch else "none", lw=0, zorder=2 if focal else 1))

    # axes cosmetics
    ax_a.set_ylim(*st["baf_ylim"])
    ax_a.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax_a.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
    ax_a.set_ylabel("B-allele frequency")
    ax_b.set_ylim(lo, hi)
    ax_b.set_yticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5])
    ax_b.set_yticklabels(["−1.5", "−1", "−0.5", "0", "0.5", "1", "1.5"])
    ax_b.set_ylabel("Copy number (log2 ratio)")
    ax_c.set_ylim(0, 1)
    ax_c.set_yticks([])
    ax_c.spines["left"].set_visible(False)
    ax_c.set_ylabel("Call", rotation=0, ha="right", va="center", labelpad=4)
    ax_a.set_xlim(0, glen)
    for ax in (ax_a, ax_b, ax_c):
        ax.tick_params(axis="x", which="both", length=0, labelbottom=False)
    ax_c.spines["bottom"].set_visible(False)
    fig.align_ylabels([ax_a, ax_b])

    # chromosome labels centred per chromosome; narrow chromosomes alternate between two rows
    base = blended_transform_factory(ax_c.transData, ax_c.transAxes)
    narrow = [c for c in chroms if chrom_len[c] / glen * pw < 3.4]
    row = {c: i % 2 for i, c in enumerate(narrow)}
    for c in chroms:
        tr = offset_copy(base, fig=fig, y=-1.5 - 6.2 * row.get(c, 0), units="points")
        ax_c.text(offsets[c] + chrom_len[c] / 2, 0, c.replace("chr", ""), transform=tr, ha="center", va="top",
                  fontsize=6)

    # one compact horizontal legend below the call track (only categories present)
    present = set(seg_in.call)
    handles = [Patch(facecolor=fc, edgecolor="#8c8c8c" if h else "none", hatch=h, lw=0.3 if h else 0, label=lab)
               for key, lab, fc, h in CALLS if key in present]
    if handles:
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=((m["left"] + pw / 2) / W, 0.4 / H),
                   ncol=len(handles), handlelength=1.4, handleheight=0.9, handletextpad=0.4, columnspacing=1.1,
                   borderaxespad=0, borderpad=0)

    paths = {"png": f"{out_prefix}.png", "pdf": f"{out_prefix}.pdf"}
    fig.savefig(paths["png"], dpi=600)
    fig.savefig(paths["pdf"])
    plt.close(fig)
    plt.style.use("default")

    pd.testing.assert_frame_equal(seg_in, segments)          # inputs untouched
    return pd.DataFrame(drawn), paths
