#!/usr/bin/env python
"""Merge mpileup allele depths with the SNP panel (alleles + 1000G superpopulation AFs) and annotate the
distance to the nearest panel target (on <= 150 bp, near <= 1 kb, else off).
Output: tmp/<sample>/<sample>.snps.tsv.gz"""
import sys

import numpy as np
import pandas as pd

from baf_common import CHROMS, parse_args, read_bed

P, _ = parse_args(__doc__)
ad = pd.read_csv(P.t("mpileup_AD.tsv.gz"), sep="\t", header=None, names=["chrom", "pos", "ref_obs", "alts_obs", "ad"],
                 dtype={"chrom": str, "pos": np.int64})
panel = pd.read_csv(P.r("snps_1kg_maf05.tsv.gz"), sep="\t", header=None,
                    names=["chrom", "pos", "ref", "alt", "AF", "EAS", "AMR", "AFR", "EUR", "SAS"],
                    dtype={"chrom": str, "pos": np.int64, "ref": str, "alt": str}, na_values=["."])
panel = panel.drop_duplicates(["chrom", "pos"], keep=False)        # multi-record positions are ambiguous
df = ad.merge(panel, on=["chrom", "pos"], how="inner")
df = df[(df.ref_obs == df.ref) & df.chrom.isin(CHROMS) & df.AF.notna()].copy()
for pop in ["EAS", "AMR", "AFR", "EUR", "SAS"]:
    df[pop] = df[pop].fillna(df.AF)

ref_n, alt_n = [], []
for alts, adv, alt in zip(df.alts_obs.values, df.ad.values, df.alt.values):
    cnt = [int(x) for x in adv.split(",")]
    al = alts.split(",")
    ref_n.append(cnt[0])
    alt_n.append(cnt[1 + al.index(alt)] if alt in al else 0)
df["ref_n"], df["alt_n"] = ref_n, alt_n
df["depth"] = df.ref_n + df.alt_n
df = df[df.depth > 0]

bed = read_bed(P.bed)
dist = np.full(len(df), np.inf)
for c, sub in df.groupby("chrom"):
    b = bed[bed.chrom == c].sort_values("s")
    if b.empty:
        continue
    s, e, p = b.s.values, b.e.values, sub.pos.values - 1
    i = np.searchsorted(s, p, side="right") - 1
    dl = np.where(i >= 0, np.maximum(p - e[np.clip(i, 0, None)] + 1, 0), np.inf)
    j = i + 1
    dr = np.where(j < len(s), s[np.clip(j, None, len(s) - 1)] - p, np.inf)
    dist[df.index.get_indexer(sub.index)] = np.minimum(dl, dr)
df["target_dist"] = dist
df["region"] = np.select([dist <= 150, dist <= 1000], ["on", "near"], "off")
df["chrom"] = pd.Categorical(df.chrom, CHROMS, ordered=True)
df = df.sort_values(["chrom", "pos"])
cols = ["chrom", "pos", "ref", "alt", "AF", "EAS", "AMR", "AFR", "EUR", "SAS", "ref_n", "alt_n", "depth",
        "target_dist", "region"]
df[cols].to_csv(P.t("snps.tsv.gz"), sep="\t", index=False)
print(df.groupby("region", observed=True).depth.describe(), file=sys.stderr)
