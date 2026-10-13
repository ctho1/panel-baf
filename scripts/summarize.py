#!/usr/bin/env python
"""QC + result summary. Output: output/<sample>/<sample>.summary.json"""
import json
import os
import subprocess

import pandas as pd

from baf_common import parse_args

P, _ = parse_args(__doc__)
fastp = json.load(open(f"{P.out}/qc/{P.sample}.fastp.json"))
params = json.load(open(P.t("params.json")))
purity = json.load(open(P.o("purity.json")))
segs = pd.read_csv(P.o("segments.tsv"), sep="\t")
genes = pd.read_csv(P.o("gene_coverage.tsv"), sep="\t")

markdup = {}
for line in open(P.t("markdup.txt")):
    if ":" in line:
        k, v = line.split(":", 1)
        markdup[k.strip()] = v.strip()

# usable reads (non-dup, primary, MAPQ >= 20) and on-target fraction (samtools from module / conda)
threads = os.environ.get("SLURM_CPUS_PER_TASK", "4")


def count(*extra):
    cmd = ["samtools", "view", "-c", "-@", threads, "-F", "0xD04", "-q", "20", *extra, P.t("bam")]
    return int(subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip())


tot, on = count(), count("-L", P.bed)
reg = pd.read_csv(P.t("target.regions.bed.gz"), sep="\t", header=None)
depth = reg.iloc[:, -1]

events = segs[~segs.call.str.startswith(("Neutral", "Not assessed", "Pericentromeric"))]
summary = dict(
    sample=P.sample,
    reads_raw_pairs=fastp["summary"]["before_filtering"]["total_reads"] // 2,
    reads_after_fastp=fastp["summary"]["after_filtering"]["total_reads"],
    duplicates=int(markdup.get("DUPLICATE TOTAL", 0)),
    reads_usable_dedup_mapq20=tot,
    on_target_fraction=round(on / tot, 4) if tot else None,
    on_target_mean_depth=round(float((depth * (reg[2] - reg[1])).sum() / (reg[2] - reg[1]).sum()), 1),
    on_target_median_depth=round(float(depth.median()), 1),
    sex=params["sex"], population_af=params["population"], het_rate_scale=params["het_rate_scale"],
    n_snps_covered=params["n_snps_covered"], n_informative_snps=params["n_inf_snps"],
    sigma_logR_1Mb=params["sigma_logR_bin"],
    purity=purity["purity"], purity_CI=purity["purity_CI_2LL"],
    events=[dict(arm=r.arm, start=int(r.start), end=int(r.end), call=r.call, confidence=r.confidence,
                 logR=round(r.logR, 3), b=r.b, b_CI=[r.b_lo, r.b_hi], state=r.fit_state) for _, r in events.iterrows()],
    focal_gene_outliers=genes[(genes.log2_rel.abs() > 0.6) & (genes.z.abs() > 4)][
        ["gene", "log2_rel", "n_het_snps", "median_minor_af"]].round(3).to_dict("records"),
)
with open(P.o("summary.json"), "w") as fh:
    json.dump(summary, fh, indent=1, default=float)
print(json.dumps({k: v for k, v in summary.items() if k not in ("events", "focal_gene_outliers")}, indent=1))
