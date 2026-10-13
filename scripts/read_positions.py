#!/usr/bin/env python
"""Fragment midpoints of non-duplicate primary reads per chromosome, used for the off-target read-depth log2
ratio: <chrom> = MAPQ >= 20, <chrom>_all = any MAPQ (to flag windows where reads cannot be placed uniquely,
e.g. regions with alternate haplotype contigs such as the MHC). Proper pairs count once (read 1, midpoint of
the fragment); other reads count with their start position (offset < read length, irrelevant for 10-kb
windows). Reads come from `samtools view` (module or conda), no pysam needed.
Output: tmp/<sample>/<sample>.fragment_midpoints.npz"""
import io
import subprocess
import sys

import numpy as np
import pandas as pd

from baf_common import CHROMS, parse_args

P, _ = parse_args(__doc__)
pos = {}
for c in CHROMS + ["chrY"]:
    cmd = f"samtools view -F 0xD04 {P.t('bam')} {c} | cut -f2,4,5,7,8,9"
    out = subprocess.run(cmd, shell=True, check=True, capture_output=True).stdout
    if not out:
        pos[c] = pos[c + "_all"] = np.zeros(0, np.int32)
        continue
    d = pd.read_csv(io.BytesIO(out), sep="\t", header=None, names=["flag", "pos", "mapq", "rnext", "pnext", "tlen"],
                    dtype={"flag": np.int32, "pos": np.int64, "mapq": np.int32, "rnext": str, "pnext": np.int64,
                           "tlen": np.int64})
    proper = (d.flag & 0x2) > 0
    read1 = (d.flag & 0x40) > 0
    paired = (d.flag & 0x1) > 0
    mate_same = paired & ((d.flag & 0x8) == 0) & (d.rnext == "=")
    mid_pair = np.minimum(d.pos, d.pnext) - 1 + d.tlen.abs() // 2
    keep = (proper & read1) | (~proper & mate_same & read1) | (~proper & ~mate_same)
    mid = np.where(proper, mid_pair, d.pos - 1)
    q20 = (d.mapq >= 20).values
    pos[c] = np.sort(mid[keep.values & q20].astype(np.int32))
    pos[c + "_all"] = np.sort(mid[keep.values].astype(np.int32))
    print(c, len(pos[c]), file=sys.stderr)
np.savez_compressed(P.t("fragment_midpoints.npz"), **pos)
