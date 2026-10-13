#!/usr/bin/env python
"""Fragment midpoints of usable reads (one per pair, non-duplicate, primary, MAPQ >= 20) per chromosome,
used for the off-target read-depth log2 ratio. Output: tmp/<sample>/<sample>.fragment_midpoints.npz"""
import sys

import numpy as np
import pysam

from baf_common import CHROMS, parse_args

P, _ = parse_args(__doc__)
bam = pysam.AlignmentFile(P.t("bam"))
pos = {}
for c in CHROMS + ["chrY"]:
    mids = []
    for r in bam.fetch(c):
        if r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail or r.mapping_quality < 20:
            continue
        if r.is_paired and r.is_proper_pair:
            if r.is_read1:
                start = min(r.reference_start, r.next_reference_start)
                mids.append(start + abs(r.template_length) // 2)
        elif r.is_paired and not r.mate_is_unmapped and r.next_reference_name == c:
            if r.is_read1:
                mids.append(r.reference_start + r.query_alignment_length // 2)
        else:
            mids.append(r.reference_start + r.query_alignment_length // 2)
    pos[c] = np.sort(np.asarray(mids, dtype=np.int32))
    print(c, len(pos[c]), file=sys.stderr)
np.savez_compressed(P.t("fragment_midpoints.npz"), **pos)
