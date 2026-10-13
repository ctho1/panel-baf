#!/usr/bin/env python3
"""Scans input/ (flat, no recursion) for Illumina R1/R2 pairs, one line per sample:
    "<sample>\\t<R1>\\t<R2>"
Pairs are matched on the text around the "_R1_001"/"_R2_001" marker, e.g.
07_SIL-3036_N3615_26_S7_R1_001.fastq.gz -> sample "07_SIL-3036_N3615_26_S7".
"""
import re
import sys
from pathlib import Path

ILLUMINA_RE = re.compile(r"^(.*)_R([12])_001(.*)$")


def strip_ext(name):
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return None


def discover(input_dir):
    r1, r2, unpaired = {}, {}, []
    for p in sorted(Path(input_dir).iterdir()):
        stem = strip_ext(p.name) if p.is_file() else None
        if stem is None:
            continue
        m = ILLUMINA_RE.match(stem)
        if not m:
            unpaired.append(p)
            continue
        prefix, read, suffix = m.groups()
        (r1 if read == "1" else r2)[prefix + suffix] = p
    pairs = [(s, r1[s], r2[s]) for s in sorted(r1) if s in r2]
    skipped = [r1[s] for s in r1 if s not in r2] + [r2[s] for s in r2 if s not in r1] + unpaired
    return pairs, skipped


def main():
    input_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "input")
    if not input_dir.is_dir():
        return
    pairs, skipped = discover(input_dir)
    for sample, a, b in pairs:
        print(f"{sample}\t{a}\t{b}")
    for p in skipped:
        print(f"SKIP {p} (no R1/R2 pair)", file=sys.stderr)


if __name__ == "__main__":
    main()
