"""Shared paths, arguments and helpers for the panel BAF pipeline scripts."""
import argparse
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
CHROMS = AUTOSOMES + ["chrX"]
BIN = 1_000_000                                   # analysis bin (logR, BAF)
BGRID = np.round(np.arange(0, 0.5001, 0.01), 2)   # grid for b (minor-allele fraction of hets)
REF_FASTA = "hg19.fa"                             # references/hg19.fa (shared PALMA hg19, symlinked)


@dataclass
class Paths:
    sample: str
    tmp: str
    out: str
    ref: str
    bed: str

    def t(self, suffix):          # tmp/<sample>/<sample>.<suffix>
        return os.path.join(self.tmp, f"{self.sample}.{suffix}")

    def o(self, suffix):          # output/<sample>/<sample>.<suffix>
        return os.path.join(self.out, f"{self.sample}.{suffix}")

    def r(self, name):            # references/<name>
        return os.path.join(self.ref, name)


def parse_args(description, extra=None):
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--tmp", required=True, help="per-sample working directory (tmp/<sample>)")
    ap.add_argument("--out", required=True, help="per-sample output directory (output/<sample>)")
    ap.add_argument("--ref", required=True, help="reference directory (references/)")
    ap.add_argument("--bed", required=True, help="panel target BED (hg19)")
    if extra:
        extra(ap)
    a = ap.parse_args()
    os.makedirs(a.tmp, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    return Paths(a.sample, a.tmp, a.out, a.ref, a.bed), a


def read_bed(path):
    bed = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2], names=["chrom", "s", "e"])
    return bed[~bed.chrom.astype(str).str.startswith(("track", "browser"))].astype({"s": int, "e": int})


def load_arms(ref):
    """Centromere (acen) start/end per chromosome from UCSC cytoBand."""
    cb = pd.read_csv(os.path.join(ref, "cytoBand_hg19.txt.gz"), sep="\t", header=None,
                     names=["chrom", "s", "e", "band", "stain"])
    acen = cb[cb.stain == "acen"].groupby("chrom").agg(s=("s", "min"), e=("e", "max"))
    return {c: (int(r.s), int(r.e)) for c, r in acen.iterrows()}


def load_wig(path):
    out, cur = {}, None
    for line in open(path):
        if line.startswith("fixedStep"):
            cur = "chr" + dict(kv.split("=") for kv in line.split()[1:])["chrom"].replace("chr", "")
            out[cur] = []
        else:
            out[cur].append(float(line))
    return {k: np.array(v) for k, v in out.items()}


def bci(ll, level=1.92):
    """MLE and profile-likelihood CI (chi2_1, 95 %) of b on BGRID."""
    i = int(np.argmax(ll))
    ok = ll >= ll[i] - level
    return BGRID[i], BGRID[ok].min(), BGRID[ok].max()


def read_fai(fasta):
    return pd.read_csv(fasta + ".fai", sep="\t", header=None, usecols=[0, 1, 2, 3, 4],
                       names=["chrom", "len", "offset", "linebases", "linewidth"]).set_index("chrom")

