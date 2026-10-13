#!/bin/bash
# Provides all reference data (run once, before the sample jobs):
#   - hg19 + aligner index: the reference shared with mNGS-pipeline (EXISTING_HG19_FASTA)
#     is symlinked; otherwise UCSC download (index is only built for ALIGNER=bwa)
#   - SNP panel: 1000 Genomes phase 3 biallelic SNVs with global MAF >= SNP_MIN_MAF,
#     incl. superpopulation AFs (AF, EAS, AMR, AFR, EUR, SAS)
#   - cytoBand, Alu elements (RepeatMasker), GC / mappability wigs (10 kb)
# Existing files are skipped.
#SBATCH --job-name=prepare_reference
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --error=./log/%x_%j.err.txt
#SBATCH --output=./log/%x_%j.out.txt

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
source scripts/config.sh
source scripts/software.sh

THREADS=${SLURM_CPUS_PER_TASK:-$THREADS_DEFAULT}
mkdir -p "$REF_DIR" log
load_software

fetch() {  # fetch <url> <dest>  (atomic)
    [ -s "$2" ] && return 0
    echo "Downloading $1 ..."
    curl -fL --retry 3 -o "$2.part" "$1"
    mv "$2.part" "$2"
}

echo "=== Preparing reference (hg19, $ALIGNER) ==="

# 1. hg19 (shared with mNGS-pipeline) + .fai
if [ ! -s "$HG19_FASTA" ]; then
    if [ -s "$EXISTING_HG19_FASTA" ]; then
        echo "Using existing hg19 FASTA: $EXISTING_HG19_FASTA"
        ln -sf "$EXISTING_HG19_FASTA" "$HG19_FASTA"
    else
        fetch "$HG19_FASTA_URL" "$HG19_FASTA.gz"
        gzip -dc "$HG19_FASTA.gz" > "$HG19_FASTA.tmp" && mv "$HG19_FASTA.tmp" "$HG19_FASTA"
    fi
fi
if [ ! -s "$HG19_FASTA.fai" ]; then
    if [ -s "$EXISTING_HG19_FASTA.fai" ]; then
        ln -sf "$EXISTING_HG19_FASTA.fai" "$HG19_FASTA.fai"
    else
        samtools faidx "$HG19_FASTA"
    fi
fi

# 2. Aligner index: link an existing index, otherwise build it
EXTS=$(aligner_index_exts)
have_index() { for e in $EXTS; do [ -s "$1.$e" ] || return 1; done; }
if ! have_index "$HG19_FASTA"; then
    if have_index "$EXISTING_HG19_FASTA"; then
        echo "Using existing $ALIGNER index: $EXISTING_HG19_FASTA.*"
        for e in $EXTS; do ln -sf "$EXISTING_HG19_FASTA.$e" "$HG19_FASTA.$e"; done
    elif [ "$ALIGNER" = "bwa" ]; then
        echo "Building bwa index (~1 h, ~5 GB RAM) ..."
        bwa index -a bwtsw "$HG19_FASTA"
    else
        echo "No bwa-mem2 index found ($EXISTING_HG19_FASTA.bwt.2bit.64)." >&2
        echo "Building one needs ~90 GB RAM: point EXISTING_HG19_FASTA to an indexed reference," >&2
        echo "build it manually (bwa-mem2 index) or use ALIGNER=bwa." >&2
        exit 1
    fi
fi

# 3. SNP panel from 1000 Genomes phase 3
if [ ! -s "$SNP_PANEL.tbi" ]; then
    KG="$REF_DIR/1000G_phase3_v5c.sites.vcf.gz"
    fetch "$KG_SITES_URL" "$KG"
    fetch "$KG_SITES_URL.tbi" "$KG.tbi"
    echo "Filtering biallelic SNVs (MAF >= $SNP_MIN_MAF) ..."
    bcftools view -v snps -m2 -M2 -i 'INFO/VT="SNP"' -t ^Y,MT "$KG" -Ou \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/EAS_AF\t%INFO/AMR_AF\t%INFO/AFR_AF\t%INFO/EUR_AF\t%INFO/SAS_AF\n' \
      | awk -F'\t' -v m="$SNP_MIN_MAF" 'BEGIN{OFS="\t"} { af=$5+0; maf=(af<0.5)?af:1-af; if (maf>=m) { $1="chr"$1; print } }' \
      | bgzip -@ "$THREADS" > "$SNP_PANEL.tmp"
    mv "$SNP_PANEL.tmp" "$SNP_PANEL"
    tabix -f -s1 -b2 -e2 "$SNP_PANEL"
    echo "SNPs in panel: $(bgzip -dc "$SNP_PANEL" | wc -l)"
    [ "${KEEP_DOWNLOADS:-false}" = "true" ] || rm -f "$KG" "$KG.tbi"
fi

# 4. Annotation
fetch "$UCSC_HG19_DB_URL/cytoBand.txt.gz" "$REF_DIR/cytoBand_hg19.txt.gz"
if [ ! -s "$REF_DIR/alu_hg19.bed.gz" ]; then
    fetch "$UCSC_HG19_DB_URL/rmsk.txt.gz" "$REF_DIR/rmsk_hg19.txt.gz"
    gzip -dc "$REF_DIR/rmsk_hg19.txt.gz" \
      | awk -F'\t' 'BEGIN{OFS="\t"} $6 !~ /_/ && $13 == "Alu" {print $6, $7, $8}' \
      | gzip > "$REF_DIR/alu_hg19.bed.gz.tmp"
    mv "$REF_DIR/alu_hg19.bed.gz.tmp" "$REF_DIR/alu_hg19.bed.gz"
    rm -f "$REF_DIR/rmsk_hg19.txt.gz"
fi
for w in gc_hg19_10kb.wig map_hg19_10kb.wig; do
    fetch "$ICHORCNA_EXTDATA_URL/$w" "$REF_DIR/$w"
done

echo "=== Reference ready ==="
