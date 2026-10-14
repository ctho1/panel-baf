#!/bin/bash
# Per-sample job (Illumina paired-end, hybrid-capture panel). Called by run_pipeline.sh.
#   [fastp, optional] -> bwa-mem2 mem (hg19) -> fixmate/sort/markdup
#   -> allele counts at 1000G SNPs (bcftools) -> BAF likelihoods + off-target log2 ratio
#   -> segmentation, purity fit, calls -> genome-wide figure (PNG + PDF)
#
# Every step is skipped if its result already exists (safe on the preemptible requeue
# partition); calls and the figure are always regenerated.
#   tmp/<sample>/     BAM, allele counts, bins, segments, purity (intermediate files)
#   output/<sample>/  <sample>.genome_BAF_logR.png and .pdf (the only deliverables)
# Trimmed FASTQs are deleted at the end.
#
# Usage: sample_pipeline.sh <sample> <R1.fastq.gz> <R2.fastq.gz>
#SBATCH --job-name=baf
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --error=./log/%x_%j.err.txt
#SBATCH --output=./log/%x_%j.out.txt

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
source scripts/config.sh
source scripts/software.sh

SAMPLE="$1"; R1="$2"; R2="$3"
THREADS=${SLURM_CPUS_PER_TASK:-$THREADS_DEFAULT}
TMP="tmp/$SAMPLE"
OUT="output/$SAMPLE"
mkdir -p "$TMP" "$OUT" log
COMMON=(--sample "$SAMPLE" --tmp "$TMP" --out "$OUT" --ref "$REF_DIR" --bed "$PANEL_BED")
BAM="$TMP/$SAMPLE.bam"
T1="$TMP/$SAMPLE.R1.trim.fastq.gz"; T2="$TMP/$SAMPLE.R2.trim.fastq.gz"
bam_done() { [ -s "$BAM" ] && [ -s "$BAM.bai" ] && [ -s "$TMP/$SAMPLE.flagstat.txt" ]; }
die_with_log() {  # die_with_log <step> <log>: copy the end of a tool log into the SLURM .err and abort
    echo "ERROR: $1 failed -- last lines of $2:" >&2
    tail -n 25 "$2" >&2
    exit 1
}

echo "=== panel-baf: $SAMPLE ==="
echo "Threads: $THREADS | aligner: $ALIGNER | reference: $HG19_FASTA"
echo "Reads:   $R1 $R2"
echo "Panel:   $PANEL_BED"

# 1. Optional adapter/polyG trimming (only if fastp is available, e.g. via FASTP_MODULES)
READS=("$R1" "$R2")
use_modules $FASTP_MODULES
if bam_done; then
    echo "--- trimming: alignment already present, skipping ---"
elif command -v fastp >/dev/null 2>&1; then
    if [ -s "$T1" ] && [ -s "$T2" ]; then
        echo "--- fastp: already present, skipping ---"
    else
        echo "--- fastp ---"
        fastp -i "$R1" -I "$R2" -o "$T1.tmp.fastq.gz" -O "$T2.tmp.fastq.gz" \
            --detect_adapter_for_pe --trim_poly_g -l 30 -w "$(( THREADS < 16 ? THREADS : 16 ))" -z 3 \
            -h "$TMP/$SAMPLE.fastp.html" -j "$TMP/$SAMPLE.fastp.json" 2> "$TMP/fastp.log" \
            || die_with_log fastp "$TMP/fastp.log"
        mv "$T1.tmp.fastq.gz" "$T1"; mv "$T2.tmp.fastq.gz" "$T2"
    fi
    READS=("$T1" "$T2")
else
    echo "--- fastp not available: no trimming (adapters are soft-clipped by the aligner) ---"
fi

# 2. Alignment (bwa-mem2 mem / bwa mem) -> fixmate -> sort -> markdup
if bam_done; then
    echo "--- alignment: $BAM already present, skipping ---"
else
    echo "--- alignment ($ALIGNER mem) ---"
    use_modules $BWA_MODULES $SAMTOOLS_MODULES
    "$ALIGNER" mem -t "$THREADS" -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA\tLB:$SAMPLE" \
        "$HG19_FASTA" "${READS[@]}" 2> "$TMP/$ALIGNER.log" \
      | samtools fixmate -m -@ 4 - - \
      | samtools sort -@ "$(( THREADS / 2 ))" -m 1500M -T "$TMP/sort_tmp" -o "$TMP/$SAMPLE.sorted.bam" - \
      || die_with_log "$ALIGNER mem" "$TMP/$ALIGNER.log"
    samtools markdup -@ "$THREADS" -s -f "$TMP/$SAMPLE.markdup.txt" "$TMP/$SAMPLE.sorted.bam" "$BAM.tmp" 2> "$TMP/markdup.log" \
      || die_with_log "samtools markdup" "$TMP/markdup.log"
    mv "$BAM.tmp" "$BAM"
    samtools index -@ "$THREADS" "$BAM"
    rm -f "$TMP/$SAMPLE.sorted.bam"
    samtools flagstat -@ "$THREADS" "$BAM" > "$TMP/$SAMPLE.flagstat.txt"
fi

# 3. Allele counts at all SNP-panel positions (on- and off-target)
AD="$TMP/$SAMPLE.mpileup_AD.tsv.gz"
if [ -s "$AD" ]; then
    echo "--- allele counts: already present, skipping ---"
else
    echo "--- allele counts (bcftools mpileup, MAPQ>=20, BQ>=20, no duplicates) ---"
    use_modules $BCFTOOLS_MODULES
    mkdir -p "$TMP/mpileup"
    count_chr() {
        c=$1
        tabix "$SNP_PANEL" "$c" | cut -f1,2 > "$TMP/mpileup/$c.sites"
        bcftools mpileup -f "$HG19_FASTA" -r "$c" -T "$TMP/mpileup/$c.sites" \
            -q 20 -Q 20 -d 100000 -I -a FORMAT/AD,FORMAT/DP -Ou "$BAM" 2>/dev/null \
          | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' > "$TMP/mpileup/$c.ad.tsv"
    }
    export -f count_chr; export TMP SNP_PANEL HG19_FASTA BAM
    for c in $(seq 1 22) X; do echo chr$c; done | xargs -P "$THREADS" -I{} bash -c 'count_chr {}'
    for c in $(seq 1 22) X; do cat "$TMP/mpileup/chr$c.ad.tsv"; done | gzip > "$AD.tmp"
    mv "$AD.tmp" "$AD"
    rm -rf "$TMP/mpileup"
fi

# 4.-8. Python steps (read_positions.py also calls samtools)
use_modules $PYTHON_MODULES $SAMTOOLS_MODULES

# 4. Fragment midpoints (off-target read depth)
if [ -s "$TMP/$SAMPLE.fragment_midpoints.npz" ]; then
    echo "--- fragment positions: already present, skipping ---"
else
    echo "--- fragment positions ---"
    python3 scripts/read_positions.py "${COMMON[@]}"
fi

# 5. SNP table (allele counts + population AFs + distance to targets)
if [ -s "$TMP/$SAMPLE.snps.tsv.gz" ]; then
    echo "--- SNP table: already present, skipping ---"
else
    echo "--- SNP table ---"
    python3 scripts/snp_table.py "${COMMON[@]}"
fi

# 6. BAF likelihoods, log2 ratio, segmentation
if [ -s "$TMP/$SAMPLE.segments.tsv" ] && [ -s "$TMP/$SAMPLE.params.json" ]; then
    echo "--- BAF/log2 profile: already present, skipping ---"
else
    echo "--- BAF/log2 profile ---"
    python3 scripts/baf_cn_profile.py "${COMMON[@]}"
fi

# 7. Purity fit + calls, 8. figure (always regenerated)
echo "--- calls and figure ---"
python3 scripts/call_segments.py "${COMMON[@]}"
python3 scripts/plot_profile.py "${COMMON[@]}" --width-mm "$FIG_WIDTH_MM"

# 9. Clean up: trimmed FASTQs are no longer needed (the BAM stays in tmp/)
rm -f "$T1" "$T2"

echo "=== done: $SAMPLE -> $OUT/$SAMPLE.genome_BAF_logR.{png,pdf} ==="
