#!/bin/bash
# Analysepfad für eine Probe (Illumina paired-end, Hybrid-Capture-Panel):
#   fastp -> minimap2 (-x sr, hg19 primär) -> fixmate/sort/markdup
#   -> Target-Coverage (mosdepth) -> Allelzählung an 1000G-SNPs (bcftools)
#   -> BAF-Likelihoods + Off-Target-logR -> Segmentierung, Purity, Calls
#   -> Gen-Coverage -> Plots + Summary
# Wird von run_pipeline.sh aufgerufen, nicht direkt.
#
# Jeder Schritt prüft zuerst, ob sein Ergebnis bereits vorhanden ist, und
# überspringt ihn dann (sicher für die requeue-Partition). Plots und
# Summary werden immer neu erzeugt.
#
# Zwischendateien (BAM, Allelzählungen, Bins) liegen in tmp/<sample>/;
# output/<sample>/ enthält die Deliverables.
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

SAMPLE="$1"; R1="$2"; R2="$3"
THREADS=${SLURM_CPUS_PER_TASK:-$THREADS_DEFAULT}
TMP="tmp/$SAMPLE"
OUT="output/$SAMPLE"
mkdir -p "$TMP" "$OUT/qc" log

export PATH="$ENV_BAF/bin:$PATH"
if command -v module >/dev/null 2>&1; then
    module purge
    ml $MINIMAP2_MODULES $SAMTOOLS_MODULES
fi
PY="$ENV_BAF/bin/python"
COMMON=(--sample "$SAMPLE" --tmp "$TMP" --out "$OUT" --ref "$REF_DIR" --bed "$PANEL_BED")

echo "=== Panel-BAF-Pipeline: $SAMPLE ==="
echo "Threads: $THREADS"
echo "Reads:   $R1 $R2"
echo "Panel:   $PANEL_BED"

# 1. Adapter-/PolyG-Trimming
T1="$TMP/$SAMPLE.R1.trim.fastq.gz"; T2="$TMP/$SAMPLE.R2.trim.fastq.gz"
if [ -s "$T1" ] && [ -s "$T2" ] && [ -s "$OUT/qc/$SAMPLE.fastp.json" ]; then
    echo "--- fastp: bereits vorhanden, überspringe ---"
else
    echo "--- fastp ---"
    fastp -i "$R1" -I "$R2" -o "$T1.tmp.fastq.gz" -O "$T2.tmp.fastq.gz" \
        --detect_adapter_for_pe --trim_poly_g -l 30 -w "$(( THREADS < 16 ? THREADS : 16 ))" -z 3 \
        -h "$OUT/qc/$SAMPLE.fastp.html" -j "$OUT/qc/$SAMPLE.fastp.json" 2> "$TMP/fastp.log"
    mv "$T1.tmp.fastq.gz" "$T1"; mv "$T2.tmp.fastq.gz" "$T2"
fi

# 2. Alignment (minimap2 -x sr) -> fixmate -> sort -> markdup
BAM="$TMP/$SAMPLE.bam"
if [ -s "$BAM" ] && [ -s "$BAM.bai" ] && [ -s "$TMP/$SAMPLE.flagstat.txt" ]; then
    echo "--- Alignment: $BAM bereits vorhanden, überspringe ---"
else
    echo "--- Alignment ---"
    minimap2 -ax sr -t "$THREADS" -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA\tLB:$SAMPLE" \
        "$HG19_MMI_SR" "$T1" "$T2" 2> "$TMP/minimap2.log" \
      | samtools fixmate -m -@ 4 - - \
      | samtools sort -@ "$(( THREADS / 2 ))" -m 1500M -T "$TMP/sort_tmp" -o "$TMP/$SAMPLE.sorted.bam" -
    samtools markdup -@ "$THREADS" -s -f "$TMP/$SAMPLE.markdup.txt" "$TMP/$SAMPLE.sorted.bam" "$BAM.tmp" 2> "$TMP/markdup.log"
    mv "$BAM.tmp" "$BAM"
    samtools index -@ "$THREADS" "$BAM"
    rm -f "$TMP/$SAMPLE.sorted.bam"
    samtools flagstat -@ "$THREADS" "$BAM" > "$TMP/$SAMPLE.flagstat.txt"
fi

# 3. Target-Coverage (dedupliziert, MAPQ >= 20)
if [ -s "$TMP/$SAMPLE.target.regions.bed.gz" ]; then
    echo "--- mosdepth: bereits vorhanden, überspringe ---"
else
    echo "--- mosdepth (Targets) ---"
    mosdepth -t 4 -Q 20 -n --by "$PANEL_BED" "$TMP/$SAMPLE.target" "$BAM"
fi

# 4. Allelzählung an allen SNP-Positionen des Panels (on- und off-target)
AD="$TMP/$SAMPLE.mpileup_AD.tsv.gz"
if [ -s "$AD" ]; then
    echo "--- Allelzählung: bereits vorhanden, überspringe ---"
else
    echo "--- Allelzählung (bcftools mpileup, MAPQ>=20, BQ>=20, ohne Duplikate) ---"
    mkdir -p "$TMP/mpileup"
    count_chr() {
        c=$1
        tabix "$SNP_PANEL" "$c" | cut -f1,2 > "$TMP/mpileup/$c.sites"
        bcftools mpileup -f "$HG19_PRIMARY" -r "$c" -T "$TMP/mpileup/$c.sites" \
            -q 20 -Q 20 -d 100000 -I -a FORMAT/AD,FORMAT/DP -Ou "$BAM" 2>/dev/null \
          | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' > "$TMP/mpileup/$c.ad.tsv"
    }
    export -f count_chr; export TMP SNP_PANEL HG19_PRIMARY BAM
    for c in $(seq 1 22) X; do echo chr$c; done | xargs -P "$THREADS" -I{} bash -c 'count_chr {}'
    for c in $(seq 1 22) X; do cat "$TMP/mpileup/chr$c.ad.tsv"; done | gzip > "$AD.tmp"
    mv "$AD.tmp" "$AD"
    rm -rf "$TMP/mpileup"
fi

# 5. Fragment-Mittelpunkte (für Off-Target-Read-Depth)
if [ -s "$TMP/$SAMPLE.fragment_midpoints.npz" ]; then
    echo "--- Fragmentpositionen: bereits vorhanden, überspringe ---"
else
    echo "--- Fragmentpositionen ---"
    "$PY" scripts/read_positions.py "${COMMON[@]}"
fi

# 6. SNP-Tabelle (Allelzählungen + Populations-AF + Target-Abstand)
if [ -s "$TMP/$SAMPLE.snps.tsv.gz" ]; then
    echo "--- SNP-Tabelle: bereits vorhanden, überspringe ---"
else
    echo "--- SNP-Tabelle ---"
    "$PY" scripts/snp_table.py "${COMMON[@]}"
fi

# 7. BAF-Likelihoods, logR, Segmentierung
if [ -s "$TMP/$SAMPLE.segments.tsv" ] && [ -s "$TMP/$SAMPLE.params.json" ]; then
    echo "--- BAF/logR-Profil: bereits vorhanden, überspringe ---"
else
    echo "--- BAF/logR-Profil ---"
    "$PY" scripts/baf_cn_profile.py "${COMMON[@]}"
fi

# 8. Purity-Fit + Calls, 9. Gen-Coverage, 10. Plots + Summary (immer neu)
echo "--- Calls, Gen-Coverage, Plots ---"
"$PY" scripts/call_segments.py "${COMMON[@]}"
"$PY" scripts/gene_coverage.py "${COMMON[@]}"
"$PY" scripts/plot_profile.py "${COMMON[@]}" --focal "$FOCAL_GENES"
"$PY" scripts/summarize.py "${COMMON[@]}"

echo "=== Fertig: $SAMPLE -> $OUT ==="
