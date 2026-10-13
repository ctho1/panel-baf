#!/bin/bash
# Stellt alle Referenzdaten bereit (einmalig, vor allen Probenjobs):
#   - hg19 nur Primärchromosomen + minimap2-Index (-x sr)
#   - SNP-Panel: 1000G Phase 3, biallelische SNVs, globale MAF >= SNP_MIN_MAF,
#     mit Superpopulations-AF (AF, EAS, AMR, AFR, EUR, SAS)
#   - cytoBand, refGene, Alu-Elemente (RepeatMasker), GC-/Mappability-Wigs (10 kb)
# Bereits vorhandene Dateien werden übersprungen.
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

THREADS=${SLURM_CPUS_PER_TASK:-$THREADS_DEFAULT}
mkdir -p "$REF_DIR" log
export PATH="$ENV_BAF/bin:$PATH"
if command -v module >/dev/null 2>&1; then
    module purge
    ml $MINIMAP2_MODULES $SAMTOOLS_MODULES
fi

fetch() {  # fetch <url> <dest>  (atomic)
    [ -s "$2" ] && return 0
    echo "Lade $1 ..."
    curl -fL --retry 3 -o "$2.part" "$1"
    mv "$2.part" "$2"
}

echo "=== Referenz vorbereiten (hg19) ==="

# 1. hg19 Primärchromosomen
if [ ! -s "$HG19_PRIMARY.fai" ]; then
    SRC="$EXISTING_HG19_FASTA"
    if [ ! -s "$SRC" ]; then
        SRC="$REF_DIR/hg19.fa"
        fetch "$HG19_FASTA_URL" "$SRC.gz"
        [ -s "$SRC" ] || gzip -dc "$SRC.gz" > "$SRC"
    fi
    echo "Extrahiere Primärchromosomen aus $SRC"
    if [ ! -s "$SRC.fai" ]; then
        # geteilte Referenz evtl. schreibgeschützt -> lokal verlinken und dort indexieren
        ln -sf "$SRC" "$REF_DIR/hg19_source.fa"
        SRC="$REF_DIR/hg19_source.fa"
        [ -s "$SRC.fai" ] || samtools faidx "$SRC"
    fi
    samtools faidx "$SRC" $(for c in $(seq 1 22) X Y M; do echo chr$c; done) > "$HG19_PRIMARY.tmp"
    mv "$HG19_PRIMARY.tmp" "$HG19_PRIMARY"
    samtools faidx "$HG19_PRIMARY"
fi

# 2. minimap2-Index (Short Reads)
if [ ! -s "$HG19_MMI_SR" ]; then
    echo "Baue minimap2-Index (sr) ..."
    minimap2 -x sr -t "$THREADS" -d "$HG19_MMI_SR.tmp" "$HG19_PRIMARY"
    mv "$HG19_MMI_SR.tmp" "$HG19_MMI_SR"
fi

# 3. SNP-Panel aus 1000G Phase 3
if [ ! -s "$SNP_PANEL.tbi" ]; then
    KG="$REF_DIR/1000G_phase3_v5c.sites.vcf.gz"
    fetch "$KG_SITES_URL" "$KG"
    fetch "$KG_SITES_URL.tbi" "$KG.tbi"
    echo "Filtere biallelische SNVs (MAF >= $SNP_MIN_MAF) ..."
    bcftools view -v snps -m2 -M2 -i 'INFO/VT="SNP"' -t ^Y,MT "$KG" -Ou \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/EAS_AF\t%INFO/AMR_AF\t%INFO/AFR_AF\t%INFO/EUR_AF\t%INFO/SAS_AF\n' \
      | awk -F'\t' -v m="$SNP_MIN_MAF" 'BEGIN{OFS="\t"} { af=$5+0; maf=(af<0.5)?af:1-af; if (maf>=m) { $1="chr"$1; print } }' \
      | bgzip -@ "$THREADS" > "$SNP_PANEL.tmp"
    mv "$SNP_PANEL.tmp" "$SNP_PANEL"
    tabix -f -s1 -b2 -e2 "$SNP_PANEL"
    echo "SNPs im Panel: $(bgzip -dc "$SNP_PANEL" | wc -l)"
    [ "${KEEP_DOWNLOADS:-false}" = "true" ] || rm -f "$KG" "$KG.tbi"
fi

# 4. Annotation
fetch "$UCSC_HG19_DB_URL/cytoBand.txt.gz" "$REF_DIR/cytoBand_hg19.txt.gz"
fetch "$UCSC_HG19_DB_URL/refGene.txt.gz" "$REF_DIR/refGene_hg19.txt.gz"
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

echo "=== Referenz bereit ==="
