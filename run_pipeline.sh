#!/bin/bash
# Entry point: scans input/ for Illumina R1/R2 pairs and creates the genome-wide BAF / copy-number
# figure for each sample. Submits one SLURM job per sample if sbatch is available, otherwise runs
# the samples serially (local testing). Run from the repository root.
#
# Usage:
#   bash run_pipeline.sh [--mail user@uni-muenster.de]
#   bash run_pipeline.sh --partition normal   # skip the automatic partition choice
#   PANEL_BED=/path/panel.bed bash run_pipeline.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"
source scripts/config.sh
source scripts/software.sh

PARTITION_OVERRIDE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mail) MAIL_USER="$2"; shift 2 ;;
        --partition) PARTITION_OVERRIDE="$2"; shift 2 ;;
        *)      echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [ ! -s "$PANEL_BED" ]; then
    echo "Panel BED not found: $PANEL_BED" >&2
    echo "Copy the target BED (hg19) there or set PANEL_BED=/path/to/panel.bed" >&2
    exit 1
fi
if ! command -v module >/dev/null 2>&1 && [ ! -x "$ENV_BAF/bin/python" ]; then
    echo "Neither a module system nor the conda env ($ENV_BAF) -- run 'bash scripts/setup_envs.sh' first (local use)." >&2
    exit 1
fi

MAIL_ARGS=()
if [[ -n "$MAIL_USER" ]]; then
    MAIL_ARGS=(--mail-user="$MAIL_USER" --mail-type=ALL)
fi

mkdir -p log output tmp input references

USE_SLURM=0
command -v sbatch >/dev/null 2>&1 && USE_SLURM=1
if [ "$USE_SLURM" = "1" ]; then
    echo "SLURM detected -- submitting jobs with sbatch."
else
    echo "No SLURM (sbatch not found) -- running samples serially."
fi

SBATCH_PARTITION_ARGS=()
if [ "$USE_SLURM" = "1" ]; then
    if [ -n "$PARTITION_OVERRIDE" ]; then
        SELECTED_PARTITIONS="$PARTITION_OVERRIDE"
        echo "SLURM partition set manually: $SELECTED_PARTITIONS"
    elif [ "$SLURM_AUTO_PARTITION" = "true" ]; then
        SELECTED_PARTITIONS=$(python3 scripts/select_slurm_partitions.py \
            --partitions "$SLURM_CPU_PARTITIONS" \
            --fallback "$SLURM_FALLBACK_PARTITION")
        echo "SLURM partitions ranked by idle CPUs: $SELECTED_PARTITIONS"
    else
        SELECTED_PARTITIONS="$SLURM_FALLBACK_PARTITION"
        echo "Automatic partition choice disabled: $SELECTED_PARTITIONS"
    fi
    SBATCH_PARTITION_ARGS=(--partition="$SELECTED_PARTITIONS")
fi

# Reference (prepared once; all sample jobs depend on it)
REF_DEP=()
REF_FILES=("$HG19_FASTA.fai" "$SNP_PANEL.tbi" "$REF_DIR/alu_hg19.bed.gz" "$REF_DIR/cytoBand_hg19.txt.gz"
           "$REF_DIR/gc_hg19_10kb.wig" "$REF_DIR/map_hg19_10kb.wig")
for e in $(aligner_index_exts); do REF_FILES+=("$HG19_FASTA.$e"); done
missing=0
for f in "${REF_FILES[@]}"; do [ -s "$f" ] || missing=1; done
if [ "$missing" = "1" ]; then
    echo ""
    echo "=== Reference incomplete -- preparing it ==="
    if [ "$USE_SLURM" = "1" ]; then
        REF_JOBID=$(sbatch --parsable "${MAIL_ARGS[@]}" "${SBATCH_PARTITION_ARGS[@]}" \
            scripts/prepare_reference.sh)
        REF_DEP=(--dependency=afterok:"$REF_JOBID")
        echo "  reference job: $REF_JOBID (all sample jobs wait for it)"
    else
        bash scripts/prepare_reference.sh
    fi
fi

submitted=0
echo ""
echo "=== Samples in input/ ==="
skip_log="$(mktemp)"
while IFS=$'\t' read -r sample r1 r2; do
    echo "  $sample"
    if [ "$USE_SLURM" = "1" ]; then
        sbatch "${MAIL_ARGS[@]}" "${REF_DEP[@]}" "${SBATCH_PARTITION_ARGS[@]}" \
            --job-name="baf_$sample" \
            scripts/sample_pipeline.sh "$sample" "$r1" "$r2"
    else
        bash scripts/sample_pipeline.sh "$sample" "$r1" "$r2"
    fi
    (( submitted++ )) || true
done < <(python3 scripts/discover_samples.py input 2>"$skip_log")
sed 's/^/  /' "$skip_log" >&2
skipped=$(wc -l < "$skip_log" | tr -d ' ')
rm -f "$skip_log"

echo ""
if [ "$USE_SLURM" = "1" ]; then
    echo "=== $submitted job(s) submitted | $skipped skipped ==="
else
    echo "=== $submitted sample(s) processed | $skipped skipped ==="
fi
