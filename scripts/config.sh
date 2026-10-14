#!/bin/bash
# Central config. Paths are relative to the repo root (portable). Override
# any value via env var, e.g.:
#   PANEL_BED=/path/to/panel.bed bash run_pipeline.sh

# SLURM partition selection: normal + requeue,
# ranked by idle CPUs; requeue jobs resume at the last finished step.
: "${SLURM_AUTO_PARTITION:=true}"
: "${SLURM_CPU_PARTITIONS:=normal,requeue}"
: "${SLURM_FALLBACK_PARTITION:=normal}"

# Panel target BED (hg19, UCSC "chr" names). Default: the NPMS v4 panel shipped with the repository
# (Schnorrenberg et al. 2026, doi:10.1111/nan.70070); set PANEL_BED for other panels.
: "${PANEL_BED:=$REPO_ROOT/panel/NPMS_v4_targets_hg19.bed}"

# ── Software ────────────────────────────────────────────────────────────────
# On PALMA everything comes from the module system (no conda needed); each step
# loads only its own module lists via `module purge && ml ...`
# (alignment: BWA + SAMTOOLS, see scripts/software.sh). Check availability with
# scripts/check_palma_modules.sh.
# Without a `module` command (local tests) the conda env from
# scripts/setup_envs.sh is used instead; USE_CONDA=true forces conda.
: "${ALIGNER:=bwa-mem2}"     # bwa-mem2 (index *.bwt.2bit.64) | bwa (classic index *.bwt/*.sa)
: "${BWA_MODULES:=palma/2024a GCC/13.3.0 bwa-mem2/2.2.1}"
: "${SAMTOOLS_MODULES:=palma/2024a GCC/13.3.0 SAMtools/1.21}"
: "${BCFTOOLS_MODULES:=palma/2024a GCC/13.3.0 BCFtools/1.21 HTSlib/1.21}"
: "${PYTHON_MODULES:=palma/2024a GCC/13.3.0 SciPy-bundle/2024.05 matplotlib/3.9.2}"
: "${FASTP_MODULES:=}"       # optional adapter/polyG trimming; empty -> no trimming (bwa soft-clips adapters)
: "${USE_CONDA:=auto}"
: "${CONDA_ENV_DIR:=$REPO_ROOT/conda_envs}"
ENV_BAF="$CONDA_ENV_DIR/baf"

# ── Reference ───────────────────────────────────────────────────────────────
# Shared hg19 + bwa-mem2 index on PALMA (/cloud, UCSC "chr" names);
# prepare_reference.sh symlinks fasta, .fai and index into references/.
: "${EXISTING_HG19_FASTA:=/cloud/wwu1/e_np_ngs/references/genomes/hg19.fa}"
: "${HG19_FASTA_URL:=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz}"
REF_DIR="$REPO_ROOT/references"
HG19_FASTA="$REF_DIR/hg19.fa"

# Annotation / population resources (downloaded once by prepare_reference.sh)
: "${KG_SITES_URL:=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz}"
: "${UCSC_HG19_DB_URL:=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database}"
: "${ICHORCNA_EXTDATA_URL:=https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata}"
SNP_PANEL="$REF_DIR/snps_1kg_maf05.tsv.gz"   # chrom pos ref alt AF EAS AMR AFR EUR SAS
: "${SNP_MIN_MAF:=0.05}"

# Threads (local override; on PALMA these come from SLURM_CPUS_PER_TASK)
: "${THREADS_DEFAULT:=8}"

# Figure: final print width in mm (174 = Springer double column, 180 = Nature)
: "${FIG_WIDTH_MM:=174}"

export PANEL_BED
export SLURM_AUTO_PARTITION SLURM_CPU_PARTITIONS SLURM_FALLBACK_PARTITION
export ALIGNER BWA_MODULES SAMTOOLS_MODULES BCFTOOLS_MODULES PYTHON_MODULES FASTP_MODULES USE_CONDA ENV_BAF
export EXISTING_HG19_FASTA HG19_FASTA_URL REF_DIR HG19_FASTA
export KG_SITES_URL UCSC_HG19_DB_URL ICHORCNA_EXTDATA_URL SNP_PANEL SNP_MIN_MAF
export THREADS_DEFAULT FIG_WIDTH_MM
