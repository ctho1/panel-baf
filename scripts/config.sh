#!/bin/bash
# Central config. Paths are relative to the repo root (portable). Override
# any value via env var, e.g.:
#   PANEL_BED=/path/to/panel.bed bash run_pipeline.sh

: "${MAIL_USER:=}"

# SLURM partition selection (same logic as mNGS-pipeline): normal + requeue,
# ranked by idle CPUs; requeue jobs resume at the last finished step.
: "${SLURM_AUTO_PARTITION:=true}"
: "${SLURM_CPU_PARTITIONS:=normal,requeue}"
: "${SLURM_FALLBACK_PARTITION:=normal}"

# Panel target BED (hg19, UCSC "chr" names). Not part of the repository --
# copy the panel BED to this path or override PANEL_BED.
: "${PANEL_BED:=$REPO_ROOT/references/panel_targets.bed}"

# Aligner + samtools from PALMA modules (conda's samtools hits a libdeflate
# conflict there); loaded via `module purge && ml $MINIMAP2_MODULES $SAMTOOLS_MODULES`.
# Locally (no `module` command) the conda env provides both.
: "${MINIMAP2_MODULES:=palma/2024a GCCcore/13.3.0 minimap2/2.29}"
: "${SAMTOOLS_MODULES:=palma/2024a GCC/13.3.0 SAMtools/1.21}"

# Conda env (see scripts/setup_envs.sh): fastp, bcftools, mosdepth, python stack
: "${CONDA_ENV_DIR:=$REPO_ROOT/conda_envs}"
ENV_BAF="$CONDA_ENV_DIR/baf"

# hg19: shared fasta on PALMA (symlinked if present), else UCSC download.
# Only primary chromosomes (chr1-22, X, Y, M) are used -- alt/hap contigs
# would push reads in MHC & co. to MAPQ 0.
: "${EXISTING_HG19_FASTA:=/cloud/wwu1/e_np_ngs/references/genomes/hg19.fa}"
: "${HG19_FASTA_URL:=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz}"
REF_DIR="$REPO_ROOT/references"
HG19_PRIMARY="$REF_DIR/hg19_primary.fa"
HG19_MMI_SR="$REF_DIR/hg19_primary.sr.mmi"

# Annotation / population resources (downloaded once by prepare_reference.sh)
: "${KG_SITES_URL:=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz}"
: "${UCSC_HG19_DB_URL:=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database}"
: "${ICHORCNA_EXTDATA_URL:=https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata}"
SNP_PANEL="$REF_DIR/snps_1kg_maf05.tsv.gz"   # chrom pos ref alt AF EAS AMR AFR EUR SAS
: "${SNP_MIN_MAF:=0.05}"

# Threads (local override; on PALMA these come from SLURM_CPUS_PER_TASK)
: "${THREADS_DEFAULT:=8}"

# Genes for the focal zoom figure: groups separated by ";", genes within a group by ","
: "${FOCAL_GENES:=EGFR;CDKN2A,CDKN2B;PTEN;PDGFRA,KIT}"

export MAIL_USER PANEL_BED
export SLURM_AUTO_PARTITION SLURM_CPU_PARTITIONS SLURM_FALLBACK_PARTITION
export MINIMAP2_MODULES SAMTOOLS_MODULES ENV_BAF
export EXISTING_HG19_FASTA HG19_FASTA_URL REF_DIR HG19_PRIMARY HG19_MMI_SR
export KG_SITES_URL UCSC_HG19_DB_URL ICHORCNA_EXTDATA_URL SNP_PANEL SNP_MIN_MAF
export THREADS_DEFAULT FOCAL_GENES
