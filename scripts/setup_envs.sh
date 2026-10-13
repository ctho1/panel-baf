#!/bin/bash
# Only for systems WITHOUT a module system (e.g. local testing on a Mac). On PALMA all
# software comes from the modules in scripts/config.sh -- no conda env is needed there.
#
# As in mNGS-pipeline: `mamba create` with --override-channels --no-channel-priority
# (works around a libmamba solver bug with channel_priority=strict).
set -euo pipefail

if [ -z "${CONDA_BASE:-}" ]; then
    if command -v conda >/dev/null 2>&1; then
        CONDA_BASE="$(conda info --base)"
    else
        CONDA_BASE="$HOME/miniforge3"
    fi
fi
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="${CONDA_ENV_DIR:-$REPO_ROOT/conda_envs}"

if [ ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    echo "Conda installation not found at: $CONDA_BASE" >&2
    exit 1
fi
source "$CONDA_BASE/etc/profile.d/conda.sh"
if command -v mamba >/dev/null 2>&1; then
    MAMBA="$(command -v mamba)"
elif [ -x "$CONDA_BASE/bin/mamba" ]; then
    MAMBA="$CONDA_BASE/bin/mamba"
else
    MAMBA="$CONDA_BASE/bin/conda"
fi

# Classic bwa instead of bwa-mem2: building a bwa-mem2 index for hg19 needs ~90 GB RAM,
# so local runs use ALIGNER=bwa.
PKGS=("python=3.12" "numpy" "pandas" "scipy" "matplotlib" "samtools>=1.20" "bcftools>=1.20" "htslib>=1.20"
      "bwa" "fastp")

mkdir -p "$ENV_DIR"
echo "=== baf (${PKGS[*]}) ==="
if [ ! -d "$ENV_DIR/baf" ]; then
    "$MAMBA" create -y --override-channels --no-channel-priority -c conda-forge -c bioconda \
        -p "$ENV_DIR/baf" "${PKGS[@]}"
else
    echo "  already present, skipping"
fi
echo "=== Done. Env at $ENV_DIR/baf -- run locally with ALIGNER=bwa ==="
