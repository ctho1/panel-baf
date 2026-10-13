#!/bin/bash
# Legt die Conda-Env an (einmalig). Linux/PALMA: minimap2 + samtools kommen
# aus dem Modulsystem; macOS (lokaler Test): zusätzlich aus Conda.
#
# Wie in mNGS-pipeline: `mamba create` mit --override-channels
# --no-channel-priority (umgeht den libmamba-Solver-Bug bei
# channel_priority=strict auf PALMA).
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
    echo "Conda-Installation nicht gefunden unter: $CONDA_BASE" >&2
    echo "Setze CONDA_BASE explizit, z.B.: CONDA_BASE=\$(conda info --base) bash scripts/setup_envs.sh" >&2
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

PKGS=("python=3.12" "numpy" "pandas" "scipy" "statsmodels" "matplotlib" "pysam"
      "bcftools>=1.20" "htslib>=1.20" "fastp" "mosdepth")
if [ "$(uname -s)" != "Linux" ]; then
    PKGS+=("minimap2>=2.28" "samtools>=1.20")
fi

mkdir -p "$ENV_DIR"
echo "=== baf (${PKGS[*]}) ==="
if [ ! -d "$ENV_DIR/baf" ]; then
    "$MAMBA" create -y --override-channels --no-channel-priority -c conda-forge -c bioconda \
        -p "$ENV_DIR/baf" "${PKGS[@]}"
else
    echo "  bereits vorhanden, überspringe"
fi
echo "=== Fertig. Env liegt unter $ENV_DIR/baf ==="
