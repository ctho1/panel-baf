#!/bin/bash
# Loads the software for a job: PALMA modules if a `module` command exists,
# otherwise (local testing) the conda env. Source after scripts/config.sh.

aligner_index_exts() {   # index file suffixes of the configured aligner
    case "$ALIGNER" in
        bwa-mem2) echo "0123 amb ann bwt.2bit.64 pac" ;;
        bwa)      echo "amb ann bwt pac sa" ;;
        *) echo "Unknown ALIGNER: $ALIGNER (expected: bwa-mem2|bwa)" >&2; return 1 ;;
    esac
}

load_software() {
    if [ "$USE_CONDA" != "true" ] && command -v module >/dev/null 2>&1; then
        local mods
        # de-duplicated, order-preserving list (palma/2024a, GCC/... appear in several lists)
        mods=$(printf '%s\n' $BWA_MODULES $SAMTOOLS_MODULES $BCFTOOLS_MODULES $PYTHON_MODULES $FASTP_MODULES | awk 'NF && !seen[$0]++')
        module purge
        # shellcheck disable=SC2086
        ml $mods
        echo "Software: PALMA modules"
    else
        if [ ! -x "$ENV_BAF/bin/python" ]; then
            echo "No module system and no conda env ($ENV_BAF) -- run 'bash scripts/setup_envs.sh'." >&2
            return 1
        fi
        export PATH="$ENV_BAF/bin:$PATH"
        echo "Software: conda env $ENV_BAF"
    fi
    PY="$(command -v python3)"
    export PY
    local missing=()
    for tool in "$ALIGNER" samtools bcftools bgzip tabix python3; do
        command -v "$tool" >/dev/null 2>&1 || missing+=("$tool")
    done
    "$PY" -c "import numpy, scipy, pandas, matplotlib" 2>/dev/null || missing+=("python: numpy/scipy/pandas/matplotlib")
    if [ ${#missing[@]} -gt 0 ]; then
        echo "Missing software: ${missing[*]} -- check the module lists in scripts/config.sh" >&2
        echo "(bash scripts/check_palma_modules.sh lists the available modules)." >&2
        return 1
    fi
}
