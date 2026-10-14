#!/bin/bash
# Software per step: right before each step `module purge && ml ...`
# loads only the modules that step needs. (Loading all lists at once makes PATH so long that
# the bwa-mem2 2.2.1 launcher aborts with "ERROR: prefix is too long!".)
# Without a `module` command (local testing) the conda env is used. Source after scripts/config.sh.

aligner_index_exts() {   # index file suffixes of the configured aligner
    case "$ALIGNER" in
        bwa-mem2) echo "0123 amb ann bwt.2bit.64 pac" ;;
        bwa)      echo "amb ann bwt pac sa" ;;
        *) echo "Unknown ALIGNER: $ALIGNER (expected: bwa-mem2|bwa)" >&2; return 1 ;;
    esac
}

use_modules() {   # use_modules $X_MODULES $Y_MODULES ...  (empty: purge only)
    if [ "$USE_CONDA" != "true" ] && command -v module >/dev/null 2>&1; then
        module purge
        [ $# -eq 0 ] || ml "$@"
        echo "Modules: ${*:-none}"
    else
        if [ ! -x "$ENV_BAF/bin/python" ]; then
            echo "No module system and no conda env ($ENV_BAF) -- run 'bash scripts/setup_envs.sh'." >&2
            return 1
        fi
        case ":$PATH:" in *":$ENV_BAF/bin:"*) ;; *) export PATH="$ENV_BAF/bin:$PATH" ;; esac
    fi
}
