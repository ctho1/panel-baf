#!/bin/bash
# Checks on PALMA which modules are available for the pipeline and tests the module lists
# configured in scripts/config.sh. Run on a login node (~1 min) and report the output:
#
#   bash scripts/check_palma_modules.sh 2>&1 | tee palma_modules.txt
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
source scripts/config.sh
source scripts/software.sh

if ! command -v module >/dev/null 2>&1; then
    echo "No 'module' command found -- please run this on PALMA." >&2
    exit 1
fi

echo "### Host: $(hostname)  date: $(date -Iseconds)"
echo "### Lmod: $(module --version 2>&1 | head -1)"
echo

echo "### 1. Available versions (module -t spider)"
for name in bwa-mem2 BWA SAMtools BCFtools HTSlib fastp Python SciPy-bundle matplotlib; do
    echo "--- $name"
    module -t spider "$name" 2>&1 | grep -i "^$name/" | sort -V | tail -8
done
echo

echo "### 2. Prerequisites of the newest versions (module spider <name>/<version>)"
for name in bwa-mem2 SAMtools BCFtools HTSlib fastp SciPy-bundle matplotlib; do
    latest=$(module -t spider "$name" 2>&1 | grep -i "^$name/" | sort -V | tail -1)
    [ -n "$latest" ] || { echo "--- $name: not found"; continue; }
    echo "--- $latest"
    module spider "$latest" 2>&1 | grep -A6 "You will need to load" | grep -v "^ *$" | head -6
done
echo

echo "### 3. Test of the module lists per step (loaded exactly as in sample_pipeline.sh)"
show() { for tool in "$@"; do printf '  %-9s %s\n' "$tool" "$(command -v "$tool" || echo MISSING)"; done; }

echo "--- alignment: $BWA_MODULES $SAMTOOLS_MODULES"
use_modules $BWA_MODULES $SAMTOOLS_MODULES
echo "  PATH length: ${#PATH} characters"
show "$ALIGNER" samtools
[ "$ALIGNER" != "bwa-mem2" ] || echo "  bwa-mem2 version: $(bwa-mem2 version 2>&1 | tail -1)"
samtools --version 2>/dev/null | head -1 | sed 's/^/  /'

echo "--- allele counts: $BCFTOOLS_MODULES"
use_modules $BCFTOOLS_MODULES
show bcftools bgzip tabix
bcftools --version 2>/dev/null | head -1 | sed 's/^/  /'

echo "--- Python steps: $PYTHON_MODULES $SAMTOOLS_MODULES"
use_modules $PYTHON_MODULES $SAMTOOLS_MODULES
show python3 samtools
python3 - <<'EOF' 2>&1
import sys
print("  python", sys.version.split()[0])
for m in ("numpy", "scipy", "pandas", "matplotlib"):
    try:
        mod = __import__(m)
        print(f"  {m:10s} {mod.__version__}")
    except Exception as exc:
        print(f"  {m:10s} MISSING ({exc})")
EOF

echo "--- fastp (optional): ${FASTP_MODULES:-not configured}"
if [ -n "$FASTP_MODULES" ]; then use_modules $FASTP_MODULES; show fastp; fi
echo
echo "### 4. Shared hg19 reference (EXISTING_HG19_FASTA)"
ls -la "$EXISTING_HG19_FASTA"* 2>&1 | head -12
