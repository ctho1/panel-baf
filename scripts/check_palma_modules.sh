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

echo "### 3. Test of the module lists in scripts/config.sh"
echo "BWA_MODULES=$BWA_MODULES"
echo "SAMTOOLS_MODULES=$SAMTOOLS_MODULES"
echo "BCFTOOLS_MODULES=$BCFTOOLS_MODULES"
echo "PYTHON_MODULES=$PYTHON_MODULES"
echo "FASTP_MODULES=$FASTP_MODULES"
if load_software; then
    echo "OK: all required programs found"
else
    echo "ERROR while loading (see above)"
fi
module list 2>&1
for tool in bwa-mem2 bwa samtools bcftools bgzip tabix fastp python3; do
    printf '%-9s %s\n' "$tool" "$(command -v "$tool" || echo '-')"
done
samtools --version 2>/dev/null | head -1
bcftools --version 2>/dev/null | head -1
python3 - <<'EOF' 2>&1
import sys
print("python", sys.version.split()[0])
for m in ("numpy", "scipy", "pandas", "matplotlib"):
    try:
        mod = __import__(m)
        print(f"{m:10s} {mod.__version__}")
    except Exception as exc:
        print(f"{m:10s} MISSING ({exc})")
EOF
echo
echo "### 4. Reference (shared with mNGS-pipeline)"
ls -la "$EXISTING_HG19_FASTA"* 2>&1 | head -12
