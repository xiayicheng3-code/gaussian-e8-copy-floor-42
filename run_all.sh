#!/usr/bin/env bash
# Runs the full Gaussian E8 verification from scratch. Exits nonzero if any check fails.
# Usage: ./run_all.sh            (Gaussian, the subject of this package)
#        ./run_all.sh eis        (optional: the superseded Eisenstein bound, see README)
set -euo pipefail
cd "$(dirname "$0")/verify"
kind="${1:-gauss}"
rm -rf ../cache
start=$(date +%s)
python3 census.py "$kind"
python3 terms.py "$kind"
python3 extras.py "$kind"
python3 thm5.py
echo "ALL CHECKS PASSED ($kind) in $(( $(date +%s) - start ))s"
