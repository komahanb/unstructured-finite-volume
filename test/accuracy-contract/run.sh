#!/bin/bash
# The accuracy contract: analytic and manufactured references, temporal and
# spatial refinement, derivative accuracy, conservation, solver residuals
# and transpose identities, each declared with its order or floor and its
# justification (README.md, cases.py). contract.py drives the application
# executable, reads its records, measures the orders and writes
# results/<set>/summary.json. The required set must pass; the rejection
# set must fail case by case with the status each declares. Pass
# --exploratory to also run the exploratory set, whose statuses are
# reported but never fail this runner.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$here/../.." && ./build.sh >/dev/null 2>&1 && ./application/build.sh >/dev/null 2>&1 )
fi

results="${UFVM_ACCURACY_RESULTS:-$here/results}"
exploratory=0
for argument in "$@"; do
    [ "$argument" = "--exploratory" ] && exploratory=1
done

python3 "$here/contract.py" --set required --results "$results/required"
python3 "$here/contract.py" --set rejection --results "$results/rejection"
if [ "$exploratory" -eq 1 ]; then
    python3 "$here/contract.py" --set exploratory --results "$results/exploratory" || true
fi
echo " accuracy contract: required declarations met and every rejection reported"
