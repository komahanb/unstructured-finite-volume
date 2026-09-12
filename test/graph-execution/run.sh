#!/bin/bash
# Verify independent incremental executions over one dependency graph.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    (cd "$here/../.." && ./build.sh >/dev/null)
fi
make -C "$here" clean >/dev/null
make -C "$here" >/dev/null
"$here/run"
declare -A reason=(
    [unpaired_advance]='a driver is paired with its data before it evaluates'
    [unpaired_advance_with]='a driver is paired with its data before it evaluates'
    [unpaired_set]='this driver is not paired with its data'
    [unpaired_clear]='this driver is not paired with its data'
    [outside_set]='a vertex is one the pairing stores'
    [outside_clear]='a vertex is one the pairing stores'
)
refusal_output="$(mktemp)"
trap 'rm -f "$refusal_output"' EXIT
for case_name in "${!reason[@]}"; do
    if "$here/refusal" "$case_name" >"$refusal_output" 2>&1; then
        echo "FAIL : invalid execution admitted: $case_name"
        exit 1
    fi
    if ! grep -Fq "${reason[$case_name]}" "$refusal_output"; then
        cat "$refusal_output"
        echo "FAIL : unrelated refusal: $case_name"
        exit 1
    fi
    echo "PASS : invalid execution refused: $case_name"
done
