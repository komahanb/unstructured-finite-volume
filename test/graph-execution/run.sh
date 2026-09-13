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
    [unpaired_advance]='advance was called before this driver was paired with its data'
    [unpaired_advance_with]='advance_with was called before this driver was paired with its data'
    [unpaired_set]='driver_set_rule was called before this driver was paired with its data'
    [unpaired_clear]='driver_clear_rule was called before this driver was paired with its data'
    [outside_set]='vertex must be one the pairing stores'
    [outside_clear]='vertex must be one the pairing stores'
    [outside_position]='must be a completed count of the stored order'
    [outside_live]='must be a completed count of the stored order'
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
