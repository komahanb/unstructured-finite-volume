#!/bin/bash
# The time integration tower's runner.
#
# LEVELS are the implementation architecture; GATES are only review
# checkpoints, and appear here as horizontal separators - never as
# directories, and never in place of a level's own result. Every
# level reports its own status, in order, and a gate line may only
# follow the levels it reviews.
#
# All ten levels are built, and after a full ladder the answer is
# read - fail closed - from Level 9's own output. The tower says
# TOWER SEALED once, here, and only after every level has reported
# for itself and the result contract has accepted the marker.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

"$here/check_imports.sh" --selftest || { echo "└── the import gate refused itself"; exit 1; }
"$here/check_imports.sh" || { echo "└── the import gate refused the tower"; exit 1; }

. "$here/check_marker.sh"
"$here/check_marker.sh" --selftest || { echo "└── the result contract refused itself"; exit 1; }

( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )

levels=(
  "level-0-carrier            0 carrier"
  "level-1-relation           1 relation"
  "level-2-relation-algebra   2 relation algebra"
  "level-3-graph              3 relational graph"
  "level-4-graph-calculus     4 graph calculus"
  "GATE                       A"
  "level-5-field-calculus     5 field calculus"
  "level-6-discretization     6 discretization"
  "level-7-minimization       7 minimization"
  "GATE                       B"
  "level-8-constitution       8 constitution"
  "level-9-statement          9 statement"
  "GATE                       C"
)

echo "TIME INTEGRATION TOWER"
echo

failed=0

for entry in "${levels[@]}"; do
    set -- $entry
    dir="$1"; shift

    if [ "$dir" = "GATE" ]; then
        echo
        echo "    ===== REVIEW GATE $1 ====="
        echo
        continue
    fi

    num="$1"; shift; label="L$num $*"
    dots="$(printf '%.*s' $((35 - ${#label})) '...................................')"

    if [ "$failed" -ne 0 ]; then
        echo "    $label $dots SKIPPED"
        continue
    fi
    if [ ! -d "$here/$dir" ]; then
        echo "    $label $dots ABSENT"
        failed=1
        continue
    fi

    make -C "$here/$dir" clean >/dev/null 2>&1 || true
    if ! make -C "$here/$dir" >/dev/null 2>&1; then
        echo "    $label $dots FAIL (build)"
        failed=1
        continue
    fi

    if ( cd "$here/$dir" && ./run >run.out 2>&1 ); then
        echo "    $label $dots PASS"
    else
        echo "    $label $dots FAIL"
        sed 's/^/        /' "$here/$dir/run.out"
        failed=1
    fi
done

if [ "$failed" -ne 0 ]; then
    echo "    the ladder stops at the first failure"
    exit 1
fi

out="$here/level-9-statement/run.out"
marks=$(grep -c 'TIME_INTEGRATION_RESULT =' "$out")
values=$(grep -o 'TIME_INTEGRATION_RESULT =.*' "$out" | sed 's/.*TIME_INTEGRATION_RESULT =//')
if ! marker_ok "$marks" 2 "$values"; then
    echo "    RUNNER FAILURE: the statement did not report one state on Q"
    exit 1
fi
echo "    TIME_INTEGRATION_RESULT = $(echo $values)"
echo
echo "    TOWER SEALED."
