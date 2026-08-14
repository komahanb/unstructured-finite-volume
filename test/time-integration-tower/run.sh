#!/bin/bash
# The time integration tower's runner.
#
# LEVELS are the implementation architecture; GATES are only review
# checkpoints, and appear here as horizontal separators - never as
# directories, and never in place of a level's own result. Every
# level reports its own status, in order, and a gate line may only
# follow the levels it reviews.
#
# THIS TOWER IS NOT SEALED, and this script must never say it is.
# Levels 0-4 are built; Levels 5-9 are not. The runner therefore
# ends by naming the frontier out loud rather than by reading a
# result marker - because Gate A produces no numerical result, and
# a result contract for a tower with no result would be a claim it
# has not earned.
#
# When Level 5 lands, whoever builds it edits the frontier lines
# below. That is honest: the runner declares where the frontier is,
# so moving the frontier is an edit somebody makes on purpose.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

"$here/check_imports.sh" --selftest || { echo "└── the import gate refused itself"; exit 1; }
"$here/check_imports.sh" || { echo "└── the import gate refused the tower"; exit 1; }

( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )

levels=(
  "level-0-carrier            0 carrier"
  "level-1-relation           1 relation"
  "level-2-relation-algebra   2 relation algebra"
  "level-3-graph              3 relational graph"
  "level-4-graph-calculus     4 graph calculus"
)

echo "TIME INTEGRATION TOWER"
echo

failed=0

for entry in "${levels[@]}"; do
    set -- $entry
    dir="$1"; shift
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

echo
echo "    ===== REVIEW GATE A ====="
echo
echo "    L5 field calculus .................. UNBUILT"
echo "    frontier stops here"
