#!/bin/bash
# The partitioned implicit PDE tower runner.
#
# LEVELS are the implementation architecture; GATES are only review
# checkpoints, and appear here as horizontal separators - never as
# directories, and never in place of a level's own result. Every
# level reports its own status, in order, and a gate line may only
# follow the levels it reviews.
#
# The frontier law still holds level by level: a failed level stops
# the ladder, and the first absent level closes the frontier.
#
# After a full ladder the solution field is read - fail closed -
# from Level 9's own output.
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
  "level-3-graph              3 graph"
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

echo "PARTITIONED IMPLICIT PDE TOWER"
echo

frontier_open=1
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
    dots="$(printf '%.*s' $((34 - ${#label})) '..................................')"

    if [ "$failed" -ne 0 ]; then
        echo "    $label $dots SKIPPED"
        continue
    fi
    if [ "$frontier_open" -eq 0 ] || [ ! -d "$here/$dir" ]; then
        echo "    $label $dots ABSENT"
        frontier_open=0
        continue
    fi

    make -C "$here/$dir" clean >/dev/null 2>&1 || true
    if ! make -C "$here/$dir" >/dev/null 2>&1; then
        echo "    $label $dots FAIL (build)"
        failed=1
        continue
    fi

    ok=1
    ( cd "$here/$dir" && ./run >run.out 2>&1 ) || ok=0
    if [ -x "$here/$dir/check_refusals.sh" ]; then
        ( cd "$here/$dir" && ./check_refusals.sh >>run.out 2>&1 ) || ok=0
    fi

    if [ "$ok" -eq 1 ]; then
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
if [ "$frontier_open" -eq 0 ]; then
    echo "    (the tower is incomplete)"
    exit 0
fi

out="$here/level-9-statement/run.out"
marks=$(grep -c 'PARTITIONED_PDE_RESULT =' "$out")
values=$(grep -o 'PARTITIONED_PDE_RESULT =.*' "$out" | sed 's/.*PARTITIONED_PDE_RESULT =//')
if ! marker_ok "$marks" 6 "$values"; then
    echo "    RUNNER FAILURE: the statement did not report one solution field"
    exit 1
fi
echo "    PARTITIONED_PDE_RESULT = $(echo $values)"
echo
echo "    TOWER SEALED."
