#!/bin/bash
# The visualization tower's runner.
#
# LEVELS are the implementation architecture; GATES are only review
# checkpoints, and appear here as horizontal separators - never as
# directories, and never in place of a level's own result. Every
# level reports its own status, in order, and a gate line may only
# follow the levels it reviews.
#
# Eight levels are built, and the eighth is RED.
#
# Level 7 asserts that the diagonal of an unchanged matrix does not
# depend on an irrelevant context graph. It does. The ladder therefore
# stops there, and REVIEW GATE B IS NOT PRINTED - a gate reviews the
# levels beneath it, and one of those does not hold.
#
# The failure is the finding, not an accident. Nothing here weakens
# the oracle to make the line green.
#
# Nothing here prints TOWER SEALED, because nothing here has sealed
# anything - the frontier is where this tower's evidence currently
# ends, and the runner's last word is where to look next.
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
  "level-4-graph-calculus     4 structural interpretation"
  "GATE                       A"
  "level-5-field-calculus     5 field calculus"
  "level-6-discretization     6 discretization"
  "level-7-minimization       7 minimization"
)

echo "VISUALIZATION TOWER"
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
    dots="$(printf '%.*s' $((34 - ${#label})) '..................................')"

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
    echo
    echo "    the ladder stops at the first failure"
    echo "    REVIEW GATE B is NOT reached"
    exit 1
fi

echo
echo "    L8 constitution ................... UNBUILT"
echo "    frontier stops here"
echo

# What the tower actually produced, replayed from the levels' own
# output. Every line of it was generated a moment ago from twelve
# occurrences and twelve coefficients; none of it is stored anywhere.
echo "    GENERATED STRUCTURE"
sed -n '/^ ----/,/^ ----/p' "$here/level-4-graph-calculus/run.out" \
    | sed '1d;$d' | sed 's/^/    /'
echo
echo "    GENERATED VALUES"
sed -n '/^ ----/,/^ ----/p' "$here/level-5-field-calculus/run.out" \
    | sed '1d;$d' | sed 's/^/    /'
echo
echo "    PRODUCTION MEASUREMENT"
sed -n '/^ ----/,/^ ----/p' "$here/level-6-discretization/run.out" \
    | sed '1d;$d' | sed 's/^/    /'
