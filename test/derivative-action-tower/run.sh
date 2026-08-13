#!/bin/bash
# The derivative action tower runner - Gate A: the structural
# rungs, and only those. The calculator's frontier law, third
# client: the FIRST absent rung reports ABSENT and closes the
# frontier; everything above reports BLOCKED, unexecuted. A genuine
# failure reports FAIL, SKIPPED above, nonzero exit. There is
# deliberately NO numerical result marker: Gate A certifies
# structure, and fabricating a number here would claim a derivative
# nothing has computed. Gates B and C are UNBUILT.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

"$here/check_imports.sh" || { echo "└── the import gate refused the tower"; exit 1; }

( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )

levels=(
  "level-0-carrier              level 0  carriers"
  "level-1-relation             level 1  relation"
  "level-2-relation-algebra     level 2  relation algebra"
  "level-3-graph                level 3  relational graph"
  "level-4-graph-calculus       level 4  graph calculus"
  "level-5-field-calculus       level 5  field calculus"
  "level-6-derivative-structure level 6  derivative structure"
)

echo "derivative action tower - Gate A"
frontier_open=1
failed=0

for entry in "${levels[@]}"; do
    dir="${entry%% *}"
    label="$(echo "${entry#* }" | sed 's/^ *//')"
    dots="$(printf '%.*s' $((33 - ${#label})) '.................................')"

    if [ "$failed" -ne 0 ]; then
        echo "├── $label $dots SKIPPED (a lower rung failed)"
        continue
    fi
    if [ "$frontier_open" -eq 0 ]; then
        echo "├── $label $dots BLOCKED (a lower rung is absent)"
        continue
    fi
    if [ ! -d "$here/$dir" ]; then
        echo "├── $label $dots ABSENT (the frontier closes here)"
        frontier_open=0
        continue
    fi

    make -C "$here/$dir" clean >/dev/null 2>&1 || true
    if ! make -C "$here/$dir" >/dev/null 2>&1; then
        echo "├── $label $dots FAIL (build)"
        failed=1
        continue
    fi

    ok=1
    ( cd "$here/$dir" && ./run >run.out 2>&1 ) || ok=0
    if [ -x "$here/$dir/check_refusals.sh" ]; then
        ( cd "$here/$dir" && ./check_refusals.sh >>run.out 2>&1 ) || ok=0
    elif [ -x "$here/$dir/refusal" ]; then
        if ( cd "$here/$dir" && ./refusal >>run.out 2>&1 ); then
            ok=0
        fi
    fi

    if [ "$ok" -eq 1 ]; then
        echo "├── $label $dots PASS"
    else
        echo "├── $label $dots FAIL"
        sed 's/^/│       /' "$here/$dir/run.out"
        failed=1
    fi
done

if [ "$failed" -ne 0 ]; then
    echo "└── the ladder stops at the first failure"
    exit 1
fi
if [ "$frontier_open" -eq 1 ]; then
    echo "└── Gate A  structural derivative ..... PASS"
else
    echo "└── every certified rung holds its truths"
fi
