#!/bin/bash
# The calculator tower: build the library once, then climb the levels
# in order, stopping the ladder at the first failure. Absent levels
# are reported absent - a missing level stays visibly unfinished
# (CALCULATOR.md 21); nothing above a failure is certified.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )

levels=(
  "level-0-carrier          level 0  carriers"
  "level-1-relation         level 1  relation"
  "level-2-relation-algebra level 2  relation algebra"
  "level-3-graph            level 3  relational graph"
  "level-4-graph-calculus   level 4  graph calculus"
  "level-5-field-calculus   level 5  field calculus"
  "level-6-discretization   level 6  discretization"
  "level-7-minimization     level 7  minimization"
  "level-8-constitution     level 8  constitution"
  "level-9-statement        level 9  statement"
)

echo "calculator tower"
broken=0
for entry in "${levels[@]}"; do
    dir="${entry%% *}"
    label="$(echo "${entry#* }" | sed 's/^ *//')"
    dots="$(printf '%.*s' $((28 - ${#label})) '............................')"
    if [ ! -d "$here/$dir" ]; then
        echo "├── $label $dots ABSENT"
        continue
    fi
    if [ "$broken" -ne 0 ]; then
        echo "├── $label $dots SKIPPED (lower level failed)"
        continue
    fi
    make -C "$here/$dir" clean >/dev/null 2>&1 || true
    if ! make -C "$here/$dir" >/dev/null 2>&1; then
        echo "├── $label $dots FAIL (build)"
        broken=1
        continue
    fi
    ok=1
    ( cd "$here/$dir" && ./run >run.out 2>&1 ) || ok=0
    if [ -x "$here/$dir/refusal" ]; then
        if ( cd "$here/$dir" && ./refusal >>run.out 2>&1 ); then
            ok=0   # a refusal that survives is a failure
        fi
    fi
    if [ "$ok" -eq 1 ]; then
        echo "├── $label $dots PASS"
    else
        echo "├── $label $dots FAIL"
        sed 's/^/│       /' "$here/$dir/run.out"
        broken=1
    fi
done
if [ "$broken" -ne 0 ]; then
    echo "└── the ladder stops at the first failure"
    exit 1
fi
echo "└── every present level holds its truths"
