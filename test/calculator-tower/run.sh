#!/bin/bash
# The calculator tower runner. The law it enforces:
#
#     L_n ABSENT or FAIL  =>  L_{m>n} cannot be certified
#
# The FIRST missing rung reports ABSENT and closes the frontier;
# everything above it reports BLOCKED and never executes, even if
# its directory physically exists. A genuine failure reports FAIL
# and everything above it reports SKIPPED. An incomplete tower is
# not an incorrect tower: only FAIL makes this runner fail.
#
# Before anything runs, the import gate audits every calculator
# source against its level's allowlist (check_imports.sh): a
# forbidden `use` is a certification failure, mechanically.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

"$here/check_imports.sh" || { echo "└── the import gate refused the tower"; exit 1; }

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
frontier_open=1   # no rung missing yet
failed=0          # a rung genuinely failed

for entry in "${levels[@]}"; do
    dir="${entry%% *}"
    label="$(echo "${entry#* }" | sed 's/^ *//')"
    dots="$(printf '%.*s' $((28 - ${#label})) '............................')"

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

    # Refusals: a level may carry its own message-checking script;
    # otherwise a bare refusal binary must simply die.
    if [ -x "$here/$dir/check_refusals.sh" ]; then
        ( cd "$here/$dir" && ./check_refusals.sh >>run.out 2>&1 ) || ok=0
    elif [ -x "$here/$dir/refusal" ]; then
        if ( cd "$here/$dir" && ./refusal >>run.out 2>&1 ); then
            ok=0   # a refusal that survives is a failure
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
if [ "$frontier_open" -eq 1 ] && [ -f "$here/level-9-statement/run.out" ]; then
    result=$(grep -o 'CALCULATOR_RESULT = .*' "$here/level-9-statement/run.out" | awk '{print $3}')
    echo "└── mathematical result ......... ${result}"
else
    echo "└── every certified rung holds its truths"
fi
