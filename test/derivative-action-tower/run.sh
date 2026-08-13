#!/bin/bash
# The derivative action tower runner. The calculator's frontier
# law, third client: the FIRST absent rung reports ABSENT and
# closes the frontier; everything above reports BLOCKED,
# unexecuted. A genuine failure reports FAIL, SKIPPED above,
# nonzero exit. Level 7 is a deliberate N/A: this orbit, as
# currently constituted, does not inhabit the minimization radial
# contract - the nucleus levels are available contracts, not a
# compulsory pipeline.
#
# After a full ladder the derivative is read - fail closed - from
# the ninth rung's own output. This tower's answer is a VECTOR, not
# a scalar: exactly one marker, carrying one value per member of
# the statement's independent domain. The runner checks the shape
# and that every token is a number; it never learns what the
# numbers should be.
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
  "N/A:                         level 7  minimization"
  "level-8-derivative-constitution level 8  derivative action"
  "level-9-statement            level 9  statement"
)

echo "derivative action tower"
frontier_open=1
failed=0

for entry in "${levels[@]}"; do
    dir="${entry%% *}"
    label="$(echo "${entry#* }" | sed 's/^ *//')"
    dots="$(printf '%.*s' $((33 - ${#label})) '.................................')"

    if [ "$dir" = "N/A:" ]; then
        echo "├── $label $dots N/A - not inhabited"
        continue
    fi
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
echo "├── Gate A  structure ................ PASS"
if [ "$frontier_open" -eq 0 ]; then
    echo "├── Gate B  numerical duality ........ ABSENT"
    echo "└── Gate C  statement ................ ABSENT"
    exit 0
fi
echo "├── Gate B  numerical duality ........ PASS"
echo "├── Gate C  statement ................ PASS"

# Fail closed: exactly one marker, one number per independent-domain
# member (the statement declares X = {y,x}: two), every token a
# number. The runner knows the SHAPE of the answer, never its value.
out="$here/level-9-statement/run.out"
marks=$(grep -c 'DERIVATIVE_RESULT =' "$out")
values=$(grep -o 'DERIVATIVE_RESULT =.*' "$out" | sed 's/.*DERIVATIVE_RESULT =//')
count=$(echo $values | wc -w)
numeric=1
for tok in $values; do
    case "$tok" in
        ''|*[!0-9-]*) numeric=0 ;;
    esac
done
if [ "$marks" -ne 1 ] || [ "$count" -ne 2 ] || [ "$numeric" -ne 1 ]; then
    echo "└── RUNNER FAILURE: the statement did not report one derivative on X"
    exit 1
fi
echo "└── derivative on X={y,x} ........... [$(echo $values | sed 's/ /, /g')]"
