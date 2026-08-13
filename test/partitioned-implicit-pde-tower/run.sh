#!/bin/bash
# The partitioned implicit PDE tower runner. Development is grouped
# by GATE, not by nucleus rung: this tower consumes most levels
# rather than reconstructing them, and the README's Rosetta table
# maps each gate's truths back onto the levels it uses.
#
# The frontier law still holds gate by gate: a failed gate stops the
# ladder, and the first absent gate closes the frontier so nothing
# above it is claimed. Gates B and C report UNBUILT until they exist.
#
# After a full ladder the solution field is read - fail closed - from
# Gate C's own output. There is no marker before then: this tower's
# answer is a six-member field, and nothing below the statement has
# computed one.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

"$here/check_imports.sh" || { echo "└── the import gate refused the tower"; exit 1; }

if [ -f "$here/check_marker.sh" ]; then
    . "$here/check_marker.sh"
    "$here/check_marker.sh" --selftest || { echo "└── the result contract refused itself"; exit 1; }
fi

( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )

gates=(
  "gate-a-partition   A partition / ownership / transport"
  "gate-b-operator    B topology-consuming action"
  "gate-c-statement   C partitioned implicit statement"
)

echo "partitioned implicit pde tower"

frontier_open=1
failed=0

for entry in "${gates[@]}"; do
    set -- $entry
    dir="$1"; letter="$2"; shift 2; label="Gate $letter · $*"
    dots="$(printf '%.*s' $((46 - ${#label})) '..............................................')"

    if [ "$failed" -ne 0 ]; then
        echo "├── $label $dots SKIPPED"
        continue
    fi
    if [ "$frontier_open" -eq 0 ] || [ ! -d "$here/$dir" ]; then
        echo "├── $label $dots UNBUILT"
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
if [ "$frontier_open" -eq 0 ]; then
    echo "└── solution field ............................. (unbuilt)"
    exit 0
fi

out="$here/gate-c-statement/run.out"
marks=$(grep -c 'PARTITIONED_PDE_RESULT =' "$out")
values=$(grep -o 'PARTITIONED_PDE_RESULT =.*' "$out" | sed 's/.*PARTITIONED_PDE_RESULT =//')
if ! marker_ok "$marks" 6 "$values"; then
    echo "└── RUNNER FAILURE: the statement did not report one solution field"
    exit 1
fi
echo "└── solution field on V(G) ..................... [$(echo $values | sed 's/  */, /g')]"
