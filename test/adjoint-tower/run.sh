#!/bin/bash
# The adjoint tower runner: the calculator's frontier law, fourth
# client - grouped by ARCHITECTURAL GATE. The nucleus ladder is
# still checked rung by rung (a failed level blocks dependent work,
# the first absent level closes the frontier and everything above
# reports BLOCKED); the gates are only how the ladder is read and
# reviewed. A gate reports PASS when every level it contains does,
# and UNBUILT while its levels are absent.
#
# There is no result marker until Gate C exists: this tower's
# answer is the total sensitivity df/dp, and nothing below the
# statement has computed one.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

"$here/check_imports.sh" || { echo "└── the import gate refused the tower"; exit 1; }

( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )

# level-directory  gate  label
levels=(
  "level-0-carrier            A 0 carrier"
  "level-1-relation           A 1 relation"
  "level-2-relation-algebra   A 2 relation algebra"
  "level-3-graph              A 3 relational graph"
  "level-4-graph-calculus     A 4 graph calculus"
  "level-5-field-calculus     A 5 field calculus"
  "level-6-discretization     A 6 discretization"
  "level-7-minimization       B 7 minimization"
  "level-8-constitution       B 8 constitution"
  "level-9-statement          C 9 statement"
)

declare -A gate_name=( [A]="structure" [B]="solve + constitution" [C]="statement" )
declare -A gate_state=( [A]="" [B]="" [C]="" )

echo "adjoint sensitivity tower"

frontier_open=1
failed=0
current=""

for idx in "${!levels[@]}"; do
    set -- ${levels[$idx]}
    dir="$1"; gate="$2"; num="$3"; shift 3; label="$num $*"
    dots="$(printf '%.*s' $((24 - ${#label})) '........................')"

    if [ "$gate" != "$current" ]; then
        [ -n "$current" ] && echo "│"
        echo "├── Gate $gate · ${gate_name[$gate]}"
        current="$gate"
    fi

    # last rung of its gate closes the branch
    next_gate=""
    [ $((idx + 1)) -lt ${#levels[@]} ] && next_gate=$(set -- ${levels[$((idx + 1))]}; echo "$2")
    if [ "$next_gate" = "$gate" ]; then tee="├"; else tee="└"; fi

    if [ "$failed" -ne 0 ]; then
        echo "│   $tee── $label $dots SKIPPED (a lower rung failed)"
        gate_state[$gate]="SKIPPED"
        continue
    fi
    if [ "$frontier_open" -eq 0 ]; then
        echo "│   $tee── $label $dots BLOCKED (a lower rung is absent)"
        [ -z "${gate_state[$gate]}" ] && gate_state[$gate]="UNBUILT"
        continue
    fi
    if [ ! -d "$here/$dir" ]; then
        echo "│   $tee── $label $dots ABSENT (the frontier closes here)"
        frontier_open=0
        [ -z "${gate_state[$gate]}" ] && gate_state[$gate]="UNBUILT"
        continue
    fi

    make -C "$here/$dir" clean >/dev/null 2>&1 || true
    if ! make -C "$here/$dir" >/dev/null 2>&1; then
        echo "│   $tee── $label $dots FAIL (build)"
        failed=1
        gate_state[$gate]="FAIL"
        continue
    fi

    ok=1
    ( cd "$here/$dir" && ./run >run.out 2>&1 ) || ok=0
    if [ -x "$here/$dir/check_refusals.sh" ]; then
        ( cd "$here/$dir" && ./check_refusals.sh >>run.out 2>&1 ) || ok=0
    fi

    if [ "$ok" -eq 1 ]; then
        echo "│   $tee── $label $dots PASS"
        [ -z "${gate_state[$gate]}" ] && gate_state[$gate]="PASS"
    else
        echo "│   $tee── $label $dots FAIL"
        sed 's/^/│       /' "$here/$dir/run.out"
        failed=1
        gate_state[$gate]="FAIL"
    fi
done

echo "│"
if [ "$failed" -ne 0 ]; then
    echo "└── the ladder stops at the first failure"
    exit 1
fi
echo "├── Gate A · structure ............... ${gate_state[A]}"
echo "├── Gate B · solve + constitution .... ${gate_state[B]}"
echo "└── Gate C · statement ............... ${gate_state[C]}"
