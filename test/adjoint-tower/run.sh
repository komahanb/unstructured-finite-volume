#!/bin/bash
# The adjoint tower runner: the calculator's frontier law, fourth
# client - grouped by ARCHITECTURAL GATE. The nucleus ladder is
# still checked rung by rung (a failed level blocks dependent work,
# the first absent level closes the frontier and everything above
# reports BLOCKED); the gates are only how the ladder is read and
# reviewed. A group reports PASS when every level it contains does,
# and UNBUILT while its levels are absent.
#
# After a full ladder the sensitivity is read - fail closed - from
# the ninth rung's own output: exactly one marker carrying one
# real token. The contract lives in check_marker.sh, which
# self-tests before the ladder runs; the runner validates shape and
# syntax only and never learns what the number should be.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

"$here/check_imports.sh" || { echo "└── the import group refused the tower"; exit 1; }

. "$here/check_marker.sh"
"$here/check_marker.sh" --selftest || { echo "└── the result contract refused itself"; exit 1; }

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
fi

# level-directory  group  label
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

declare -A group_name=( [A]="structure" [B]="solve + constitution" [C]="statement" )
declare -A group_state=( [A]="" [B]="" [C]="" )

echo "adjoint sensitivity tower"

frontier_open=1
failed=0
current=""

for idx in "${!levels[@]}"; do
    set -- ${levels[$idx]}
    dir="$1"; group="$2"; num="$3"; shift 3; label="$num $*"
    dots="$(printf '%.*s' $((24 - ${#label})) '........................')"

    if [ "$group" != "$current" ]; then
        [ -n "$current" ] && echo "│"
        echo "├── Gate $group · ${group_name[$group]}"
        current="$group"
    fi

    # last rung of its group closes the branch
    next_group=""
    [ $((idx + 1)) -lt ${#levels[@]} ] && next_group=$(set -- ${levels[$((idx + 1))]}; echo "$2")
    if [ "$next_group" = "$group" ]; then tee="├"; else tee="└"; fi

    if [ "$failed" -ne 0 ]; then
        echo "│   $tee── $label $dots SKIPPED (a lower rung failed)"
        group_state[$group]="SKIPPED"
        continue
    fi
    if [ "$frontier_open" -eq 0 ]; then
        echo "│   $tee── $label $dots BLOCKED (a lower rung is absent)"
        [ -z "${group_state[$group]}" ] && group_state[$group]="UNBUILT"
        continue
    fi
    if [ ! -d "$here/$dir" ]; then
        echo "│   $tee── $label $dots ABSENT (the frontier closes here)"
        frontier_open=0
        [ -z "${group_state[$group]}" ] && group_state[$group]="UNBUILT"
        continue
    fi

    make -C "$here/$dir" clean >/dev/null 2>&1 || true
    if ! make -C "$here/$dir" >/dev/null 2>&1; then
        echo "│   $tee── $label $dots FAIL (build)"
        failed=1
        group_state[$group]="FAIL"
        continue
    fi

    valid=1
    ( cd "$here/$dir" && ./run >run.out 2>&1 ) || valid=0
    if [ -x "$here/$dir/check_refusals.sh" ]; then
        ( cd "$here/$dir" && ./check_refusals.sh >>run.out 2>&1 ) || valid=0
    fi

    if [ "$valid" -eq 1 ]; then
        echo "│   $tee── $label $dots PASS"
        [ -z "${group_state[$group]}" ] && group_state[$group]="PASS"
    else
        echo "│   $tee── $label $dots FAIL"
        sed 's/^/│       /' "$here/$dir/run.out"
        failed=1
        group_state[$group]="FAIL"
    fi
done

echo "│"
if [ "$failed" -ne 0 ]; then
    echo "└── the ladder stops at the first failure"
    exit 1
fi
echo "├── Gate A · structure ............... ${group_state[A]}"
echo "├── Gate B · solve + constitution .... ${group_state[B]}"
echo "├── Gate C · statement ............... ${group_state[C]}"

if [ "${group_state[C]}" != "PASS" ]; then
    echo "└── total sensitivity .............. (unbuilt)"
    exit 0
fi

out="$here/level-9-statement/run.out"
marks=$(grep -c 'ADJOINT_RESULT =' "$out")
values=$(grep -o 'ADJOINT_RESULT =.*' "$out" | sed 's/.*ADJOINT_RESULT =//')
if ! marker_valid "$marks" 1 "$values"; then
    echo "└── RUNNER FAILURE: the statement did not report one sensitivity"
    exit 1
fi
echo "└── total sensitivity df/dp ......... $(echo $values)"
