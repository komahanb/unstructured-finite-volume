#!/bin/bash
# build the library and the graph minimization suite, run it.
# # rung 2: the mesh on the tower
set -e

here="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$here/../.." && ./build.sh >/dev/null )
fi

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

if ./refusal >refusal.out 2>&1; then
    echo " FAIL : a lopsided system was accepted"
    exit 1
fi
grep -q "equal" refusal.out && echo " PASS : the square family refuses a lopsided system" || { cat refusal.out; exit 1; }
rm -f refusal.out

# solver restriction refusals: one invalid selection per case
declare -A reason=(
    [outside]='a restriction selects unknowns of the whole domain'
    [repeated]='a restriction selects each unknown once'
    [empty]='a restriction selects an unknown at least'
    [split_block]='a restriction selects whole blocks'
    [aggregates]='a restriction selects unknowns of the stated aggregates'
    [flags]='a restriction selects unknowns of the stated flags'
    [no_retained]='a restriction retains an unknown at least'
    [partition]='a restriction selects unknowns of the partition'
    [unstated]='the operator is stated on the solver domain before it is evaluated'
    [state_domain]='the state is defined on the unknown domain'
    [design_domain]='the design is defined on the point domain with one value per point'
    [design_count]='the design is defined on the point domain with one value per point'
    [state_count]='a value vector must fill its domain exactly'
    [pairing_domain]='an inner product pairs fields on the same domain'
    [host]='the residual is applied on its own unknown graph'
    [direction_domain]='a direction in the state is defined on the unknown domain'
    [stored_domain]='the design is defined on the point domain with one value per point'
)
for case_name in "${!reason[@]}"; do
    if ./restriction_refusal "$case_name" >refusal.out 2>&1; then
        echo " FAIL : an invalid restriction was admitted: $case_name"
        exit 1
    fi
    if ! grep -Fq "${reason[$case_name]}" refusal.out; then
        cat refusal.out
        echo " FAIL : an unrelated refusal: $case_name"
        exit 1
    fi
    echo " PASS : an invalid restriction is refused: $case_name"
done
rm -f refusal.out
