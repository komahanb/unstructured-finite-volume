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
    [outside]='every selected index must belong to the whole domain'
    [repeated]='a restriction must select each unknown once'
    [empty]='a restriction must select an unknown at least'
    [split_block]='a restriction must select whole blocks'
    [aggregates]='a restriction must select unknowns of the stated aggregates'
    [flags]='a restriction must select unknowns of the stated flags'
    [no_retained]='a restriction must retain an unknown at least'
    [partition]='a restriction must select unknowns of the partition'
    [unstated]='evaluate was called before the operator was stated on the solver domain'
    [state_domain]='the bound state must be defined on the unknown domain'
    [design_domain]='the design must be defined on the point domain with one value per point'
    [design_count]='the design must be defined on the point domain with one value per point'
    [state_count]='requires the values to fill num_entries * num_components exactly'
    [pairing_domain]='requires fields on the same domain'
    [other_placement]='the design must be defined on the point domain with one value per point'
    [host]="received a graph that is not this residual's own unknown graph"
    [direction_domain]='a direction in the state must be defined on the unknown domain'
    [stored_domain]='the design must be defined on the point domain with one value per point'
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
