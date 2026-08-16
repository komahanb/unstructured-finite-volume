#!/bin/bash
# Build the library and the spike; compile the language fixtures and
# the immutability candidates; run the law suite, the mutation
# experiments and every refusal; measure the kernel.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

F90="$(make -C "$here" -s print-f90)"
FSTD="$(make -C "$here" -s print-std)"
KERNEL="-fcoarray=single -I$here -I$here/../../lib"
LINK="$here/../../lib/libufvm.a"

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

admit () {
    if $F90 -std=$FSTD -I"$here" -I"$here/../../lib" -fsyntax-only "$1.f90" 2>probe.out; then
        echo " PASS : $1 compiles"
    else
        echo " FAIL : $1 does not compile"; cat probe.out; exit 1
    fi
}

reject () {
    if $F90 -std=$FSTD -I"$here" -I"$here/../../lib" -fsyntax-only "$1.f90" 2>probe.out; then
        echo " FAIL : $1 compiled; the reported constraint does not hold"
        exit 1
    fi
    if grep -q "$2" probe.out; then
        echo " PASS : $1 is rejected, for the reported reason"
    else
        echo " FAIL : $1 is rejected for a different reason"
        cat probe.out; exit 1
    fi
}

execute () {
    if ! $F90 -std=$FSTD $KERNEL "$1.f90" "${@:2}" $LINK -o candidate 2>probe.out; then
        echo " FAIL : $1 does not build"; cat probe.out; exit 1
    fi
    ./candidate
    echo " PASS : $1 builds and runs"
}

#---------------------------------------------------------------------
# Declaration order, navigation, and the branch invariant. Compiled
# against the production kernel in src, not a copy of it.
#---------------------------------------------------------------------

echo " FIXTURES ($F90 -std=$FSTD)"
cd "$here/fortran-recursion"
rm -f ./*.mod

admit branch_before_graph
admit polymorphic_known
admit encapsulated_navigation
reject graph_before_branch          "has not been previously defined"
reject function_result_base         "leftmost part-ref in a data-ref cannot be a function reference"
reject chained_navigation           "Unclassifiable statement"
reject private_status_write         "status_.* is a PRIVATE component"
reject private_reference_write      "known_.* is a PRIVATE component"
reject private_structure_constructor "status_.* is a PRIVATE component"

rm -f probe.out ./*.mod
echo ''

#---------------------------------------------------------------------
# Candidate mechanisms for publication immutability. Each is either
# rejected by the compiler, or built and run to show what it fails to
# prevent.
#---------------------------------------------------------------------

echo " IMMUTABILITY CANDIDATES ($F90 -std=$FSTD)"
cd "$here/fortran-mutation"
rm -f ./*.mod ./*.o candidate

$F90 -std=$FSTD -c accessor_candidates.f90 -o accessor_candidates.o

execute seal_parent_assignment
execute seal_branch_assignment
execute pointer_accessor_assignment accessor_candidates.o
execute transitive_mutation
reject  value_accessor_assignment "must have the pointer attribute"
reject  accessor_navigation       "Invalid character in name"

rm -f probe.out ./*.mod ./*.o candidate
cd "$here"
echo ''

#---------------------------------------------------------------------

./run
echo ''
./mutation

declare -A message=(
  [twice]="graph identity is assigned once"
  [undeclared]="KNOWN requires a graph with assigned identity"
  [nonumber]="no number bound to this identity"
  [nosymbol]="no symbol bound to this identity"
  [noindex]="no index bound to this identity"
  [nooperator]="no operator defined for this symbol"
  [unknownvalue]="UNKNOWN branch has no value"
  [unknownmember]="UNKNOWN member is not NULL"
)

echo ''
for case in twice undeclared nonumber nosymbol noindex nooperator \
            unknownvalue unknownmember; do
    if ./refusal "$case" >refusal.out 2>&1; then
        echo " FAIL : $case was admitted"
        exit 1
    fi
    if grep -q "${message[$case]}" refusal.out; then
        echo " PASS : $case is refused"
    else
        echo " FAIL : $case is refused for a different reason"
        cat refusal.out
        exit 1
    fi
done
rm -f refusal.out

echo ''
echo " KERNEL SIZE (src/fractal_graph.f90)"
printf '   code    : %s lines\n' "$(grep -c -v -E '^[[:space:]]*(!|$)' ../../src/fractal_graph.f90)"
printf '   comment : %s lines\n' "$(grep -c -E '^[[:space:]]*!'        ../../src/fractal_graph.f90)"
printf '   blank   : %s lines\n' "$(grep -c -E '^[[:space:]]*$'        ../../src/fractal_graph.f90)"
