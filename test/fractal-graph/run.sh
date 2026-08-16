#!/bin/bash
# Build the library and the spike, compile the declaration-order
# probes, run the test suite, run every refusal, and measure the
# kernel.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

#---------------------------------------------------------------------
# Declaration order. Two representations must compile and one must be
# rejected, with the same compiler and standard the kernel is built
# with.
#---------------------------------------------------------------------

F90="$(make -C "$here" -s print-f90)"
FSTD="$(make -C "$here" -s print-std)"

echo " DECLARATION ORDER ($F90 -std=$FSTD)"
cd "$here/fortran-recursion"
rm -f ./*.mod
for admitted in branch_before_graph polymorphic_known; do
    if $F90 -std=$FSTD -fsyntax-only $admitted.f90 2>probe.out; then
        echo " PASS : $admitted compiles"
    else
        echo " FAIL : $admitted does not compile"; cat probe.out; exit 1
    fi
done
if $F90 -std=$FSTD -fsyntax-only graph_before_branch.f90 2>probe.out; then
    echo " FAIL : graph_before_branch compiled; the order constraint is not as reported"
    exit 1
fi
if grep -q "has not been previously defined" probe.out; then
    echo " PASS : graph_before_branch is rejected, for the reported reason"
else
    echo " FAIL : graph_before_branch is rejected for a different reason"
    cat probe.out; exit 1
fi
rm -f probe.out ./*.mod
cd "$here"
echo ''

#---------------------------------------------------------------------

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

./run

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
echo " KERNEL SIZE (fractal_graph.f90)"
printf '   code    : %s lines\n' "$(grep -c -v -E '^[[:space:]]*(!|$)' fractal_graph.f90)"
printf '   comment : %s lines\n' "$(grep -c -E '^[[:space:]]*!'        fractal_graph.f90)"
printf '   blank   : %s lines\n' "$(grep -c -E '^[[:space:]]*$'        fractal_graph.f90)"
