#!/bin/bash
# Build the library and the relational-view suite; compile the
# assignment-prohibition candidates; run the laws, the lifetime cases
# and every refusal.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

F90="$(make -C "$here" -s print-f90)"
FSTD="$(make -C "$here" -s print-std)"

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

#---------------------------------------------------------------------
# Can Fortran prohibit assignment of a derived type from client code?
# Each candidate is either rejected by the compiler, or built and run
# to show what it fails to prevent.
#---------------------------------------------------------------------

reject () {
    if $F90 -std=$FSTD -fcoarray=single -Wall -pedantic -fsyntax-only "$1.f90" 2>probe.out; then
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
    if ! $F90 -std=$FSTD -fcoarray=single -Wall -pedantic "$1.f90" -o candidate 2>probe.out; then
        echo " FAIL : $1 does not build"; cat probe.out; exit 1
    fi
    ./candidate
    echo " PASS : $1 admits the assignment it meant to prohibit"
}

echo " ASSIGNMENT-PROHIBITION CANDIDATES ($F90 -std=$FSTD)"
cd "$here/fortran-assignment"
rm -f ./*.mod ./*.o candidate

execute coarray_component
execute private_generic
execute unmatched_specific
reject  coarray_pointer "shall be a nonpointer, nonallocatable scalar"
reject  lock_component   "must have a codimension or be a subcomponent of a coarray"

# The rule the refusal itself rests on: an INTENT(OUT) dummy of a
# finalizable type is finalized before the body runs.
if ! $F90 -std=$FSTD -fcoarray=single "intent_out_finalizes.f90" -o candidate 2>probe.out; then
    echo " FAIL : intent_out_finalizes does not build"; cat probe.out; exit 1
fi
./candidate
echo " PASS : intent_out_finalizes shows why the refusal takes INTENT(INOUT)"

rm -f probe.out ./*.mod ./*.o candidate
cd "$here"
echo ''

#---------------------------------------------------------------------

./run
echo ''
./lifetime

declare -A message=(
  [nosetbound]="no member set is bound to that element"
  [norelationbound]="no relation is bound to that element"
  [assign]="a relational_binding is not assignable"
)

echo ''
for case in nosetbound norelationbound assign; do
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
