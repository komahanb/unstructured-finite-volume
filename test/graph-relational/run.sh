#!/bin/bash
# Build the library and the relational-view suite; compile the
# assignment-prohibition candidates; run the laws, the lifetime cases
# and every refusal.
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null )
fi

F90="$(make -C "$suite_dir" -s print-f90)"
FSTD="$(make -C "$suite_dir" -s print-std)"

make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null

#---------------------------------------------------------------------
# Can Fortran prohibit assignment of a derived type from client code?
# Each candidate is either rejected by the compiler, or built and run
# to show what it fails to prevent.
#---------------------------------------------------------------------

reject () {
    if $F90 -std=$FSTD -fcoarray=single -Wall -pedantic -fsyntax-only "$1.f90" 2>compiler.out; then
        echo " FAIL : $1 compiled; the reported constraint does not hold"
        exit 1
    fi
    if grep -q "$2" compiler.out; then
        echo " PASS : $1 is rejected, for the reported reason"
    else
        echo " FAIL : $1 is rejected for a different reason"
        cat compiler.out; exit 1
    fi
}

execute () {
    if ! $F90 -std=$FSTD -fcoarray=single -Wall -pedantic "$1.f90" -o candidate 2>compiler.out; then
        echo " FAIL : $1 does not build"; cat compiler.out; exit 1
    fi
    ./candidate
    echo " PASS : $1 admits the assignment it meant to prohibit"
}

echo " ASSIGNMENT-PROHIBITION CANDIDATES ($F90 -std=$FSTD)"
cd "$suite_dir/fortran-assignment"
rm -f ./*.mod ./*.o candidate

execute coarray_component
execute private_generic
execute unmatched_specific
reject  coarray_pointer "shall be a nonpointer, nonallocatable scalar"
reject  lock_component   "must have a codimension or be a subcomponent of a coarray"

# The rule the refusal itself rests on: an INTENT(OUT) dummy of a
# finalizable type is finalized before the body runs.
if ! $F90 -std=$FSTD -fcoarray=single "intent_out_finalizes.f90" -o candidate 2>compiler.out; then
    echo " FAIL : intent_out_finalizes does not build"; cat compiler.out; exit 1
fi
./candidate
echo " PASS : intent_out_finalizes shows why the refusal takes INTENT(INOUT)"

rm -f compiler.out ./*.mod ./*.o candidate
cd "$suite_dir"
echo ''

#---------------------------------------------------------------------

./run
echo ''
./lifetime

declare -A message=(
  [nosetbound]="no member set is bound to that element"
  [norelationbound]="no relation is bound to that element"
  [sharedbind]="a binding is extended by its sole owner"
  [releasedtwin]="this binding's objects have been released"
  [unsignedset]="a binding stores identified objects"
  [unsignedrelation]="a binding stores identified objects"
  [boundview]="a view cannot be bound"
)

echo ''
for case in nosetbound norelationbound sharedbind releasedtwin \
            unsignedset unsignedrelation boundview; do
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
