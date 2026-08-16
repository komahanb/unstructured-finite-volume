#!/bin/bash
# Build the library and the sequence-view suite; compare the two
# candidate API forms; run the laws and every refusal.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

F90="$(make -C "$here" -s print-f90)"
FSTD="$(make -C "$here" -s print-std)"

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

#---------------------------------------------------------------------
# Section 3. Two candidate representations of a sequence: as a branch,
# and as the first cell. Both compile and agree; the branch form is
# smaller and needs no artificial empty graph. Kept so the rejected
# form stays measurable rather than remembered.
#---------------------------------------------------------------------

echo " CANDIDATE API FORMS"
cd "$here/candidates"
rm -f ./*.mod ./*.o compare
$F90 -std=$FSTD -fcoarray=single -I"$here/../../lib" \
     branch_form.f90 graph_form.f90 compare.f90 \
     "$here/../../lib/libufvm.a" -o compare
./compare
printf '   branch form : %s module code lines, 1 line per call site\n' \
       "$($F90 --version >/dev/null; grep -c -v -E '^[[:space:]]*(!|$)' branch_form.f90)"
printf '   graph  form : %s module code lines, 5 lines per call site\n' \
       "$(grep -c -v -E '^[[:space:]]*(!|$)' graph_form.f90)"
echo " PASS : the branch form is the smaller exact one; it is what src ships"
rm -f ./*.mod ./*.o compare
cd "$here"
echo ''

#---------------------------------------------------------------------

./run

declare -A message=(
  [cellnull]="a sequence cell holds a KNOWN element"
  [cellunknown]="a sequence cell holds a KNOWN element"
  [sizeunknownholder]="the extent depends on an unknown tail"
  [sizeunknowntail]="the extent depends on an unknown tail"
  [containsunknowntail]="membership depends on an unknown tail"
  [indexzero]="a sequence is indexed from one"
  [pastend]="the sequence has no such element"
  [pastunknown]="that element lies beyond an unknown tail"
  [emptyindexed]="the sequence has no such element"
)

echo ''
for case in cellnull cellunknown sizeunknownholder sizeunknowntail \
            containsunknowntail indexzero pastend pastunknown emptyindexed; do
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
