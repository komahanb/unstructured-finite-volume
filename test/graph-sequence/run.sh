#!/bin/bash
# Build the library and the sequence-view suite; compare the two
# candidate API forms; run the laws and every refusal.
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
# Section 3. Two candidate representations of a sequence: as a branch,
# and as the first cell. Both compile and agree; the branch form is
# smaller and needs no artificial empty graph. Retained so the rejected
# form stays measurable rather than remembered.
#---------------------------------------------------------------------

echo " CANDIDATE API FORMS"
cd "$suite_dir/candidates"
rm -f ./*.mod ./*.o compare
$F90 -std=$FSTD -fcoarray=single -I"$suite_dir/../../lib" \
     branch_form.f90 graph_form.f90 compare.f90 \
     "$suite_dir/../../lib/libufvm.a" -o compare
./compare
printf '   branch form : %s module code lines, 1 line per call site\n' \
       "$($F90 --version >/dev/null; grep -c -v -E '^[[:space:]]*(!|$)' branch_form.f90)"
printf '   graph  form : %s module code lines, 5 lines per call site\n' \
       "$(grep -c -v -E '^[[:space:]]*(!|$)' graph_form.f90)"
echo " PASS : the branch form is the smaller exact one; it is what src ships"
rm -f ./*.mod ./*.o compare
cd "$suite_dir"
echo ''

#---------------------------------------------------------------------

./run

declare -A message=(
  [cellnull]="a sequence cell contains a KNOWN element"
  [cellunknown]="a sequence cell contains a KNOWN element"
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
