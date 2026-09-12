#!/bin/bash
# build the library and the algebra suite, run the laws, then run
# every refusal and assert each terminates for its stated reason.
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null )
fi

make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null

cd "$suite_dir" && ./run

declare -A reason=(
  [slot]="a slot index must name a slot of the relation"
  [embed]="a restriction domain must embed in the slot it restricts"
  [none]="a projection selects at least one slot"
  [range]="a slot index must name a slot of the relation"
  [repeat]="a projection selects each slot at most once"
  [binary]="composition takes two binary relations"
  [middle]="composition requires one shared middle domain"
)

for case in slot embed none range repeat binary middle; do
    if ./refusal "$case" >refusal.out 2>&1; then
        echo " FAIL : '$case' was accepted"
        exit 1
    fi
    if grep -q "${reason[$case]}" refusal.out; then
        echo " PASS : '$case' is refused, loudly"
    else
        echo " FAIL : '$case' terminated for the wrong reason"
        cat refusal.out
        exit 1
    fi
done
rm -f refusal.out
