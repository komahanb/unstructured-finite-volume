#!/bin/bash
# build the library and the relation suite, run the laws, then run
# every refusal and assert that each dies for its stated reason.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

declare -A reason=(
  [member]="a tuple names a member its domain does not hold"
  [arity]="each tuple has exactly one part per slot"
  [undeclared]="a signature refers to declared domains only"
  [empty]="a relation relates at least one domain"
  [twice]="a relation never signs twice"
)

for case in member arity undeclared empty twice; do
    if ./refusal "$case" >refusal.out 2>&1; then
        echo " FAIL : '$case' was accepted"
        exit 1
    fi
    if grep -q "${reason[$case]}" refusal.out; then
        echo " PASS : '$case' is refused, loudly"
    else
        echo " FAIL : '$case' died for the wrong reason"
        cat refusal.out
        exit 1
    fi
done
rm -f refusal.out
