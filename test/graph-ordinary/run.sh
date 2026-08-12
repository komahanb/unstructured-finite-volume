#!/bin/bash
# build the library and the ordinary profile suite, run the laws,
# then run every refusal and assert each dies for its stated reason.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

declare -A reason=(
  [tailless]="every edge has exactly one tail"
  [twotailed]="every edge has exactly one tail"
  [twoheaded]="no edge has two heads"
  [ternary]="must be binary"
  [mismatched]="must share the tail's domains"
)

for case in tailless twotailed twoheaded ternary mismatched; do
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
