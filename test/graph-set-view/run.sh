#!/bin/bash
# Build the library and the set foundation suite, run the laws, then
# run every refusal and assert each dies for its stated reason.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

declare -A reason=(
  [unsigned]="keyed on assigned identity"
  [twice]="a set is described once"
  [undescribed]="no representation describes that set"
)

echo ''
for case in unsigned twice undescribed; do
    if ./refusal "$case" >refusal.out 2>&1; then
        echo " FAIL : '$case' was admitted"
        exit 1
    fi
    if grep -q "${reason[$case]}" refusal.out; then
        echo " PASS : '$case' is refused"
    else
        echo " FAIL : '$case' is refused for a different reason"
        cat refusal.out
        exit 1
    fi
done
rm -f refusal.out
