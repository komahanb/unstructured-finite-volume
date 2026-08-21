#!/bin/bash
# Build the library and the inclusion suite, run the laws, then run
# every refusal and assert each dies for its stated reason.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

declare -A reason=(
  [unsigned]="keyed on assigned identity"
  [selfsame]="not declared into itself"
  [twoambients]="declared into one ambient"
  [cycle]="an inclusion chain is finite"
)

echo ''
for case in unsigned selfsame twoambients cycle; do
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
