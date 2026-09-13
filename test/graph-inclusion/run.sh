#!/bin/bash
# Build the library and the inclusion suite, run the laws, then run
# every refusal and assert each stops for its stated reason.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$here/../.." && ./build.sh >/dev/null )
fi

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

declare -A reason=(
  [unsigned]="requires the added set to have an assigned identity"
  [selfsame]="a set cannot be declared into itself"
  [twoambients]="it cannot be declared twice"
  [cycle]="so it revisits a set - a cycle"
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
