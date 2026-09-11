#!/bin/bash
# build the library and the epistemic-view suite, run the laws, then
# run every refusal and assert that each terminates for its stated reason.
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null )
fi

make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null

cd "$suite_dir" && ./run

declare -A reason=(
  [dataunknown]="Q is not KNOWN"
  [residualunknown]="R is not KNOWN"
  [nullname]="NULL has no epistemic name"
)

for case in dataunknown residualunknown nullname; do
    if ./refusal "$case" >refusal.out 2>&1; then
        echo " FAIL : '$case' was accepted"
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
