#!/bin/bash
# build the library and the binary relation suite, run the laws, then
# run every refusal and assert that each stops for its stated reason.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$here/../.." && ./build.sh >/dev/null )
fi

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

declare -A reason=(
  [member]="a tuple names a member its domain does not contain"
  [arity]="each tuple has exactly one part per position"
  [undeclared]="a signature refers to declared domains only"
)

for case in member arity undeclared; do
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
