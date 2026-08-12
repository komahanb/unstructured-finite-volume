#!/bin/bash
# build the library and the carrier suite, run the laws, then run the
# refusal and assert that it dies for the stated reason.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

if ./refusal >refusal.out 2>&1; then
    echo " FAIL : a second declare was accepted"
    exit 1
fi
if grep -q "a domain never signs twice" refusal.out; then
    echo " PASS : a second declare is refused, loudly"
else
    echo " FAIL : the refusal died for the wrong reason"
    cat refusal.out
    exit 1
fi
rm -f refusal.out
