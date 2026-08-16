#!/bin/bash
# Build the library and the identity-map lifetime suite, run the law
# natively, then run it again under valgrind and assert that a map
# outliving its binders reads nothing it does not own.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./lifetime

echo ''
if command -v valgrind >/dev/null 2>&1; then
    # The law is about STORAGE, so the allocator must agree with the
    # native run. Before the gate this same probe reported 170 invalid
    # reads from 48 contexts and still printed ALL PROPOSITIONS HOLD.
    if valgrind --error-exitcode=9 --leak-check=full \
                --errors-for-leak-kinds=definite \
                ./lifetime >valgrind.out 2>&1; then
        echo " PASS : no map reads storage it does not own"
    else
        echo " FAIL : the maps still borrow their keys"
        grep -E "ERROR SUMMARY|Invalid read|definitely lost" valgrind.out || true
        exit 1
    fi
    rm -f valgrind.out
else
    echo " SKIP : valgrind is not installed; the native run stands alone"
fi
