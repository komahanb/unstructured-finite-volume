#!/bin/bash
# build the library and the graph minimization suite, run it.
# # rung 2: the mesh on the tower
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

if ./refusal >refusal.out 2>&1; then
    echo " FAIL : a lopsided system was accepted"
    exit 1
fi
grep -q "equal" refusal.out && echo " PASS : the square family refuses a lopsided system" || { cat refusal.out; exit 1; }
rm -f refusal.out
