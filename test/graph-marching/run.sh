#!/bin/bash
# build the library and the graph marching suite, run it.
# # rung 2: the mesh on the tower
set -e

here="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$here/../.." && ./build.sh >/dev/null )
fi

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run
