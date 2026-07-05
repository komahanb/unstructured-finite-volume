#!/bin/bash
# build the library and the standalone graph suite, run it.
# (no mesh needed - the graph is exercised by itself)
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run
