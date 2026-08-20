#!/bin/bash
# build the library and the traversal benchmark, run it.
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run

# The section-66 reference: printed beside every run so a regression
# reads as a before/after table.
echo ""
echo " --- baseline (test/graph-benchmark/baseline) ---"
cat "$here/baseline"
