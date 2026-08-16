#!/bin/bash
# Build the library and the map prototypes, then run them. Analysis
# evidence for doc/member-set-fractal-map.md - not production laws.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here"
./set
echo ''
./scale
