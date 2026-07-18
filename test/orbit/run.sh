#!/bin/bash
# build the library and the orbit painting, run it.
# (needs the square mesh - the painting is quantized onto its cells)
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

# the canvas: the suite reads ../square-80.msh, which is generated, not tracked
bash "$here/../../meshgen/ensure.sh" "$here/../square-80.msh"

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run
