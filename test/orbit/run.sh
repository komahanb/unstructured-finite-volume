#!/bin/bash
# build the library and the orbit painting, run it.
# (needs the square mesh - the painting is quantized onto its cells)
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run
