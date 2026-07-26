#!/bin/bash
# build the library and the mandelbrot-physics suite, run it.
# (the fractal is painted on the same mesh the orbit painting used)
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

bash "$here/../../meshgen/ensure.sh" "$here/../square-tri-40.msh"

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run
