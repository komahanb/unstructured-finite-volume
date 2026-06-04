#!/bin/bash
# build the library and the BDF integrator order-of-accuracy test, run it.
# (no mesh needed - the test marches a scalar ode)
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run
