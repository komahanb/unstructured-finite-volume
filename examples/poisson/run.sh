#!/bin/bash
# build, generate the square mesh on the fly (msh 4.1, gmsh), and solve the
# poisson convergence driver (uses square-40 by default).
set -e

here="$(cd "$(dirname "$0")" && pwd)"

# clean rebuild so the driver relinks against the current libufvm
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

bash "$here/../../meshgen/ensure.sh" "$here/square-40.msh"

cd "$here" && ./run
