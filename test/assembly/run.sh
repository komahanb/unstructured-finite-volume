#!/bin/bash
# build, generate the mesh on the fly (msh 4.1, gmsh), and run the
# assembly test (A = L + U + D on box-36).
set -e

here="$(cd "$(dirname "$0")" && pwd)"

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

bash "$here/../../meshgen/ensure.sh" "$here/../box-36.msh"

cd "$here" && ./run
