#!/bin/bash
# build, generate the mesh on the fly (msh 4.1, gmsh), and run the
# unsteady/2d solve test (rectangle).
set -e

here="$(cd "$(dirname "$0")" && pwd)"

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

bash "$here/../../meshgen/ensure.sh" "$here/../rectangle.msh"

cd "$here" && ./run
