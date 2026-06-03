#!/bin/bash
# build, generate the meshes on the fly (msh 4.1, gmsh), and run the
# loader test (loads box-36 + cylinder-coarse, reads rectangle).
set -e

here="$(cd "$(dirname "$0")" && pwd)"

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

for m in box-36 cylinder-coarse rectangle; do
   bash "$here/../../meshgen/ensure.sh" "$here/../$m.msh"
done

cd "$here" && ./run
