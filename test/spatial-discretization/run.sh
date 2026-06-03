#!/bin/bash
# build, generate the refinement meshes on the fly (msh 4.1, gmsh), and
# run the spatial order-of-accuracy study (2d poisson, square-10..80).
set -e

here="$(cd "$(dirname "$0")" && pwd)"

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

for n in 10 20 40 80; do
   bash "$here/../../meshgen/ensure.sh" "$here/../square-$n.msh"
done

cd "$here" && ./run
