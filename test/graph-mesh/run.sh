#!/bin/bash
# build the library and the graph mesh suite, run it.
# # rung 2: the mesh on the tower
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

# The pre-tower graph stack was deleted after view_mesh_geometry
# reproduced its measurements bitwise on all eleven sample meshes.
# This guard keeps it deleted: the files must not return and nothing
# may reference them.
src="$here/../../src"
for dead in interface_graph.f90 class_stored_graph.f90 class_mesh.f90 \
            class_array_mesh_loader.f90; do
    if [ -e "$src/$dead" ]; then
        echo " FAIL : the deleted pre-tower file $dead has returned"
        exit 1
    fi
done
grep -q "use interface_graph\|use class_stored_graph\|use class_mesh \|use class_mesh," "$src"/*.f90 \
    && { echo " FAIL : a reference to the deleted pre-tower stack exists"; exit 1; } \
    || echo " PASS : the pre-tower graph stack stays deleted"
