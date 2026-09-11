#!/bin/bash
# build the library and the graph mesh suite, run it.
# # rung 2: the mesh on the tower
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null )
fi

# the gmsh path reads test/square-10.msh, generated from its recipe.
[ -e "$suite_dir/../square-10.msh" ] || "$suite_dir/../../meshgen/ensure.sh" "$suite_dir/../square-10.msh" >/dev/null

make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null

cd "$suite_dir" && ./run

# The pre-tower graph stack was deleted after view_mesh_geometry
# reproduced its measurements bitwise on all eleven sample meshes.
# The check: the files must not return and nothing may reference them.
src="$suite_dir/../../src"
for deleted_source in interface_graph.f90 class_stored_graph.f90 class_mesh.f90 \
            class_array_mesh_loader.f90; do
    if [ -e "$src/$deleted_source" ]; then
        echo " FAIL : the deleted pre-tower file $deleted_source has returned"
        exit 1
    fi
done
grep -q "use interface_graph\|use class_stored_graph\|use class_mesh \|use class_mesh," "$src"/*.f90 \
    && { echo " FAIL : a reference to the deleted pre-tower stack exists"; exit 1; } \
    || echo " PASS : the pre-tower graph stack stays deleted"
