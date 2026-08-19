#!/bin/bash
# build the library and the graph mesh suite, run it.
# # rung 2: the mesh on the tower
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

# Quarantine guard: the pre-tower graph stack (interface_graph,
# class_stored_graph, class_mesh) may be imported only by itself and
# by the mesh builder bridge, which exists to retire it. Any new
# importer fails here.
src="$here/../../src"
allowed="class_stored_graph.f90 class_mesh.f90 class_mesh_builder.f90"
offenders=""
for f in "$src"/*.f90; do
    base=$(basename "$f")
    case " $allowed " in *" $base "*) continue;; esac
    grep -q "use interface_graph\|use class_stored_graph\|use class_mesh\b" "$f" \
        && offenders="$offenders $base"
done
if [ -n "$offenders" ]; then
    echo " FAIL : the pre-tower graph stack gained an importer:$offenders"
    exit 1
else
    echo " PASS : the pre-tower graph stack has no importer outside its quarantine"
fi
