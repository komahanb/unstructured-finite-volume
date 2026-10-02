#!/bin/bash
# Build and run the programs of the Taylor-Green vortex:
#
#   taylor_green_vortex  the Navier-Stokes Taylor-Green vortex on the
#                        periodic box of box.geo, 16 x 16 cells, ten
#                        instants, the chain dirk(2), bdf(2), adams(2);
#                        the error against the exact solution, second
#                        order in the cell width (8.93e-2 at 8 x 8,
#                        2.33e-2 at 16 x 16)
#   taylor_green_vortex_3d
#                        the vortex extended along z on the periodic
#                        8 x 8 x 8 box of box_3d.py, five instants,
#                        finite volumes of order 2
#
# gmsh meshes box.geo into box.msh, and box_3d.py writes box_3d.msh,
# when the mesh is absent.
#
#     run.sh             on one process
#     IMAGES=n run.sh    on n coarray images: the library and the
#                        programs built with COARRAY=lib, each program
#                        run by mpirun -np n (MPIRUN_FLAGS adds flags)
set -e

dir="$(cd "$(dirname "$0")" && pwd)"
images=${IMAGES:-0}

if [ "$images" -eq 0 ]; then
    if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
        ( cd "$dir/../.." && ./build.sh >/dev/null )
    fi
    make -C "$dir" clean >/dev/null 2>&1 || true
    make -C "$dir" >/dev/null
else
    # the coarray library is rebuilt by the Makefile when a source of
    # src is newer than it
    make -C "$dir" clean >/dev/null 2>&1 || true
    make -C "$dir" COARRAY=lib >/dev/null
fi

make -C "$dir" box.msh box_3d.msh >/dev/null

cd "$dir"
launch() {
    echo "== $*"
    if [ "$images" -eq 0 ]; then
        "./$1" "${@:2}"
    elif [ "$images" -eq 1 ]; then
        OMPI_MCA_osc=pt2pt mpirun ${MPIRUN_FLAGS:-} -np 1 "./$1" "${@:2}"
    else
        mpirun ${MPIRUN_FLAGS:-} -np "$images" "./$1" "${@:2}"
    fi
}

launch taylor_green_vortex
launch taylor_green_vortex_3d
