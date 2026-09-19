#!/bin/bash
# Build the library and run the three programs of this suite:
#
#   dirk_coefficients    the two-stage third-order SDIRK tableau as the
#                        zero of a continuous_residual on a point
#                        manifold; W11 = h (3 + sqrt 3)/6 = 0.078867...
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
#   run                  the same tableau on the expression and the
#                        residual_operator directly
#
# gmsh meshes box.geo into box.msh when the mesh is absent.
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null )
fi

make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null
make -C "$suite_dir" box.msh box_3d.msh >/dev/null

cd "$suite_dir"
echo "== dirk_coefficients"
./dirk_coefficients
echo "== taylor_green_vortex"
./taylor_green_vortex
echo "== taylor_green_vortex_3d"
./taylor_green_vortex_3d
echo "== run"
./run
