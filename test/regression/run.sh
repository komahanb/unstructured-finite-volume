#!/bin/bash
# Build the library and the regression suites (3d and 2d), then run them.
set -e

root="$(cd "$(dirname "$0")/../.." && pwd)"
flags="-std=f2018 -fcoarray=single -cpp -fPIC -Wno-line-truncation"

make -C "$root/src" F90=gfortran         >/dev/null
make -C "$root/src" install F90=gfortran >/dev/null

cd "$root/test/regression"

# generate the meshes the suites need (msh 4.1, on the fly - nothing committed)
for m in sphere box-3 box-36 square-10; do
   bash "$root/meshgen/ensure.sh" "$root/test/$m.msh"
done

gfortran $flags -I "$root/lib/" -c test.f90 -o test.o
gfortran test.o -o run "$root/lib/libufvm.a" -llapack -fcoarray=single

gfortran $flags -I "$root/lib/" -c test_2d.f90 -o test_2d.o
gfortran test_2d.o -o run2d "$root/lib/libufvm.a" -llapack -fcoarray=single

echo "===================== 3D suite ====================="
./run

echo "===================== 2D suite ====================="
./run2d

echo "================= order of accuracy ================"
bash "$root/test/spatial-discretization/run.sh"
