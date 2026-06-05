#!/bin/bash
# build the library and the advection verification suite, run it.
# (generates the square meshes the suite needs via meshgen)
set -e

root="$(cd "$(dirname "$0")/../.." && pwd)"
here="$root/test/advection"
flags="-std=f2018 -fcoarray=single -cpp -fPIC -Wno-line-truncation"

( cd "$root" && ./build.sh >/dev/null )

for m in square-10 square-20 square-40; do
   bash "$root/meshgen/ensure.sh" "$here/$m.msh"
done

cd "$here"
gfortran $flags -I "$root/lib/" -c test.f90 -o test.o
gfortran test.o -o run "$root/lib/libufvm.a" -fcoarray=single

./run
