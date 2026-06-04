#!/bin/bash
# build the library and the AMG verification suite, run it.
# (generates the meshes the suite needs via meshgen)
set -e

root="$(cd "$(dirname "$0")/../.." && pwd)"
here="$root/test/amg"
flags="-std=f2018 -fcoarray=single -cpp -fPIC -Wno-line-truncation"

( cd "$root" && ./build.sh >/dev/null )

for m in box-36 sphere square-10 square-20 square-40 square-80; do
   bash "$root/meshgen/ensure.sh" "$here/$m.msh"
done

cd "$here"
gfortran $flags -I "$root/lib/" -c test.f90 -o test.o
gfortran test.o -o run "$root/lib/libufvm.a" -fcoarray=single

./run
