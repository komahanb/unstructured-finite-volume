#!/bin/bash
# build the library and the nonsymmetric-krylov suite, run it.
# (pure serial - the test operator is a synthetic advection-diffusion csr,
# no mesh / coarrays needed)
set -e

root="$(cd "$(dirname "$0")/../.." && pwd)"
here="$root/test/krylov"
flags="-std=f2018 -fcoarray=single -cpp -fPIC -Wno-line-truncation"

( cd "$root" && ./build.sh >/dev/null )

cd "$here"
gfortran $flags -I "$root/lib/" -I "$root/src/" -c "$root/test/fixtures/class_csr_system.f90" -o csr_system.o
gfortran $flags -I "$root/lib/" -I . -c test.f90 -o test.o
gfortran csr_system.o test.o -o run "$root/lib/libufvm.a" -fcoarray=single

./run
