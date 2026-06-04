#!/bin/bash
# build the library and the adjoint gradient verifications (steady +
# transient), run them. (generates the small box mesh on the fly)
set -e

root="$(cd "$(dirname "$0")/../.." && pwd)"
here="$root/test/adjoint"
flags="-std=f2018 -fcoarray=single -cpp -fPIC -Wno-line-truncation"

( cd "$root" && ./build.sh >/dev/null )

# mesh the adjoint tests need (msh 4.1, on the fly - nothing committed)
bash "$root/meshgen/ensure.sh" "$here/box-36.msh"

cd "$here"

gfortran $flags -I "$root/lib/" -I "$root/src/" -c test.f90 -o test.o
gfortran test.o -o run "$root/lib/libufvm.a" -llapack -fcoarray=single

gfortran $flags -I "$root/lib/" -I "$root/src/" -c test_transient.f90 -o test_transient.o
gfortran test_transient.o -o run_transient "$root/lib/libufvm.a" -llapack -fcoarray=single

echo "==================== steady adjoint ===================="
./run

echo "=================== transient adjoint =================="
./run_transient
