#!/bin/bash
# build the library (serial + coarray) and the distributed-CG suite, run it
# serially (1 image, the reduction-to-serial check) and on 2 and 4 images.
#
# uses the Ubuntu-packaged OpenCoarrays (/usr/bin/caf.openmpi); the bare caf
# on PATH is broken here. override with CAF= / CAFRUN=.
set -e

root="$(cd "$(dirname "$0")/../.." && pwd)"
here="$root/test/parallel"
sflags="-std=f2018 -fcoarray=single -cpp -fPIC -Wno-line-truncation"
pflags="-std=f2018 -cpp -fPIC -Wno-line-truncation -O2 -fbounds-check"
CAF="${CAF:-/usr/bin/caf.openmpi}"
CAFRUN="${CAFRUN:-/usr/bin/cafrun.openmpi}"

# serial library (lib/) and parallel coarray library (lib_par/)
( cd "$root" && ./build.sh >/dev/null )
( cd "$root" && bash build_parallel.sh >/dev/null )

# meshes the suite needs
for m in box-36 square-20 square-40; do
   bash "$root/meshgen/ensure.sh" "$here/$m.msh"
done

cd "$here"
unset DISPLAY    # drop harmless X11 "MIT-MAGIC-COOKIE" noise from MPI startup

# --- serial: 1 image, distributed_cg must reduce to the serial CG ---
echo "### serial build (1 image, -fcoarray=single) ###"
gfortran $sflags -I "$root/lib/" -c test.f90 -o test_s.o
gfortran test_s.o -o run_serial "$root/lib/libufvm.a" -fcoarray=single
./run_serial

# --- distributed: real coarray runs on 2 and 4 images ---
$CAF $pflags -I "$root/lib_par/" -c test.f90 -o test_p.o
$CAF test_p.o -o run "$root/lib_par/libufvm.a"
for np in 2 4; do
   echo "### cafrun -np $np ###"
   $CAFRUN -np $np ./run
done
