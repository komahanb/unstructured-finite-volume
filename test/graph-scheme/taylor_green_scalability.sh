#!/bin/bash
# The Taylor-Green study on a periodic box of n x n cells over 1, 2,
# 4, ... images: the coarray library and the study built, the box
# meshed by gmsh, and for each image count the wall seconds, the
# forward sweep's seconds and the error against the exact solution
# printed on one line, the error being the same in every digit at
# every image count when the exchange is exact.
#
#     taylor_green_scalability.sh [n] [max_images] [study arguments...]
#
# n defaults to 32 and max_images to the processor count; the study
# arguments follow those of taylor_green_study.f90 after the mesh
# (T, instants, order, dirk, bdf, adams, from2, from3 default to its
# values). MPIRUN_FLAGS adds flags to mpirun (on a shared node
# --bind-to none --mca mpi_yield_when_idle 1 stops the busy-wait of
# idle images from contending with other users).
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"
n=${1:-32}
max_images=${2:-$(nproc)}
shift $(( $# < 2 ? $# : 2 ))
mesh="box_${n}.msh"

cd "$suite_dir"
if [ ! -e "$mesh" ]; then
    sed "s/^n = .*/n = $n;/" box.geo > "box_${n}.geo"
    gmsh -2 "box_${n}.geo" -o "$mesh" > "gmsh_${n}.log"
fi

make -C ../../src clean COARRAY=lib > /dev/null
make -C ../../src COARRAY=lib OPTIMIZE=yes > /dev/null
make -C ../../src install COARRAY=lib OPTIMIZE=yes > /dev/null
make clean > /dev/null
make COARRAY=lib taylor_green_study > /dev/null

echo "box $n x $n, arguments: 1.0 10 2 2 2 2 5 8 $mesh ${*}"
images=1
while [ "$images" -le "$max_images" ]; do
    log="taylor_green_scalability_${n}_np${images}.log"
    t0=$(date +%s.%N)
    if [ "$images" = 1 ]; then
        OMPI_MCA_osc=pt2pt mpirun ${MPIRUN_FLAGS} -np 1 ./taylor_green_study 1.0 10 2 2 2 2 5 8 "$mesh" "$@" > "$log" 2>&1
    else
        mpirun ${MPIRUN_FLAGS} -np "$images" ./taylor_green_study 1.0 10 2 2 2 2 5 8 "$mesh" "$@" > "$log" 2>&1
    fi
    t1=$(date +%s.%N)
    forward=$(grep -a -m1 'took' "$log" | sed 's/.*took *//; s/ s;.*//')
    error=$(grep -a -m1 'error against' "$log" | awk '{print $6}')
    printf "images %3d  wall %8.2f s  forward %8s s  error %s\n" "$images" "$(echo "$t1 - $t0" | bc)" "$forward" "$error"
    images=$(( 2 * images ))
done
