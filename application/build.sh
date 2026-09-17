#!/bin/bash
# Build the packed graph-time-integrator application.
#
# src/ is built by ../build.sh into ../lib, which holds the library modules
# and objects this executable links against. The application layer is one
# source file: modules first, main program last.
set -e

cd "$(dirname "$0")"
F90=${F90:-gfortran-15}
FLAGS="-std=f2023 -fcoarray=single -cpp -Wall -fbounds-check -O2"

# OPENMP=yes compiles and links with -fopenmp, as the library must
# have been built (../build.sh with the same setting).
OPENMP=${OPENMP:-no}
if [ "$OPENMP" = yes ]; then
   FLAGS="$FLAGS -fopenmp"
fi

# PRECISION=quad links against lib_quad and puts its binaries in
# quad/, so a double and a quadruple build coexist.
PRECISION=${PRECISION:-double}
LIB=../lib
OBJ=.build
OUT=.
if [ "$PRECISION" = quad ]; then
   FLAGS="$FLAGS -DPRECISION_QUAD"
   LIB=../lib_quad
   OBJ=.build_quad
   OUT=quad
fi
mkdir -p $OBJ $OUT

SOURCE=module_graph_time_integrator.f90
PROGRAM=graph_time_integrator

$F90 $FLAGS -I$LIB -J$OBJ -o $OUT/$PROGRAM $SOURCE $LIB/*.o

echo "built: $PROGRAM"

# The two programs on the continuous/discrete residual interface: the
# SDIRK tableau of order three and the Taylor-Green vortex, whose
# periodic box is meshed from box.geo by gmsh when box.msh is absent.
[ -f box.msh ] || gmsh -2 box.geo -o box.msh > gmsh.log
for PROGRAM in dirk_coefficients taylor_green_vortex; do
   $F90 $FLAGS -I$LIB -J$OBJ -o $OUT/$PROGRAM $PROGRAM.f90 $LIB/*.o
   echo "built: $PROGRAM"
done
