#!/bin/bash
# build the gti modules and every driver in application/.
#
# src/ is built by ../build.sh into ../lib, which holds the .mod and .o
# files this links against. the module compile order below is the
# dependency order of the `use` statements; objects and .mod files go
# to .build/ so the source directory holds only sources and binaries.
set -e

cd "$(dirname "$0")"
F90=${F90:-gfortran-15}
FLAGS="-std=f2023 -fcoarray=single -cpp -Wall -fbounds-check -O2"

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

MODULES="../physics/physics_vanderpol
         gti_configuration gti_sweeps gti_expansion gti_block
         gti_march gti_adaptive gti_space gti_field gti_chain gti_driver"

for m in $MODULES; do
   $F90 $FLAGS -I$LIB -J$OBJ -c $m.f90 -o $OBJ/$(basename $m).o
done

DRIVERS=${*:-$(ls *.f90 | grep -v '^gti_' | sed 's/\.f90$//')}

for d in $DRIVERS; do
   $F90 $FLAGS -I$LIB -I$OBJ -J$OBJ -o $OUT/$d $d.f90 $OBJ/*.o $LIB/*.o
done

echo "built: $DRIVERS"
