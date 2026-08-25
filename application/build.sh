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
FLAGS="-std=f2023 -fcoarray=single -cpp -Wall -fbounds-check"
OBJ=.build
mkdir -p $OBJ

MODULES="../physics/physics_integrand ../physics/physics_vanderpol
         gti_configuration gti_sweeps gti_expansion gti_block gti_stage
         gti_march gti_taylor gti_chain"

for m in $MODULES; do
   $F90 $FLAGS -I../lib -J$OBJ -c $m.f90 -o $OBJ/$(basename $m).o
done

DRIVERS=${*:-$(ls *.f90 | grep -v '^gti_' | sed 's/\.f90$//')}

for d in $DRIVERS; do
   $F90 $FLAGS -I../lib -I$OBJ -J$OBJ -o $d $d.f90 $OBJ/*.o ../lib/*.o
done

echo "built: $DRIVERS"
