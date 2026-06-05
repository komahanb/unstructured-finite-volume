#!/bin/bash
# build the ufvm library for DISTRIBUTED (coarray) runs, into lib_par/.
#
# every source is compiled with -fcoarray=lib so the coarray code in
# class_distributed_cg gets a real runtime; the rest of the library has no
# coarray syntax and is unaffected by the flag. the serial build (build.sh,
# -fcoarray=single, into lib/) is left completely alone.
#
# NOTE: uses the Ubuntu-packaged OpenCoarrays wrapper /usr/bin/caf.openmpi.
# the bare `caf` on PATH (/usr/local, 2.9.2) is mismatched against the system
# OpenMPI here and silently breaks co_sum - do not use it. override with CAF=.
set -e

root="$(cd "$(dirname "$0")" && pwd)"
FC="${CAF:-/usr/bin/caf.openmpi}"
flags="-fcoarray=lib -cpp -fPIC -std=f2018 -Wno-line-truncation -O2 -g -fbacktrace -fbounds-check"

mkdir -p "$root/lib_par"
cd "$root/src"

# OBJECTS lists the .o files in dependency order; compile each .f90 with caf,
# placing objects + .mod files in lib_par/.
objs=$(grep -v '^[[:space:]]*#' OBJECTS | sed 's/objects[[:space:]]*=//; s/\\//g')
for o in $objs; do
   f="${o%.o}.f90"
   echo "caf  $f"
   $FC $flags -J "$root/lib_par" -I "$root/lib_par" -c "$f" -o "$root/lib_par/$o"
done

cd "$root/lib_par"
ar rcs libufvm.a $objs
echo "parallel (coarray) library built in lib_par/libufvm.a  (FC=$FC)"
