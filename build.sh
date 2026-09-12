#!/bin/bash
# build the ufvm library into lib/.
# the compiler is auto-detected in Makefile.in, so a bare make is sufficient -
# here, in every example and in the tests. override with `make F90=...`.
set -e

# PRECISION=quad builds the whole library at real128 into lib_quad/;
# the default is real64 into lib/.
PRECISION=${PRECISION:-double}

# OPENMP=yes compiles the library with -fopenmp; the default serial
# build treats every !$omp directive as a comment.
OPENMP=${OPENMP:-no}

mkdir -p lib
make -C src clean PRECISION=$PRECISION OPTIMIZE=yes OPENMP=$OPENMP
make -C src PRECISION=$PRECISION OPTIMIZE=yes OPENMP=$OPENMP
make -C src install PRECISION=$PRECISION OPTIMIZE=yes OPENMP=$OPENMP

echo "library built in lib/ - now 'make' and run an example (e.g. examples/solver)"
