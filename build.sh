#!/bin/bash
# build the ufvm library into lib/.
# the compiler is auto-detected in Makefile.in, so a bare make just works -
# here, in every example and in the tests. override with `make F90=...`.
set -e

# PRECISION=quad builds the whole tower at real128 into lib_quad/;
# the default is real64 into lib/.
PRECISION=${PRECISION:-double}

mkdir -p lib
make -C src clean PRECISION=$PRECISION
make -C src PRECISION=$PRECISION
make -C src install PRECISION=$PRECISION

echo "library built in lib/ - now 'make' and run an example (e.g. examples/solver)"
