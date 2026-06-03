#!/bin/bash
# build the ufvm library into lib/.
# the compiler is auto-detected in Makefile.in, so a bare make just works -
# here, in every example and in the tests. override with `make F90=...`.
set -e

mkdir -p lib
make -C src clean
make -C src
make -C src install

echo "library built in lib/ - now 'make' and run an example (e.g. examples/solver)"
