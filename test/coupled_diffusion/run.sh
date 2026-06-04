#!/bin/bash
# Coupled-diffusion adjoint gradient check (self-contained: needs only a
# fortran compiler + lapack; not tied to the ufvm library).
set -e

here="$(cd "$(dirname "$0")" && pwd)"
flags="-std=f2018 -O2 -fbacktrace"

cd "$here"
gfortran $flags -c test.f90 -o test.o
gfortran test.o -o run -llapack

./run
