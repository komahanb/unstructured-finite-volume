#!/bin/bash
# Build the library and the regression suite, then run it.
set -e
root="$(cd "$(dirname "$0")/../.." && pwd)"
make -C "$root/src"  F90=gfortran        >/dev/null
make -C "$root/src"  install F90=gfortran >/dev/null
make -C "$root/test/regression" F90=gfortran >/dev/null
cd "$root/test/regression" && ./run
