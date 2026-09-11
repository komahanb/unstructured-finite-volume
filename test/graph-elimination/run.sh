#!/bin/bash
# Sparse elimination must store coefficients, not dependency paths.
set -e
here="$(cd "$(dirname "$0")" && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$here/../.." && ./build.sh >/dev/null )
fi
make -C "$here" clean >/dev/null
make -C "$here" >/dev/null
# The forty-row case has only two retained columns. Enumerating all
# substitution paths would require gigabytes; the coefficients do not.
( ulimit -v 262144; cd "$here"; ./run )
