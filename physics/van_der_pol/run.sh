#!/bin/bash
# Build and run the program of the van der Pol oscillator, stated once
# on [0, 2] x {nu}, every law r_l paired with its multiplier lambda_l,
# and four problems as Lagrangians of that statement: (a) the physics,
# without J; (b) the adjoint and (c) the sensitivity, J added, that
# Lagrangian printed; (d) the design, the pin of nu removed. Two runs:
#
#   van_der_pol          eleven instants, order 1 along nu: u(2) =
#                        0.2912 (0.3233 converged), nu = 1.5269 found,
#                        and the estimate of the error of J at nu_0
#                        and at the design found
#   van_der_pol 21 2 0 3 21 instants, order 3: the derivatives of u(T)
#                        and of J along nu to order 3, each of J by the
#                        state and by the adjoint
#
#     run.sh             on one process
#     IMAGES=n run.sh    on n coarray images: the library and the
#                        programs built with COARRAY=lib, each program
#                        run by mpirun -np n (MPIRUN_FLAGS adds flags)
set -e

dir="$(cd "$(dirname "$0")" && pwd)"
images=${IMAGES:-0}

if [ "$images" -eq 0 ]; then
    if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
        ( cd "$dir/../.." && ./build.sh >/dev/null )
    fi
    make -C "$dir" clean >/dev/null 2>&1 || true
    make -C "$dir" >/dev/null
else
    # the coarray library is rebuilt by the Makefile when a source of
    # src is newer than it
    make -C "$dir" clean >/dev/null 2>&1 || true
    make -C "$dir" COARRAY=lib >/dev/null
fi

cd "$dir"
launch() {
    echo "== $*"
    if [ "$images" -eq 0 ]; then
        "./$1" "${@:2}"
    elif [ "$images" -eq 1 ]; then
        OMPI_MCA_osc=pt2pt mpirun ${MPIRUN_FLAGS:-} -np 1 "./$1" "${@:2}"
    else
        mpirun ${MPIRUN_FLAGS:-} -np "$images" "./$1" "${@:2}"
    fi
}

launch van_der_pol
launch van_der_pol 21 2 0 3
