#!/bin/bash
# Build and run the programs of the van der Pol oscillator:
#
#   van_der_pol          the oscillator on [0, 2] by the chain dirk(2),
#                        bdf(2), adams(2) over eleven instants; u(2) =
#                        0.2912 (0.3233 converged)
#   van_der_pol_sensitivity
#                        the oscillator with nu and mu as design
#                        coordinates, the energy as the objective, the
#                        adjoint, the reactions of the initial data,
#                        the multipliers of the design conditions
#                        (the gradient of the energy) with a row of
#                        its Hessian, and the derivatives of u(T) and
#                        of the energy along each coordinate to order
#                        3 by the expansion along the design
#   van_der_pol_philosophical_template
#                        the oscillator stated once; order 0 of the
#                        expansion along nu is van_der_pol in every
#                        digit
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
launch van_der_pol_sensitivity
launch van_der_pol_philosophical_template
