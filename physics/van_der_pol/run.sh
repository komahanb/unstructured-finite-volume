#!/bin/bash
# Build and run the oscillator to follow the energy question from
# motion to sensitivity, design and the accuracy of the computation.
# We state the motion on [0, 2] x {nu}, pairing each law r_l with its
# multiplier lambda_l. Starting without J gives the physics (a).
# Adding J gives the adjoint (b) and sensitivity (c); this Lagrangian
# is printed. Removing the pin of nu lets us choose the design (d).
# We then place instants by its estimated error (e), and choose their
# count to meet a design tolerance (f). Returning to J, we choose
# instants and families from a pool to meet an energy tolerance (g)
# or to use a prescribed number of solves (h).
#
# The first run shows the motion and design with eleven instants and
# order 1 along nu. Its recorded u(2) = 0.2912 (0.3233 converged),
# and nu = 1.5269. It estimates the error of J both at nu_0 and at
# the design found:
#
#     van_der_pol
#
# The second run uses 21 instants and order 3 to obtain derivatives
# of u(T) and J along nu through order 3. Comparing the derivatives
# of J from the state and adjoint expansions audits the sensitivity:
#
#     van_der_pol 21 2 0 3
#
# Use run.sh for one process. With IMAGES=n, the library and programs
# are built with COARRAY=lib and run by mpirun -np n on n coarray
# images. MPIRUN_FLAGS supplies any additional launcher flags.
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
    # Let the Makefile rebuild the coarray library when a source in
    # src is newer, so the program uses the current library.
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
