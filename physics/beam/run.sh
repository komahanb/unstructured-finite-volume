#!/bin/bash
# Build and run the programs of the Euler-Bernoulli beam:
#
#   euler_bernoulli_beam the beam rhoA w_tt + d_xx(EI w_xx) = f on
#                        [0, 2] x [0, 1], 41 instants, 17 points, no
#                        mesh file: simply supported against its mode
#                        (2.80e-4), J against T EI (pi/L)^4 L/4; clamped
#                        with EI(x), rhoA(x) and a manufactured load
#                        (6.87e-4), and the same on the open line, the
#                        trapezoid rule on the segments; case 1 on the
#                        open line, and in both the error of J
#                        estimated by mu_x: 1.009 of it on the open
#                        line, 0.78 under the fitted differences; case
#                        2 with the Kelvin-Voigt damping c(x) w_xxt, the
#                        mixed components along x of w_t (3.73e-4)
#   euler_bernoulli_beam_optimization
#                        EI(x) = EI_0 (1 + s1 x/L + s2 (x/L)^2) found
#                        where dJ/ds = 0, no condition fixing the
#                        design: minimize iterates on it, the gradient
#                        and hessian from the multiplier of the
#                        design; s = (0.5118, -0.3104) against (0.5,
#                        -0.3) on 21 x 17; with the bound s1 <= 0.45,
#                        the third argument, s = (0.45, -0.2240) and
#                        the bound's multiplier 2.4205e-5
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

launch euler_bernoulli_beam
launch euler_bernoulli_beam_optimization
launch euler_bernoulli_beam_optimization 21 17 0.45
