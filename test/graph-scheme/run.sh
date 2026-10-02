#!/bin/bash
# Build the library and run the programs of this suite: the studies,
# the audits and the derivations of schemes; the problems they read
# are stated in physics/.
#
#   dirk_coefficients    the two-stage third-order SDIRK tableau as the
#                        zero of a continuous_residual on a point
#                        manifold; W11 = h (3 + sqrt 3)/6 = 0.078867...
#   van_der_pol_continuous_adjoint
#                        the continuous adjoint of the oscillator,
#                        formed as the stationarity of the Lagrangian
#                        density and solved with the state as one
#                        two-point block, against the discrete adjoint
#   van_der_pol_study    the Van der Pol case of the graph time
#                        integrators paper: its Table 4 to degree four
#                        at h = 0.02, and its three-window chain
#   van_der_pol_design_audit
#                        the expansion along the design audited on two
#                        design coordinates, the initial velocity and
#                        the objective reading them: the derivatives of
#                        J by the expansion of the state against the
#                        adjoint's along nu, mu and nu + mu, the mixed
#                        derivative read from both rows of the Hessian
#                        and by polarization, and the conditions summed
#                        bare; every difference is rounding, 1.4e-14 at
#                        most, and the program stops above 1e-10
#   euler_bernoulli_beam_grid
#                        the points of the beamline graded by one design
#                        coordinate s, moved by Newton on the sum of the
#                        squares of the estimate's contributions: s =
#                        0.372, the error of J 2.34e-2 to -7.4e-4
#   run                  the same tableau on the expression and the
#                        residual_operator directly
#
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null )
fi

make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null

cd "$suite_dir"
echo "== dirk_coefficients"
./dirk_coefficients
echo "== van_der_pol_continuous_adjoint"
./van_der_pol_continuous_adjoint
echo "== van_der_pol_study"
./van_der_pol_study
echo "== van_der_pol_design_audit"
./van_der_pol_design_audit
echo "== euler_bernoulli_beam_grid"
./euler_bernoulli_beam_grid
echo "== run"
./run
