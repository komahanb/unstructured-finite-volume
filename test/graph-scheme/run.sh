#!/bin/bash
# Build the library and run the three programs of this suite:
#
#   dirk_coefficients    the two-stage third-order SDIRK tableau as the
#                        zero of a continuous_residual on a point
#                        manifold; W11 = h (3 + sqrt 3)/6 = 0.078867...
#   taylor_green_vortex  the Navier-Stokes Taylor-Green vortex on the
#                        periodic box of box.geo, 16 x 16 cells, ten
#                        instants, the chain dirk(2), bdf(2), adams(2);
#                        the error against the exact solution, second
#                        order in the cell width (8.93e-2 at 8 x 8,
#                        2.33e-2 at 16 x 16)
#   taylor_green_vortex_3d
#                        the vortex extended along z on the periodic
#                        8 x 8 x 8 box of box_3d.py, five instants,
#                        finite volumes of order 2
#   van_der_pol          the van der Pol oscillator on [0, 2] by the
#                        chain dirk(2), bdf(2), adams(2) over eleven
#                        instants; u(2) = 0.2912 (0.3233 converged)
#   van_der_pol_sensitivity
#                        the oscillator with nu and mu as design
#                        coordinates, the energy as the objective, the
#                        adjoint, the reactions of the initial data,
#                        the multipliers of the design conditions
#                        (the gradient of the energy) with a row of
#                        its Hessian, and the derivatives of u(T) and
#                        of the energy along each coordinate to order
#                        3 by the expansion along the design
#   van_der_pol_continuous_adjoint
#                        the continuous adjoint of the oscillator,
#                        formed as the stationarity of the Lagrangian
#                        density and solved with the state as one
#                        two-point block, against the discrete adjoint
#   van_der_pol_study    the Van der Pol case of the graph time
#                        integrators paper: its Table 4 to degree four
#                        at h = 0.02, and its three-window chain
#   euler_bernoulli_beam the beam rhoA w_tt + d_xx(EI w_xx) = f on
#                        [0, 2] x [0, 1], 41 instants, 17 points, no
#                        mesh file: simply supported against its mode
#                        (2.80e-4), J against T EI (pi/L)^4 L/4; clamped
#                        with EI(x), rhoA(x) and a manufactured load
#                        (6.87e-4), and the same on the open line, the
#                        trapezoid rule on the segments
#   euler_bernoulli_beam_optimization
#                        EI(x) = EI_0 (1 + s1 x/L + s2 (x/L)^2) found
#                        by Newton on the design, gradient and hessian
#                        from the multiplier of the design condition:
#                        s = (0.5096, -0.3107) against (0.5, -0.3) on
#                        21 x 17, the gradient 8.3e-12
#   run                  the same tableau on the expression and the
#                        residual_operator directly
#
# gmsh meshes box.geo into box.msh when the mesh is absent.
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null )
fi

make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null
make -C "$suite_dir" box.msh box_3d.msh >/dev/null

cd "$suite_dir"
echo "== dirk_coefficients"
./dirk_coefficients
echo "== taylor_green_vortex"
./taylor_green_vortex
echo "== taylor_green_vortex_3d"
./taylor_green_vortex_3d
echo "== van_der_pol"
./van_der_pol
echo "== van_der_pol_sensitivity"
./van_der_pol_sensitivity
echo "== van_der_pol_continuous_adjoint"
./van_der_pol_continuous_adjoint
echo "== van_der_pol_study"
./van_der_pol_study
echo "== euler_bernoulli_beam"
./euler_bernoulli_beam
echo "== euler_bernoulli_beam_optimization"
./euler_bernoulli_beam_optimization
echo "== run"
./run
