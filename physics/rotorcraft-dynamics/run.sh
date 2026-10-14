#!/bin/bash
# Build and run the programs of the dynamics of a rotorcraft:
#
#   triple_pendulum      the triple pendulum of the dissertation,
#                        section 5.1, over two seconds by bdf(2): the
#                        energy at every tenth of the interval and the
#                        largest condition of the joints, below 1e-15;
#                        the instant count from the command line
#                        (default 201)
#   triple_pendulum_sensitivity
#                        the pendulum with a design in a mass and in a
#                        joint: dJ/ds by the adjoint and by the
#                        expansion, agreeing to 2.9e-13
#   blade_on_hub         the blade of the rotor on a turning hub, as a
#                        beam in the hub's axes or as a rigid plate:
#                        the driver's reactions and the tip's motion;
#                        the beam approaches the plate as E grows
#   rotorcraft_dynamics  the rotor of komahanb/rotorcraft-optimization,
#                        26 bodies and 35 joints, one revolution by
#                        bdf(2), the setup read from rotorcraft.nml:
#                        the blade model, the scenario, the steps, the
#                        meshes the bodies are drawn with
#
# Each program prints on image 1 the convergence of every solve, and of
# the design's iteration, before its results, as van_der_pol does.
#
#     run.sh             on one process
set -e

dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$dir/../.." && ./build.sh >/dev/null )
fi
make -C "$dir" clean >/dev/null 2>&1 || true
make -C "$dir" >/dev/null

cd "$dir"
echo "== triple_pendulum"
./triple_pendulum
echo "== triple_pendulum_sensitivity"
./triple_pendulum_sensitivity
echo "== blade_on_hub"
./blade_on_hub
echo "== rotorcraft_dynamics"
./rotorcraft_dynamics
