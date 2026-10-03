#!/bin/bash
# Build and run the programs of the dynamics of a rotorcraft:
#
#   triple_pendulum      the triple pendulum of the dissertation,
#                        section 5.1, over two seconds by bdf(2): the
#                        energy at every tenth of the interval and the
#                        largest condition of the joints, below 1e-15;
#                        the instant count from the command line
#                        (default 201)
#   rotorcraft_dynamics  the rotor of komahanb/rotorcraft-optimization,
#                        26 bodies and 34 joints, one revolution by
#                        bdf(1), the setup read from rotorcraft.nml:
#                        the blade model, the scenario, the steps, the
#                        meshes the bodies are drawn with
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
echo "== rotorcraft_dynamics"
./rotorcraft_dynamics
