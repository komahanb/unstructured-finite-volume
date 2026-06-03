#!/bin/bash
# build, generate the config's mesh on the fly (msh 4.1, gmsh), and march.
#   bash run.sh heat.cfg     (default: heat.cfg)
set -e

here="$(cd "$(dirname "$0")" && pwd)"
cfg="${1:-heat.cfg}"

# clean rebuild so the driver relinks against the current libufvm
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

mesh="$(awk '$1=="mesh"{print $2}' "$here/$cfg")"
bash "$here/../../meshgen/ensure.sh" "$here/$mesh"

cd "$here" && ./run "$cfg"
