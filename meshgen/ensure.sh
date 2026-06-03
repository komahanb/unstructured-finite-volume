#!/bin/bash
# generate the mesh at <path> (msh 4.1, via gmsh) from its basename's recipe.
# the basename must be a known recipe in generate.py, e.g.
#   ensure.sh ../../test/box-3.msh   ->   generate.py box-3 ../../test/box-3.msh
set -e

[ -n "$1" ] || { echo "usage: ensure.sh <meshpath>" >&2; exit 1; }

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "$here/generate.py" "$(basename "$1" .msh)" "$1"
