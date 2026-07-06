#!/bin/bash
# build, generate the mesh on the fly (msh 4.1, gmsh), and run the
# assembly test (A = L + U + D on box-36).
set -e

here="$(cd "$(dirname "$0")" && pwd)"

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

bash "$here/../../meshgen/ensure.sh" "$here/../box-36.msh"

cd "$here" && ./run

# expected failure: the REVERSE product on an undeclared, un-overridden
# system must refuse (nonzero exit)
if ./refusal >/dev/null 2>&1; then
  echo " FAIL : transpose refusal did not fire"
  exit 1
else
  echo " PASS : transpose refusal fires on an undeclared system"
fi
