#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [negsing]="singular tolerance is positive"
  [noqstar]="q star has values"
  [noeta]="design direction has values"
  [wide]="tangent rhs is a vector"
  [rhssize]="rhs size matches unknown size"
  [singular]="dense_direct: pivot is nonsingular"
)
for case in negsing noqstar noeta wide rhssize singular; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

migrated="$here/../../src/gti_time_local_tangent_driver.f90"

# the tower migration, statically: the driver solves through the
# existing graph minimization tower - the dense Jacobian rides
# class_graph_stencil inside the adapter, dense_direct eliminates
# it - and no private helper, GTI solver module, or backend
# language survives.
grep -q "class_graph_dense_direct" "$migrated" \
    && echo " PASS : the driver solves through class_graph_dense_direct" \
    || { echo " FAIL : the dense direct import is missing"; exit 1; }
grep -q "class_graph_stencil" "$migrated" \
    && echo " PASS : the stencil representation is named at the seat" \
    || { echo " FAIL : the stencil representation is unnamed"; exit 1; }
grep -q "subroutine dense_solve" "$migrated" \
    && { echo " FAIL : a private dense_solve helper survives"; exit 1; } \
    || echo " PASS : no private dense_solve helper remains"
grep -qE "gti_linear_solve_backend|gti_dense_linear_solve_backend|gti_dense_block_solver" "$migrated" \
    && { echo " FAIL : a GTI solver module is named"; exit 1; } \
    || echo " PASS : no GTI solver module exists"
grep -qi "backend" "$migrated" \
    && { echo " FAIL : backend language appears in the driver"; exit 1; } \
    || echo " PASS : no backend language exists"
