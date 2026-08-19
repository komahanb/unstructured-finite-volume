#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [negiter]="max_iterations is positive"
  [negrtol]="residual tolerance is nonnegative"
  [negstep]="step tolerance is nonnegative"
  [negsing]="singular tolerance is positive"
  [novalues]="newton_driver: trial q has values"
  [dqshape]="direction q shape matches trial q"
  [wide]="residual is a vector"
  [short]="residual size matches unknown size"
  [shortjac]="newton_driver: residual size matches unknown size"
  [singular]="dense_direct: pivot is nonsingular"
)
for case in negiter negrtol negstep negsing novalues dqshape wide short shortjac singular; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

migrated="$here/../../src/gti_time_local_newton_driver.f90"

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
grep -q "subroutine driver_jacobian_action" "$migrated" \
    && grep -q "subroutine driver_residual" "$migrated" \
    && echo " PASS : residual and Jacobian action still live in the driver" \
    || { echo " FAIL : Newton lost more than its linear solve"; exit 1; }

# zero duplication, statically: the Jacobian probe loop exists ONCE
# in all of src - inside this driver's dense_jacobian verb - and
# every seat that eliminates J whole calls it.
srcdir="$here/../../src"
probes=$(grep -l "jacobian(:, j) = column_values" "$srcdir"/*.f90 | wc -l)
[ "$probes" -eq 1 ] && grep -q "jacobian(:, j) = column_values" "$srcdir/gti_time_local_newton_driver.f90" \
    && echo " PASS : the Jacobian probe loop exists once, in dense_jacobian" \
    || { echo " FAIL : the probe loop is duplicated ($probes files)"; exit 1; }
grep -q "subroutine driver_dense_jacobian" "$srcdir/gti_time_local_newton_driver.f90" \
    && echo " PASS : dense_jacobian is the one assembly verb" \
    || { echo " FAIL : the assembly verb is missing"; exit 1; }
