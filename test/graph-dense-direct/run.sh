#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [tolzero]="dense_direct: singular tolerance is positive"
  [sizemismatch]="dense_direct: solution size matches rhs"
  [singular]="dense_direct: pivot is nonsingular"
  [nonsquare]="stencil: a dense matrix is square"
  [badwidth]="stencil: the width carries a whole number per member"
  [badresult]="stencil: the operation result matches the width"
)
for case in tolzero sizemismatch singular nonsquare badwidth badresult; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

src="$here/../../src/operation_dense_direct.f90"

# Static checks on src/operation_dense_direct.f90: dense_direct
# must extend minimizer, must import operation_minimization, and
# must not reference gti_ modules, the word "backend", or define a
# solve_transpose method. The dense-array adapters stay deleted:
# neither the solver nor the marcher names them outside a comment.
grep -q "extends(minimizer) :: dense_direct" "$src" \
    && echo " PASS : dense_direct extends operation_minimization's minimizer" \
    || { echo " FAIL : dense_direct is not a minimizer concretion"; exit 1; }
grep -q "use operation_minimization" "$src" \
    && echo " PASS : the solver imports operation_minimization" \
    || { echo " FAIL : a required import is missing"; exit 1; }
grep -hv '^ *!' "$src" "$here/../../src/operation_marching.f90" \
    | grep -v '% transpose()' \
    | grep -qE "dense_matrix_of|solve_dense_matrix|jstep|transpose\(" \
    && { echo " FAIL : a dense-array adapter or a hand-built step matrix is referenced"; exit 1; } \
    || echo " PASS : the dense-array adapters stay deleted"
grep -v '^ *!' "$src" | grep -qE "gti_" \
    && { echo " FAIL : a gti_ reference appears in the solver source"; exit 1; } \
    || echo " PASS : no gti_ reference exists in the solver source"
grep -qi "backend" "$src" \
    && { echo " FAIL : the word 'backend' appears in the solver source"; exit 1; } \
    || echo " PASS : no backend language exists"
grep -q "solve_transpose" "$src" \
    && { echo " FAIL : a solve_transpose method exists"; exit 1; } \
    || echo " PASS : no solve_transpose method exists"

# Banned-word checks; 'Bell' is matched case-sensitively so that
# words merely containing 'bell' do not trigger it.
grep -q "fiber" "$src" \
    && { echo " FAIL : the word 'fiber' appears in the solver source"; exit 1; } \
    || echo " PASS : no 'fiber' in the solver source"
grep -q "Bell" "$src" \
    && { echo " FAIL : the word 'Bell' appears in the solver source"; exit 1; } \
    || echo " PASS : no 'Bell' in the solver source"
grep -q "jet" "$src" \
    && { echo " FAIL : the word 'jet' appears in the solver source"; exit 1; } \
    || echo " PASS : no 'jet' in the solver source"
