#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
fi
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [zero_tolerance]="dense_direct: singular_tolerance must be positive"
  [size_mismatch]="dense_direct: size(x) must equal size(rhs)"
  [singular]="dense_direct: factorisation found a singular pivot"
  [nonsquare]="stencil: a dense matrix must be square"
  [nonintegral_width]="stencil: the width must be a positive whole multiple of n_dom"
  [incompatible_components]="operation: the bound field does not satisfy its argument contract"
)
for case in zero_tolerance size_mismatch singular nonsquare nonintegral_width incompatible_components; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused with the expected diagnostic" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

src="$here/../../src/operation_dense_direct.f90"

# Static checks on src/operation_dense_direct.f90: dense_direct
# must extend minimizer, must import operation_minimization, and
# must not reference gti_ modules, the word "backend", or define a
# solve_transpose method. The dense-array adapters stay deleted:
# the solver does not name them outside a comment.
grep -q "extends(minimizer) :: dense_direct" "$src" \
    && echo " PASS : dense_direct extends operation_minimization's minimizer" \
    || { echo " FAIL : dense_direct is not a minimizer concretion"; exit 1; }
grep -q "use operation_minimization" "$src" \
    && echo " PASS : the solver imports operation_minimization" \
    || { echo " FAIL : a required import is missing"; exit 1; }
grep -hv '^ *!' "$src" \
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
