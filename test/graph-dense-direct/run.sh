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
  [nonsquare]="dense_direct: dense matrix is square"
)
for case in tolzero sizemismatch singular nonsquare; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

src="$here/../../src/class_graph_dense_direct.f90"

# the seat's place in the tower, statically: dense_direct extends
# the minimizer family and rides the stencil representation - no
# GTI vocabulary, no second hierarchy, no backend language.
grep -q "extends(minimizer) :: dense_direct" "$src" \
    && echo " PASS : dense_direct extends graph_minimization's minimizer" \
    || { echo " FAIL : dense_direct is not a minimizer concretion"; exit 1; }
grep -q "use graph_minimization" "$src" && grep -q "use class_graph_stencil" "$src" \
    && echo " PASS : the seat imports the tower it joins" \
    || { echo " FAIL : the tower imports are missing"; exit 1; }
grep -v '^ *!' "$src" | grep -qE "gti_" \
    && { echo " FAIL : a GTI seat reaches the solver tower"; exit 1; } \
    || echo " PASS : no GTI import or name exists in the solver seat"
grep -qi "backend" "$src" \
    && { echo " FAIL : backend language appears in the solver seat"; exit 1; } \
    || echo " PASS : no backend language exists"
grep -q "solve_transpose" "$src" \
    && { echo " FAIL : a solve_transpose method exists"; exit 1; } \
    || echo " PASS : no solve_transpose method - a transpose is a stencil"

# word sweeps: 'Bell' checked case-sensitively as the proper noun.
grep -q "fiber" "$src" \
    && { echo " FAIL : the word 'fiber' appears in the seat"; exit 1; } \
    || echo " PASS : no 'fiber' in the seat"
grep -q "Bell" "$src" \
    && { echo " FAIL : the word 'Bell' appears in the seat"; exit 1; } \
    || echo " PASS : no 'Bell' in the seat"
grep -q "jet" "$src" \
    && { echo " FAIL : the word 'jet' appears in the seat"; exit 1; } \
    || echo " PASS : no 'jet' in the seat"
