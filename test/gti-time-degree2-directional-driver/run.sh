#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [negsing]="degree2_directional_driver: singular tolerance is positive"
  [norelations]="degree2_directional_driver: at least one relation is required"
  [idx0]="degree2_directional_driver: relation index is in range"
  [idxhigh]="degree2_directional_driver: relation index is in range"
  [firstsize]="derivative array size matches graph vertices"
  [secondsize]="derivative array size matches graph vertices"
  [noeta]="degree2_directional_driver: design direction has values"
  [unsolved]="degree2_directional_driver: unknown vertex has solution"
  [unsolvedhistory]="degree2_directional_driver: history vertex has solution"
  [noq]="degree2_directional_driver: unknown q has values"
  [badseed]="derivative shape matches vertex q"
  [shortres]="degree residual size matches unknown size"
  [singular]="dense_direct: pivot is nonsingular"
)
for case in negsing norelations idx0 idxhigh firstsize secondsize noeta unsolved unsolvedhistory noq badseed shortres singular; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
# check 17: the higher-order chain-rule seat is genuinely the
# assembler of B1 and B2 - the import is in the production source.
grep -q "use gti_chain_rule_assemblies" \
    "$here/../../src/gti_time_degree2_directional_driver.f90" \
    && echo " PASS : the production source imports gti_chain_rule_assemblies" \
    || { echo " FAIL : chain-rule assembler import missing"; exit 1; }

migrated="$here/../../src/gti_time_degree2_directional_driver.f90"

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
