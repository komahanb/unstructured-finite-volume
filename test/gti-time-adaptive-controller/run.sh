#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [attempts0]="gti_time_adaptive_controller: attempt budget is positive"
  [novertex]="gti_time_adaptive_controller: the graph carries an initial vertex"
  [badstep]="gti_time_adaptive_controller: policy proposes a positive step"
  [badorder]="gti_time_adaptive_controller: policy proposes a supported order"
  [earlyorder2]="gti_time_adaptive_controller: order two needs two history vertices"
  [polinit]="gti_time_halving_step_policy: initial step is positive"
  [poltol]="gti_time_halving_step_policy: error tolerance is positive"
  [polorder]="gti_time_halving_step_policy: preferred order is supported"
)
for case in attempts0 novertex badstep badorder earlyorder2 polinit poltol polorder; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

src="$here/../../src/gti_time_adaptive_controller.f90"

# the roadmap's sharp criterion, statically: the controller is a
# client of the growth seat and the motif builders, and of nothing
# that executes - no traversal driver, no local calculus, no chain
# rule, no change protocol, no attached values, no kernel, no mesh.
grep -q "use gti_time_adaptive_growth_drivers" "$src" \
    && echo " PASS : the controller is a client of the growth seat" \
    || { echo " FAIL : the growth import is missing"; exit 1; }
grep -q "use gti_time_motif_builders" "$src" \
    && echo " PASS : the controller mints candidates from the motif builders" \
    || { echo " FAIL : the motif builder import is missing"; exit 1; }
grep -q "use gti_time_forward_drivers, only : gti_time_forward_options" "$src" \
    && echo " PASS : the forward import carries only the options type" \
    || { echo " FAIL : the forward import is wider than the options type"; exit 1; }
grep -E '^ *use ' "$src" | grep -vE 'iso_fortran_env|gti_value_buffers|gti_state_bundles|gti_design_bundles|gti_form_interface|gti_time_graphs|gti_time_motif_builders|gti_time_forward_drivers, only : gti_time_forward_options|gti_time_adaptive_growth_drivers' \
    && { echo " FAIL : the controller imports beyond its lawful seats"; exit 1; } \
    || echo " PASS : no traversal, calculus, chain-rule, protocol, map, or mesh import"

# word sweeps: 'Bell' checked case-sensitively as the proper noun.
grep -q "fiber" "$src" \
    && { echo " FAIL : the word 'fiber' appears in the controller"; exit 1; } \
    || echo " PASS : no 'fiber' in the controller"
grep -q "Bell" "$src" \
    && { echo " FAIL : the word 'Bell' appears in the controller"; exit 1; } \
    || echo " PASS : no 'Bell' in the controller"
grep -q "jet" "$src" \
    && { echo " FAIL : the word 'jet' appears in the controller"; exit 1; } \
    || echo " PASS : no 'jet' in the controller"
