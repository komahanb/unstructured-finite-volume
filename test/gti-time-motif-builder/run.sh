#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [bdf0]="BDF order is supported"
  [bdf3]="BDF order is supported"
  [bdfh0]="step size is positive"
  [bdfhneg]="step size is positive"
  [dirkh0]="step size is positive"
  [dirkg0]="DIRK diagonal is nonzero"
  [abm0]="ABM order is supported"
  [abm3]="ABM order is supported"
  [abmh0]="step size is positive"
  [vbdf0]="BDF order is supported"
  [vbdf3]="BDF order is supported"
  [vbdfcount1]="BDF step count matches order"
  [vbdfcount2]="BDF step count matches order"
  [vbdfh0]="time step is positive"
  [vbdfhneg]="time step is positive"
  [vabm3]="ABM order is supported"
  [vabmh0]="time step is positive"
  [ftcount]="BDF takes two or three times"
  [ftorder]="times are increasing"
)
for case in bdf0 bdf3 bdfh0 bdfhneg dirkh0 dirkg0 abm0 abm3 abmh0 \
            vbdf0 vbdf3 vbdfcount1 vbdfcount2 vbdfh0 vbdfhneg \
            vabm3 vabmh0 ftcount ftorder; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

src="$here/../../src/gti_time_motif_builder.f90"

# the builder's boundary, statically: rows out of names and steps,
# nothing else - no graph, no traversal, no chain rule, no change
# protocol, no attached values, no mesh, no policy.
grep -E '^ *use ' "$src" | grep -vE 'iso_fortran_env|gti_form_interface|gti_time_local_schemes' \
    && { echo " FAIL : the builder imports beyond its three seats"; exit 1; } \
    || echo " PASS : the builder imports only its three lawful seats"

# word sweeps: 'Bell' checked case-sensitively as the proper noun.
grep -q "fiber" "$src" \
    && { echo " FAIL : the word 'fiber' appears in the builder"; exit 1; } \
    || echo " PASS : no 'fiber' in the builder"
grep -q "Bell" "$src" \
    && { echo " FAIL : the word 'Bell' appears in the builder"; exit 1; } \
    || echo " PASS : no 'Bell' in the builder"
grep -q "jet" "$src" \
    && { echo " FAIL : the word 'jet' appears in the builder"; exit 1; } \
    || echo " PASS : no 'jet' in the builder"
