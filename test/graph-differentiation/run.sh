#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [bdforder]="step: only orders one and two carry tables so far"
  [bdfcount]="step: the step count matches the order"
  [bdfstep]="step: a time step is positive"
  [dupslot]="chain_rule: duplicate slot path is refused"
  [badslot]="chain_rule: a path names an input slot"
  [negdegree]="chain_rule: degree is supported"
  [pastcalculus]="chain_rule: the statement supports the requested order"
  [hugedegree]="chain_rule: partition coefficient is representable"
  [statepath]="march_directional: the state path is computed, not supplied"
  [orderzero]="march_directional: the order is positive"
  [unfrozen]="linearization: the tangent is taken at a frozen state"
  [eagerclock]="march_adaptive: the duration is positive"
  [flatcalculus]="operation: the requested order is within max_degree"
  [shallowcalculus]="march_directional: the action's max_degree covers the requested order"
)
for case in bdforder bdfcount bdfstep dupslot badslot negdegree \
            pastcalculus hugedegree statepath orderzero \
            unfrozen eagerclock flatcalculus shallowcalculus; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

# Regression guard: the test sources must not reference the
# deleted gti_* modules.
grep -q "gti_" "$here"/*.f90 \
    && { echo " FAIL : a gti_ reference appears in the suite sources"; exit 1; } \
    || echo " PASS : no gti_ reference in the suite sources"
