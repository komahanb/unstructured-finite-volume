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
  [dupslot]="chain_rule: duplicate slot channel is refused"
  [badslot]="chain_rule: a channel names an input slot"
  [negdegree]="chain_rule: degree is supported"
  [pastcalculus]="chain_rule: the statement supports the requested order"
  [hugedegree]="chain_rule: pattern coefficient is representable"
  [statepath]="march_directional: the state path is the walk's own"
  [orderzero]="march_directional: the order is positive"
  [blindreverse]="march_adjoint: the implicit reverse walk needs the action"
  [unfrozen]="linearization: an exact tangent is taken at a frozen state"
  [eagerclock]="march_adaptive: the duration is positive"
)
for case in bdforder bdfcount bdfstep dupslot badslot negdegree \
            pastcalculus hugedegree statepath orderzero blindreverse \
            unfrozen eagerclock; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

# the book stays closed, statically: this suite speaks only the
# living tower - no GTI vocabulary anywhere in its sources.
grep -q "gti_" "$here"/*.f90 \
    && { echo " FAIL : a GTI name appears in the suite"; exit 1; } \
    || echo " PASS : the suite speaks only the living tower"
