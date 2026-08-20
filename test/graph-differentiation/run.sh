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
  [dupslot]="chain_rule: duplicate argument path is refused"
  [badslot]="operation: the argument is declared"
  [foreignpath]="chain_rule: a path names an argument of the statement"
  [foreignvariation]="operation: a variation names an argument of the operation"
  [undeclared]="operation: the argument space is declared before an argument is named"
  [historyreach]="step: the history argument is within the reach"
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
for case in bdforder bdfcount bdfstep dupslot badslot foreignpath \
            foreignvariation undeclared historyreach negdegree \
            pastcalculus hugedegree statepath orderzero \
            unfrozen eagerclock flatcalculus shallowcalculus; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

# Constructor-bypass audit: every concrete operation, in src/ and in
# the suites, declares its argument space in the file that defines
# it - a module constructor or an attach entry calls
# declare_arguments. An operation built by default initialization
# owns no arguments and is refused by the 'undeclared' case above.
root="$here/../.."
missing=$(grep -rlE "extends\((operation|discretization)\)" "$root/src" "$root/test" --include=*.f90 \
    | while read -r f; do grep -qE "^ *type *, *extends\((operation|discretization)\) *::" "$f" || continue; \
        grep -q "declare_arguments(" "$f" || echo "$f"; done)
if [ -n "$missing" ]; then
    echo " FAIL : operations without a declared argument space:"; echo "$missing"; exit 1
else
    echo " PASS : every concrete operation declares its argument space"
fi

# Regression guard: the test sources must not reference the
# deleted gti_* modules.
grep -q "gti_" "$here"/*.f90 \
    && { echo " FAIL : a gti_ reference appears in the suite sources"; exit 1; } \
    || echo " PASS : no gti_ reference in the suite sources"
