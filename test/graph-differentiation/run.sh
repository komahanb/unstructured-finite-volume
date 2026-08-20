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
  [historyshape]="step: a history state matches the state's storage shape"
  [historymissing]="step: the history state is given"
  [unsupplied]="marcher: every auxiliary argument of the action is supplied as a parameter"
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
            foreignvariation undeclared historyreach historyshape historymissing \
            unsupplied negdegree \
            pastcalculus hugedegree statepath orderzero \
            unfrozen eagerclock flatcalculus shallowcalculus; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

# Constructor-bypass audit, per concrete type: for every
# `type, extends(operation|discretization) :: NAME` (abstract types
# excepted) some procedure in the same file must both mention
# type(NAME)/class(NAME) and call declare_arguments - a module
# constructor, or an attach entry. An operation built by default
# initialization owns no arguments and is refused by the 'undeclared'
# case above.
root="$here/../.."
missing=$(grep -rlE "extends\((operation|discretization)\)" "$root/src" "$root/test" --include=*.f90 \
  | while read -r f; do
      awk '
        /^[ \t]*type[ \t]*,/ && /extends\((operation|discretization)\)[ \t]*::/ && !/abstract/ {
          n = split($0, a, "::"); t = a[2]; sub(/^[ \t]+/, "", t); sub(/[ \t!].*$/, "", t); types[t] = 1 }
        /^[ \t]*(pure[ \t]+|impure[ \t]+|recursive[ \t]+)*(type\([a-z_0-9]+\)[ \t]+|integer[ \t]+|logical[ \t]+)?(function|subroutine)[ \t]/ { inproc = 1; body = "" }
        inproc { body = body "\n" $0 }
        /^[ \t]*end[ \t]+(function|subroutine)/ {
          if (inproc && body ~ /declare_arguments\(/)
            for (t in types) if (body ~ ("(type|class)\\(" t "\\)")) declared[t] = 1
          inproc = 0 }
        END { for (t in types) if (!(t in declared)) print FILENAME ": " t }
      ' "$f"
    done)
if [ -n "$missing" ]; then
    echo " FAIL : concrete operations without a declaring constructor:"; echo "$missing"; exit 1
else
    echo " PASS : every concrete operation type has a procedure that declares its arguments"
fi

# Regression guard: the test sources must not reference the
# deleted gti_* modules.
grep -q "gti_" "$here"/*.f90 \
    && { echo " FAIL : a gti_ reference appears in the suite sources"; exit 1; } \
    || echo " PASS : no gti_ reference in the suite sources"
