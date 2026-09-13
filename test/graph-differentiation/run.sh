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
  [dupslot]="a duplicate argument path is rejected"
  [badslot]="argument() was called with k outside 1..declared_arguments"
  [foreignpath]="a path must name one of its declared arguments"
  [foreignvariation]="found a variation naming an argument of another operation"
  [undeclared]="was called before declare_arguments()"
  [negdegree]="does not support a negative degree"
  [pastcalculus]="the statement does not support the requested order"
  [hugedegree]="the partition coefficient is not representable"
  [unfrozen]="was requested before freeze_inputs was called"
  [flatcalculus]="partial_action() is not implemented by this operation type"
)
for case in dupslot badslot foreignpath foreignvariation undeclared \
            negdegree pastcalculus hugedegree unfrozen flatcalculus; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused" || { cat refusal.out; exit 1; }
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
