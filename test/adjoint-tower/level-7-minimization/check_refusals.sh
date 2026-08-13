#!/bin/bash
# Level 7's refusals: every offered field has the RIGHT SIZE and the
# WRONG IDENTITY. Each case must die, and die for its stated reason -
# the message is checked, not merely the nonzero exit.
set -e

here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"

declare -A reason=(
  [primal-rhs-on-Q]="a right-hand side lives on the residual domain"
  [adjoint-rhs-on-Y]="a right-hand side lives on the residual domain"
  [primal-state-on-Y]="the state must live on the state domain"
  [adjoint-covector-on-Q]="the covector must live on the residual-row domain"
)

for case in primal-rhs-on-Q adjoint-rhs-on-Y primal-state-on-Y adjoint-covector-on-Q; do
    if ./refusal "$case" >refusal.out 2>&1; then
        echo " FAIL : '$case' was accepted"
        exit 1
    fi
    if grep -q "${reason[$case]}" refusal.out; then
        echo " PASS : '$case' is refused, loudly"
    else
        echo " FAIL : '$case' died for the wrong reason"
        cat refusal.out
        exit 1
    fi
done
rm -f refusal.out
