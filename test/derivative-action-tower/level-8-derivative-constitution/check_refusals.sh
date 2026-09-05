#!/bin/bash
# Level 8 refusals: each case must die, and die for its stated
# reason - the message is checked, not merely the nonzero exit.
set -e

here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"

declare -A reason=(
  [unbound-law]="no law binds this operation symbol"
  [missing-port]="a port names exactly one slot"
  [primal-starvation]="an operation was scheduled before its primal inputs exist"
  [tangent-starvation]="a tangent was demanded before its input tangents exist"
)

for case in unbound-law missing-port primal-starvation tangent-starvation; do
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
