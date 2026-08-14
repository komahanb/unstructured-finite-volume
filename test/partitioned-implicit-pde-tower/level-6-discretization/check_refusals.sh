#!/bin/bash
# Gate B's refusal: a state of the RIGHT SIZE on the WRONG carrier
# must die, and die for its stated reason - the message is checked,
# not merely the nonzero exit.
set -e

here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"

declare -A reason=(
  [foreign-carrier]="the state must live on this graph's vertex carrier"
)

for case in foreign-carrier; do
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
