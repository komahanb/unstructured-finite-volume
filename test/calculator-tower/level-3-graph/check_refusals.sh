#!/bin/bash
# Level 3 refusals: each case must die, and die for its stated
# reason - the message is checked, not merely the nonzero exit.
set -e

here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"

# A foreign carrier, a domain seated twice and a relation seated twice
# are now answered .false. by relational_valid, in test.f90: they are
# properties of a representation already built. What remains a refusal
# is the storage law.
declare -A reason=(
  [view]="a view cannot be bound"
)

for case in view; do
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
